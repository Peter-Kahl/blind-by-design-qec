#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ozguler_prop1_check.py
======================

Version:
    0.2.0 (2026-10-07)

Supplementary research code for:

    Peter Kahl, 'Blind by Design: Quantum Error Correction, Functional
    Incompatibility and the Value of Evidence' (2026), §6.1, and the
    accompanying reliance note (docs/reliance_note_ozguler2026.md).

Author:
    Peter Kahl
    Independent researcher, Lex et Ratio
    https://www.lexetratio.com
    ORCID: 0009-0003-1616-4843

Repository:
    https://github.com/Peter-Kahl/blind-by-design-qec

Copyright (c) 2026 Peter Kahl. MIT License (see LICENSE).
SPDX-License-Identifier: MIT

Purpose
-------
Independently reconstruct the periodic toric-code instrument of Özgüler (2026,
arXiv:2609.19090), Proposition 1, from the preprint's own definitions, and
check numerically the steps of the proposition on which the paper relies.
No code from the preprint is used.

The instrument: an L x L periodic toric code (2 L^2 edge qubits, two logical
qubits), coherent X rotations on every edge, ideal extraction of the Z-type
(plaquette) syndrome with one dependent check removed, and a deterministic
minimum-weight X recovery with a lexicographic tie rule.

For each lattice size the script checks:

1. F_s(-theta) = F_s(theta) for every syndrome effect F_s = K_s^dagger K_s;
2. the effects are real and of the parity form F_s = q_s I + b_s P, where P is
   the product of the two logical X operators (expected for odd L only);
3. completeness, sum_s F_s = I;
4. equality of all multi-round syndrome-history probabilities at +theta and
   -theta, for random complex logical inputs; and, at L = 3,
5. persistence of sign symmetry, and K_s(-theta) = conj(K_s(theta)), for
   unequal per-edge angles reversed together.

Sanity checks confirm that the constructed code states are orthonormal,
stabilised by every plaquette and star check, and that the logical operators
commute with the plaquette checks.

Scope
-----
Proposition 1 is stated for odd L >= 3. L = 2 is included as a contrast
outside that scope: there the parity form is not expected to hold, although
sign symmetry does. That observation is not the preprint's and is not used
in the paper. State-vector simulation scales as 2^(2 L^2), so only L = 2
(8 qubits) and L = 3 (18 qubits) are tractable. The checks are finite
numerical checks; the general result rests on the proof (see the reliance
note).

Requirements
------------
Python 3.8 or later; NumPy 1.17 or later.

Run
---
    python ozguler_prop1_check.py

Random inputs use fixed seeds; the run is deterministic.
"""

from __future__ import annotations

import itertools
import platform

import numpy as np

SCRIPT_VERSION = "0.2.0"
THETA = 0.1
INPUT_SEED = 1
UNEQUAL_ANGLE_SEED = 7
N_RANDOM_INPUTS = 3


def toric_code(L: int):
    """Plaquette (Z) and star (X) check matrices of an L x L periodic toric code.

    Edges are indexed horizontal(x, y) = x + L y and
    vertical(x, y) = L^2 + x + L y, with periodic coordinates.
    """
    n = 2 * L * L

    def horizontal(x, y):
        return (x % L) + L * (y % L)

    def vertical(x, y):
        return L * L + (x % L) + L * (y % L)

    hz = np.zeros((L * L, n), dtype=np.uint8)
    hx = np.zeros((L * L, n), dtype=np.uint8)
    for x in range(L):
        for y in range(L):
            row = x + L * y
            for e in (horizontal(x, y), horizontal(x, y + 1), vertical(x, y), vertical(x + 1, y)):
                hz[row, e] ^= 1          # plaquette
            for e in (horizontal(x, y), horizontal(x - 1, y), vertical(x, y), vertical(x, y - 1)):
                hx[row, e] ^= 1          # star
    return n, hz, hx, horizontal, vertical


def mask_of(row) -> int:
    return int(sum(1 << int(i) for i in np.nonzero(row)[0]))


def parity_of_masked(indices: np.ndarray, mask: int) -> np.ndarray:
    """Bitwise parity of (indices & mask) for every basis index."""
    masked = indices & mask
    parity = np.zeros_like(indices)
    while masked.any():
        parity ^= masked & 1
        masked >>= 1
    return parity


def build_instrument(L: int):
    """Return the code isometry, syndrome map, recoveries and logical masks."""
    n, hz, hx, horizontal, vertical = toric_code(L)
    dim = 1 << n
    idx = np.arange(dim, dtype=np.int64)

    # Reduced Z syndrome: the last plaquette check is the product of the others.
    syndrome = np.zeros(dim, dtype=np.int64)
    for k, row in enumerate(hz[:-1]):
        syndrome |= parity_of_masked(idx, mask_of(row)) << k
    n_syndromes = 1 << (L * L - 1)

    # Minimum-weight X recovery; ties broken by the lexicographically first
    # support (itertools.combinations order).
    recovery = {}
    for weight in range(n + 1):
        for support in itertools.combinations(range(n), weight):
            e = sum(1 << i for i in support)
            recovery.setdefault(int(syndrome[e]), e)
        if len(recovery) == n_syndromes:
            break

    # |0_L 0_L>: uniform superposition over the X-star group applied to |0...0>.
    star_group = {0}
    for row in hx:
        m = mask_of(row)
        star_group |= {v ^ m for v in star_group}
    zero = np.zeros(dim, dtype=complex)
    zero[sorted(star_group)] = 1 / np.sqrt(len(star_group))

    # Two logical X operators: non-contractible cycles of the dual lattice.
    logical_x1 = sum(1 << vertical(x, 0) for x in range(L))
    logical_x2 = sum(1 << horizontal(0, y) for y in range(L))

    codewords = []
    for a in (0, 1):
        for b in (0, 1):
            m = (logical_x1 if a else 0) ^ (logical_x2 if b else 0)
            codewords.append(zero[idx ^ m])
    C = np.array(codewords)        # rows: |a b>_L, a and b in {0, 1}
    return dict(n=n, hz=hz, hx=hx, idx=idx, syndrome=syndrome, recovery=recovery,
                C=C, logical_x1=logical_x1, logical_x2=logical_x2)


def sanity_checks(inst) -> float:
    """Largest violation of orthonormality, stabiliser and commutation conditions."""
    C, idx = inst["C"], inst["idx"]
    worst = np.abs(C.conj() @ C.T - np.eye(4)).max()
    for row in inst["hz"]:                       # Z plaquettes: eigenvalue +1
        sign = 1 - 2 * parity_of_masked(idx, mask_of(row))
        worst = max(worst, np.abs(C * sign - C).max())
    for row in inst["hx"]:                       # X stars: eigenvalue +1
        worst = max(worst, np.abs(C[:, idx ^ mask_of(row)] - C).max())
    for m in (inst["logical_x1"], inst["logical_x2"]):   # commute with plaquettes
        for row in inst["hz"]:
            worst = max(worst, float(bin(m & mask_of(row)).count("1") % 2))
    return float(worst)


def kraus_operators(inst, angles):
    """Logical Kraus operators K_s[i, j] = <c_i| R_s P_s U |c_j> for one round."""
    C, idx, n = inst["C"], inst["idx"], inst["n"]
    angles = np.broadcast_to(np.asarray(angles, dtype=float), (n,))
    evolved = C.astype(complex).copy()
    for q in range(n):
        evolved = np.cos(angles[q] / 2) * evolved - 1j * np.sin(angles[q] / 2) * evolved[:, idx ^ (1 << q)]
    kraus = {}
    for s, r in inst["recovery"].items():
        projected = evolved * (inst["syndrome"] == s)
        corrected = projected[:, idx ^ r]
        kraus[s] = C.conj() @ corrected.T
    return kraus


def max_history_difference(plus, minus, rounds, rng, n_inputs):
    worst = 0.0
    keys = list(plus)
    for _ in range(n_inputs):
        psi = rng.normal(size=4) + 1j * rng.normal(size=4)
        psi /= np.linalg.norm(psi)
        rho = np.outer(psi, psi.conj())
        for history in itertools.product(keys, repeat=rounds):
            a = np.eye(4, dtype=complex)
            b = np.eye(4, dtype=complex)
            for s in history:
                a = plus[s] @ a
                b = minus[s] @ b
            worst = max(worst, abs(np.trace(a @ rho @ a.conj().T) - np.trace(b @ rho @ b.conj().T)))
    return float(worst)


def check_lattice(L: int) -> None:
    inst = build_instrument(L)
    plus = kraus_operators(inst, THETA)
    minus = kraus_operators(inst, -THETA)
    f_plus = {s: k.conj().T @ k for s, k in plus.items()}
    f_minus = {s: k.conj().T @ k for s, k in minus.items()}
    P = np.kron(np.array([[0, 1], [1, 0]]), np.array([[0, 1], [1, 0]]))

    sign_gap = max(np.abs(f_plus[s] - f_minus[s]).max() for s in f_plus)
    imag = max(np.abs(f.imag).max() for f in f_plus.values())
    parity_gap = 0.0
    for f in f_plus.values():
        q = np.trace(f).real / 4
        b = np.trace(f @ P).real / 4
        parity_gap = max(parity_gap, np.abs(f - (q * np.eye(4) + b * P)).max())
    completeness = np.abs(sum(f_plus.values()) - np.eye(4)).max()
    rounds = 2 if L >= 3 else 3
    history = max_history_difference(plus, minus, rounds, np.random.default_rng(INPUT_SEED), N_RANDOM_INPUTS)

    scope = "within the proposition's scope (odd L)" if L % 2 else "outside the proposition's scope (even L): parity form not expected"
    print(f"\n=== L = {L}: {inst['n']} qubits, {len(plus)} syndromes; {scope} ===")
    print(f"code-space sanity checks, largest violation   {sanity_checks(inst):.1e}")
    print(f"max |F_s(+theta) - F_s(-theta)|                {sign_gap:.1e}")
    print(f"max |Im F_s|                                   {imag:.1e}")
    print(f"max |F_s - (q_s I + b_s P)|                    {parity_gap:.1e}")
    print(f"completeness |sum F_s - I|                     {completeness:.1e}")
    print(f"max history difference ({rounds} rounds, {N_RANDOM_INPUTS} random inputs) {history:.1e}")


def check_unequal_angles(L: int = 3) -> None:
    inst = build_instrument(L)
    angles = np.random.default_rng(UNEQUAL_ANGLE_SEED).uniform(0.02, 0.3, inst["n"])
    plus = kraus_operators(inst, angles)
    minus = kraus_operators(inst, -angles)
    sign_gap = max(np.abs(plus[s].conj().T @ plus[s] - minus[s].conj().T @ minus[s]).max() for s in plus)
    conj_gap = max(np.abs(minus[s] - plus[s].conj()).max() for s in plus)
    print(f"\n=== L = {L}: unequal angles on all {inst['n']} edges, reversed together ===")
    print(f"max |F_s(+) - F_s(-)|                          {sign_gap:.1e}")
    print(f"max |K_s(-) - conj K_s(+)|                     {conj_gap:.1e}")


def main() -> None:
    print("=== Run configuration ===")
    print(f"script version = {SCRIPT_VERSION}")
    print(f"Python version = {platform.python_version()}")
    print(f"NumPy version = {np.__version__}")
    print(f"theta = {THETA}; input seed = {INPUT_SEED}; unequal-angle seed = {UNEQUAL_ANGLE_SEED}")
    for L in (2, 3):
        check_lattice(L)
    check_unequal_angles(3)


if __name__ == "__main__":
    main()
