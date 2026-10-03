#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ozguler_prop1_check.py
======================

Independent numerical check accompanying:

    Peter Kahl, "Blind by Design: Quantum Error Correction, Functional
    Incompatibility and the Value of Evidence" (2026).

Author
------
Peter Kahl
Independent Researcher, Lex et Ratio
https://www.lexetratio.com
ORCID: 0009-0003-1616-4843

Purpose
-------
Reconstruct the periodic toric-code instrument from Özgüler (2026),
Proposition 1, and numerically test the sign-symmetry and parity-form claims
used in the paper. The implementation follows the cited construction rather
than importing code from the source under examination.

The script checks:

1. F_s(-theta) = F_s(theta) for every syndrome effect;
2. reality and parity form F_s = q_s I + b_s P (for odd L);
3. equality of finite syndrome-history laws at +theta and -theta for complex
   logical inputs;
4. K_s(-theta) = conjugate(K_s(theta)); and
5. persistence of the sign symmetry for unequal per-edge rotations reversed
   globally.

State-vector simulation scales exponentially. The default checks therefore use
L=2 (8 physical qubits) and L=3 (18 physical qubits).

Requirements
------------
Python 3.10+ and NumPy.

Run
---
    python ozguler_prop1_check.py

Licence
-------
MIT License. See LICENSE in the repository root.
SPDX-License-Identifier: MIT

Copyright (c) 2026 Peter Kahl.
"""

from __future__ import annotations

import itertools
import platform

import numpy as np


def build_toric_code(L: int):
    """Build binary X- and Z-check matrices for an L x L periodic toric code.

    Returns the number of edge qubits, the Z-check matrix, the X-check matrix,
    and helper functions mapping horizontal/vertical lattice edges to qubits.
    """
    n = 2 * L * L
    horizontal = lambda x, y: (x % L) + L * (y % L)
    vertical = lambda x, y: L * L + (x % L) + L * (y % L)
    hz = np.zeros((L * L, n), dtype=np.uint8)
    hx = np.zeros((L * L, n), dtype=np.uint8)

    for x in range(L):
        for y in range(L):
            row = x + L * y
            for edge in (
                horizontal(x, y),
                horizontal(x, y + 1),
                vertical(x, y),
                vertical(x + 1, y),
            ):
                hz[row, edge] ^= 1
            for edge in (
                horizontal(x, y),
                horizontal(x - 1, y),
                vertical(x, y),
                vertical(x, y - 1),
            ):
                hx[row, edge] ^= 1

    return n, hz, hx, horizontal, vertical


def binary_rowspace(matrix: np.ndarray) -> list[int]:
    """Return the GF(2) row space as integer bit masks."""
    vectors = {0}
    for row in matrix:
        mask = int("".join(str(bit) for bit in row[::-1]), 2)
        vectors |= {vector ^ mask for vector in vectors}
    return sorted(vectors)


def instrument(L: int, thetas):
    """Construct logical syndrome instruments for specified X-rotation angles.

    ``thetas`` may contain scalar uniform angles or an n-component array of
    per-edge angles. The syndrome basis is reduced by one dependent Z check,
    and each syndrome is assigned a minimum-weight X recovery with a
    lexicographic tie rule.
    """
    n, hz, hx, horizontal, vertical = build_toric_code(L)
    dimension = 1 << n
    hz_reduced = hz[:-1]
    basis_indices = np.arange(dimension, dtype=np.int64)

    syndrome = np.zeros(dimension, dtype=np.int64)
    for k, row in enumerate(hz_reduced):
        mask = int(sum(1 << i for i in np.nonzero(row)[0]))
        parity = np.zeros(dimension, dtype=np.int64)
        masked = basis_indices & mask
        while masked.any():
            parity ^= masked & 1
            masked >>= 1
        syndrome |= parity << k

    recovery: dict[int, int] = {}
    for weight in range(n + 1):
        for combination in itertools.combinations(range(n), weight):
            error = sum(1 << i for i in combination)
            syn = int(syndrome[error])
            recovery.setdefault(syn, error)
        if len(recovery) == 1 << (L * L - 1):
            break

    rowspace = binary_rowspace(hx)
    logical_zero = np.zeros(dimension, dtype=complex)
    logical_zero[rowspace] = 1 / np.sqrt(len(rowspace))

    logical_x = sum(1 << vertical(x, 0) for x in range(L))
    logical_y = sum(1 << horizontal(0, y) for y in range(L))

    codewords = []
    for a in (0, 1):
        for b in (0, 1):
            mask = (logical_x if a else 0) ^ (logical_y if b else 0)
            codewords.append(logical_zero[basis_indices ^ mask])
    code_isometry = np.array(codewords)

    output = {}
    for theta in thetas:
        angle_vector = np.full(n, theta) if np.isscalar(theta) else np.asarray(theta)
        evolved = code_isometry.astype(complex).copy()
        for qubit in range(n):
            cosine = np.cos(angle_vector[qubit] / 2)
            sine = -1j * np.sin(angle_vector[qubit] / 2)
            evolved = cosine * evolved + sine * evolved[:, basis_indices ^ (1 << qubit)]

        kraus = {}
        for syn, rec in recovery.items():
            projected = evolved * (syndrome == syn)
            corrected = projected[:, basis_indices ^ rec]
            kraus[syn] = code_isometry.conj() @ corrected.T

        key = theta if np.isscalar(theta) else "vec"
        output[key] = kraus

    return output, logical_x, logical_y


def check_uniform_angles(L: int, theta: float = 0.1) -> None:
    """Run effect, completeness, parity-form and history checks for one L."""
    logical_x = np.array([[0, 1], [1, 0]])
    parity_operator = np.kron(logical_x, logical_x)

    results, _, _ = instrument(L, [theta, -theta])
    plus, minus = results[theta], results[-theta]

    effect_difference = max(
        np.abs(plus[s].conj().T @ plus[s] - minus[s].conj().T @ minus[s]).max()
        for s in plus
    )
    imaginary_effect = max(
        np.abs((plus[s].conj().T @ plus[s]).imag).max() for s in plus
    )

    parity_form_difference = 0.0
    for syndrome in plus:
        effect = plus[syndrome].conj().T @ plus[syndrome]
        q = np.trace(effect).real / 4
        b = np.trace(effect @ parity_operator).real / 4
        parity_form_difference = max(
            parity_form_difference,
            np.abs(effect - (q * np.eye(4) + b * parity_operator)).max(),
        )

    completeness = np.abs(
        sum(plus[s].conj().T @ plus[s] for s in plus) - np.eye(4)
    ).max()

    rng = np.random.default_rng(1)
    worst_history_difference = 0.0
    rounds = 2 if L == 3 else 3
    for _ in range(3):
        psi = rng.normal(size=4) + 1j * rng.normal(size=4)
        psi /= np.linalg.norm(psi)
        rho = np.outer(psi, psi.conj())
        for history in itertools.product(list(plus), repeat=rounds):
            a = np.eye(4, dtype=complex)
            b = np.eye(4, dtype=complex)
            for syndrome in history:
                a = plus[syndrome] @ a
                b = minus[syndrome] @ b
            p_plus = np.trace(a @ rho @ a.conj().T)
            p_minus = np.trace(b @ rho @ b.conj().T)
            worst_history_difference = max(
                worst_history_difference, abs(p_plus - p_minus)
            )

    print(
        f"L={L}: {len(plus)} syndromes | "
        f"max|F(+)-F(-)| {effect_difference:.1e} | "
        f"max|Im F| {imaginary_effect:.1e} | "
        f"max|F-(qI+bP)| {parity_form_difference:.1e} | "
        f"completeness {completeness:.1e} | "
        f"{rounds}-round history max diff {worst_history_difference:.1e}"
    )


def check_unequal_angles() -> None:
    """Test global sign reversal for unequal rotations on all 18 L=3 edges."""
    rng = np.random.default_rng(7)
    thetas = rng.uniform(0.02, 0.3, 18)
    plus = instrument(3, [thetas])[0]["vec"]
    minus = instrument(3, [-thetas])[0]["vec"]

    effect_difference = max(
        np.abs(plus[s].conj().T @ plus[s] - minus[s].conj().T @ minus[s]).max()
        for s in plus
    )
    conjugacy_difference = max(
        np.abs(minus[s] - plus[s].conj()).max() for s in plus
    )
    print(
        "L=3, unequal angles, global sign flip: "
        f"max|F(+)-F(-)| {effect_difference:.1e} | "
        f"max|K(-)-conj K(+)| {conjugacy_difference:.1e}"
    )


def main() -> None:
    """Run all independent checks of the cited proposition."""
    print(f"Python {platform.python_version()} | numpy {np.__version__}")
    for L in (2, 3):
        check_uniform_angles(L)
    check_unequal_angles()


if __name__ == "__main__":
    main()
