#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
blind_by_design_toy_repetition_code.py
======================================

Version:
    0.1.0 (2026-10-04)

Supplementary research code for:

    Peter Kahl, 'Blind by Design: Quantum Error Correction, Functional
    Incompatibility and the Value of Evidence' (2026), §§6.2-6.3.

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
Construct the three-qubit repetition-code instrument under coherent X
rotations and check the numerical claims of §6.2 of the paper:

1. the syndrome probabilities are the same for every encoded logical input
   (blindness to the stored state) and the same at +theta and -theta
   (blindness to the sign of the rotation);
2. each syndrome effect F_s = K_s^dagger K_s, restricted to the code space, is
   a multiple of the identity, and the effects sum to the identity;
3. the sign matters: a syndrome-conditioned correction tuned to +theta
   restores the memory at +theta and damages it at -theta (entanglement
   fidelity after 50 rounds);
4. an outside 'sentinel' qubit measured in the Y basis reveals the sign
   without touching the protected memory; and
5. every four-round syndrome history has the same probability at +theta and
   -theta, for complex logical inputs and unequal angles on the three qubits
   reversed together.

The numerical agreement in check 5 is a check on the implementation. The
equality itself follows from a real-structure symmetry of the instrument
(§6.2 of the paper): the code states, parity projectors and corrections are
real in the computational basis, and reversing the rotation's sign is complex
conjugation.

Requirements
------------
Python 3.8 or later; NumPy 1.17 or later.

Run
---
    python blind_by_design_toy_repetition_code.py

The script is deterministic and needs no input files.
"""

from __future__ import annotations

import itertools
import platform
from typing import Dict, Tuple

import numpy as np

SCRIPT_VERSION = "0.1.0"

THETA = 0.2              # rotation angle used in §6.2
FIDELITY_ROUNDS = 50     # rounds for the entanglement-fidelity comparison
HISTORY_ROUNDS = 4       # length of the enumerated syndrome histories
UNEQUAL_THETAS = (0.13, 0.27, 0.05)
TOLERANCE = 1e-12        # numerical tolerance for the pass/fail summary

X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.diag([1, -1]).astype(complex)
I2 = np.eye(2, dtype=complex)

Syndrome = Tuple[int, int]


def kron(*operators: np.ndarray) -> np.ndarray:
    """Kronecker product of a sequence of operators."""
    result = np.array([[1]], dtype=complex)
    for operator in operators:
        result = np.kron(result, operator)
    return result


def rx(theta: float) -> np.ndarray:
    """Single-qubit rotation R_x(theta) = exp(-i theta X / 2)."""
    return np.cos(theta / 2) * I2 - 1j * np.sin(theta / 2) * X


# Encoding isometry: |0_L> = |000>, |1_L> = |111>.
V = np.zeros((8, 2), dtype=complex)
V[0, 0] = 1
V[7, 1] = 1

Z1Z2 = kron(Z, Z, I2)
Z2Z3 = kron(I2, Z, Z)
IDENTITY8 = np.eye(8, dtype=complex)

# Syndrome (s1, s2): s1 = 1 when qubits 1 and 2 disagree, s2 = 1 when qubits 2
# and 3 disagree.
PROJECTORS: Dict[Syndrome, np.ndarray] = {
    (s1, s2): ((IDENTITY8 + (-1) ** s1 * Z1Z2) / 2) @ ((IDENTITY8 + (-1) ** s2 * Z2Z3) / 2)
    for s1, s2 in itertools.product((0, 1), repeat=2)
}

# Majority-vote recovery for each syndrome.
RECOVERY: Dict[Syndrome, np.ndarray] = {
    (0, 0): kron(I2, I2, I2),
    (1, 0): kron(X, I2, I2),
    (1, 1): kron(I2, X, I2),
    (0, 1): kron(I2, I2, X),
}

LOGICAL_INPUTS = {
    "0L": np.array([1, 0], dtype=complex),
    "1L": np.array([0, 1], dtype=complex),
    "+L": np.array([1, 1], dtype=complex) / np.sqrt(2),
    "+iL": np.array([1, 1j], dtype=complex) / np.sqrt(2),
}


def logical_kraus(thetas) -> Dict[Syndrome, np.ndarray]:
    """Logical Kraus operators K_s = V^dagger R_s P_s U V for one round.

    ``thetas`` is a scalar (uniform rotation) or three per-qubit angles.
    Every syndrome sector is mapped back into the code space by its recovery,
    so V^dagger loses no probability.
    """
    angles = np.broadcast_to(np.asarray(thetas, dtype=float), (3,))
    unitary = kron(*(rx(float(a)) for a in angles))
    return {
        s: V.conj().T @ RECOVERY[s] @ PROJECTORS[s] @ unitary @ V
        for s in PROJECTORS
    }


def effects(kraus: Dict[Syndrome, np.ndarray]) -> Dict[Syndrome, np.ndarray]:
    """Syndrome effects F_s = K_s^dagger K_s on the code space."""
    return {s: k.conj().T @ k for s, k in kraus.items()}


def polar_corrections(theta: float) -> Dict[Syndrome, np.ndarray]:
    """Syndrome-conditioned logical corrections tuned to +theta.

    Each correction is the inverse of the unitary factor in the polar
    decomposition K_s = W_s |K_s|.
    """
    table = {}
    for s, k in logical_kraus(theta).items():
        u, _, vh = np.linalg.svd(k)
        table[s] = (u @ vh).conj().T
    return table


def entanglement_fidelity(theta: float, corrections, rounds: int) -> float:
    """Entanglement fidelity of the corrected memory after ``rounds`` rounds."""
    kraus = logical_kraus(theta)
    phi = np.zeros(4, dtype=complex)
    phi[0] = phi[3] = 1 / np.sqrt(2)
    state = np.outer(phi, phi.conj())
    for _ in range(rounds):
        state = sum(
            np.kron(corrections[s] @ kraus[s], I2)
            @ state
            @ np.kron(corrections[s] @ kraus[s], I2).conj().T
            for s in kraus
        )
    return float(np.real(phi.conj() @ state @ phi))


def max_history_difference(plus, minus, psi: np.ndarray, rounds: int) -> float:
    """Largest |P(h | +) - P(h | -)| over all syndrome histories of length ``rounds``."""
    rho = np.outer(psi, psi.conj())
    worst = 0.0
    for history in itertools.product(sorted(plus), repeat=rounds):
        a = np.eye(2, dtype=complex)
        b = np.eye(2, dtype=complex)
        for s in history:
            a = plus[s] @ a
            b = minus[s] @ b
        p_plus = np.real(np.trace(a @ rho @ a.conj().T))
        p_minus = np.real(np.trace(b @ rho @ b.conj().T))
        worst = max(worst, abs(p_plus - p_minus))
    return worst


def main() -> None:
    print("=== Run configuration ===")
    print(f"script version = {SCRIPT_VERSION}")
    print(f"Python version = {platform.python_version()}")
    print(f"NumPy version = {np.__version__}")
    print(f"theta = {THETA}; fidelity rounds = {FIDELITY_ROUNDS}; "
          f"history rounds = {HISTORY_ROUNDS}; unequal angles = {UNEQUAL_THETAS}")

    passed = []

    # 1-2. Syndrome laws and effects -------------------------------------------
    print("\n=== 1. Syndrome probabilities for each logical input ===")
    print("syndrome order: (0,0) (0,1) (1,0) (1,1)")
    laws = {}
    for sign in (+1, -1):
        kraus = logical_kraus(sign * THETA)
        for name, psi in LOGICAL_INPUTS.items():
            rho = np.outer(psi, psi.conj())
            law = np.array([np.real(np.trace(kraus[s] @ rho @ kraus[s].conj().T))
                            for s in sorted(kraus)])
            laws[(sign, name)] = law
            print(f"theta = {sign * THETA:+.1f}  input {name:>3}: "
                  + "  ".join(f"{p:.8f}" for p in law))
    reference = laws[(+1, "0L")]
    spread = max(np.abs(law - reference).max() for law in laws.values())
    print(f"largest difference between any two laws: {spread:.1e}")
    passed.append(("syndrome law independent of input and sign", spread < TOLERANCE))

    print("\n=== 2. Syndrome effects at +theta ===")
    eff = effects(logical_kraus(THETA))
    scalar_gap = 0.0
    for s in sorted(eff):
        f = eff[s]
        c = np.real(np.trace(f)) / 2
        scalar_gap = max(scalar_gap, np.abs(f - c * np.eye(2)).max())
        print(f"syndrome {s}: F_s = {c:.10f} x identity")
    completeness = np.abs(sum(eff.values()) - np.eye(2)).max()
    print(f"largest deviation from a multiple of the identity: {scalar_gap:.1e}")
    print(f"completeness |sum F_s - I|: {completeness:.1e}")
    passed.append(("effects are multiples of the identity", scalar_gap < TOLERANCE))
    passed.append(("effects sum to the identity", completeness < TOLERANCE))

    # 3. Decision relevance ------------------------------------------------------
    print(f"\n=== 3. Entanglement fidelity after {FIDELITY_ROUNDS} rounds ===")
    no_extra = {s: np.eye(2, dtype=complex) for s in PROJECTORS}
    tuned = polar_corrections(THETA)
    for sign in (+1, -1):
        theta = sign * THETA
        print(f"true theta = {theta:+.1f}:  no extra correction {entanglement_fidelity(theta, no_extra, FIDELITY_ROUNDS):.6f}"
              f"   correction tuned to +theta {entanglement_fidelity(theta, tuned, FIDELITY_ROUNDS):.6f}")

    # 4. Sentinel ---------------------------------------------------------------
    print("\n=== 4. Sentinel qubit prepared in |0>, rotated, measured in the Y basis ===")
    ket0 = np.array([1, 0], dtype=complex)
    for sign in (+1, -1):
        psi = rx(sign * THETA) @ ket0
        print(f"true theta = {sign * THETA:+.1f}:  <Y> = {np.real(psi.conj() @ Y @ psi):+.6f}")

    # 5. History-level check ----------------------------------------------------
    print(f"\n=== 5. All {HISTORY_ROUNDS}-round syndrome histories, unequal angles reversed together ===")
    plus = logical_kraus(np.array(UNEQUAL_THETAS))
    minus = logical_kraus(-np.array(UNEQUAL_THETAS))
    worst = 0.0
    for label, psi in (("(|0> + i|1>)/sqrt2", np.array([1, 1j], dtype=complex) / np.sqrt(2)),
                       ("0.6|0> + 0.8i|1>", np.array([0.6, 0.8j], dtype=complex))):
        d = max_history_difference(plus, minus, psi, HISTORY_ROUNDS)
        worst = max(worst, d)
        print(f"input {label}: max |P(+) - P(-)| over {4 ** HISTORY_ROUNDS} histories = {d:.1e}")
    passed.append(("history laws equal at +theta and -theta", worst < TOLERANCE))

    print("\n=== Summary ===")
    for label, ok in passed:
        print(f"{'PASS' if ok else 'FAIL'}  {label}")


if __name__ == "__main__":
    main()
