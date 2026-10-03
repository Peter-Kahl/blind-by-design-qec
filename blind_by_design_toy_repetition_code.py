#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
blind_by_design_toy_repetition_code.py
======================================

Supplementary numerical experiment for:

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
Construct a three-qubit repetition-code instrument under coherent X rotations
and verify the sign-blindness claims used in the paper. The script checks:

1. syndrome probabilities for +theta and -theta;
2. whether syndrome effects are scalar on the logical subspace;
3. logical entanglement fidelity under an incumbent recovery and a recovery
   table tuned at +theta;
4. a simple external sentinel measurement that distinguishes the sign; and
5. equality of all four-round syndrome-history probabilities under a global
   sign reversal of unequal physical rotation angles.

Requirements
------------
Python 3.10+ and NumPy.

Run
---
    python blind_by_design_toy_repetition_code.py

The script is deterministic and requires no input files.

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


X = np.array([[0, 1], [1, 0]], dtype=complex)
I = np.eye(2, dtype=complex)
Z = np.diag([1, -1]).astype(complex)

N_PHYSICAL_QUBITS = 3
HISTORY_ROUNDS = 4
THETA = 0.2


def kron(*operators: np.ndarray) -> np.ndarray:
    """Return the Kronecker product of an arbitrary sequence of operators."""
    result = np.array([[1]], dtype=complex)
    for operator in operators:
        result = np.kron(result, operator)
    return result


def rx(theta: float) -> np.ndarray:
    """Return the single-qubit X-axis rotation R_x(theta)."""
    return np.cos(theta / 2) * I - 1j * np.sin(theta / 2) * X


# Logical encoding |0_L> = |000>, |1_L> = |111>.
_e0 = np.zeros(2**N_PHYSICAL_QUBITS, dtype=complex)
_e1 = np.zeros(2**N_PHYSICAL_QUBITS, dtype=complex)
_e0[0] = 1
_e1[-1] = 1
V = np.stack([_e0, _e1], axis=1)

Z1Z2 = kron(Z, Z, I)
Z2Z3 = kron(I, Z, Z)


def build_syndrome_projectors() -> dict[tuple[int, int], np.ndarray]:
    """Construct projectors for the two repetition-code stabiliser syndromes."""
    projectors: dict[tuple[int, int], np.ndarray] = {}
    identity = np.eye(2**N_PHYSICAL_QUBITS, dtype=complex)
    for s1, s2 in itertools.product((0, 1), repeat=2):
        p1 = (identity + (-1) ** s1 * Z1Z2) / 2
        p2 = (identity + (-1) ** s2 * Z2Z3) / 2
        projectors[(s1, s2)] = p1 @ p2
    return projectors


PROJECTORS = build_syndrome_projectors()
RECOVERY = {
    (0, 0): kron(I, I, I),
    (1, 0): kron(X, I, I),
    (1, 1): kron(I, X, I),
    (0, 1): kron(I, I, X),
}

LOGICAL_INPUTS = {
    "0L": np.array([1, 0], dtype=complex),
    "1L": np.array([0, 1], dtype=complex),
    "+L": np.array([1, 1], dtype=complex) / np.sqrt(2),
    "+iL": np.array([1, 1j], dtype=complex) / np.sqrt(2),
}


def kraus(theta: float) -> dict[tuple[int, int], np.ndarray]:
    """Return logical Kraus operators for a uniform physical X rotation."""
    unitary = kron(rx(theta), rx(theta), rx(theta))
    return {
        syndrome: V.conj().T @ RECOVERY[syndrome] @ projector @ unitary @ V
        for syndrome, projector in PROJECTORS.items()
    }


def kraus_unequal(thetas: np.ndarray) -> dict[tuple[int, int], np.ndarray]:
    """Return logical Kraus operators for independent rotations on the three qubits."""
    if len(thetas) != N_PHYSICAL_QUBITS:
        raise ValueError(f"Expected {N_PHYSICAL_QUBITS} rotation angles.")
    unitary = kron(*(rx(float(theta)) for theta in thetas))
    return {
        syndrome: V.conj().T @ RECOVERY[syndrome] @ projector @ unitary @ V
        for syndrome, projector in PROJECTORS.items()
    }


def polar_tables(theta: float) -> dict[tuple[int, int], np.ndarray]:
    """Construct syndrome-conditioned logical corrections from polar unitaries."""
    table: dict[tuple[int, int], np.ndarray] = {}
    for syndrome, operator in kraus(theta).items():
        u, _singular_values, vh = np.linalg.svd(operator)
        table[syndrome] = (u @ vh).conj().T
    return table


def entanglement_fidelity(
    theta: float,
    correction_table: dict[tuple[int, int], np.ndarray],
    rounds: int = 50,
) -> float:
    """Return logical entanglement fidelity after repeated corrected rounds."""
    operators = kraus(theta)
    phi = np.zeros(4, dtype=complex)
    phi[0] = phi[3] = 1 / np.sqrt(2)
    state = np.outer(phi, phi.conj())

    for _ in range(rounds):
        state = sum(
            np.kron(correction_table[s] @ operators[s], I)
            @ state
            @ np.kron(correction_table[s] @ operators[s], I).conj().T
            for s in operators
        )

    return float(np.real(phi.conj() @ state @ phi))


def print_single_round_checks(theta: float) -> None:
    """Print syndrome laws and logical effects for the basic sign-reversal test."""
    for sign in (1, -1):
        operators = kraus(sign * theta)
        rows = []
        for name, psi in LOGICAL_INPUTS.items():
            rho = np.outer(psi, psi.conj())
            probabilities = [
                np.real(np.trace(operators[s] @ rho @ operators[s].conj().T))
                for s in sorted(operators)
            ]
            rows.append((name, np.round(probabilities, 8)))
        print("theta =", sign * theta, rows)

    print("\nLogical effects F_s = K_s^† K_s at +theta:")
    for syndrome, operator in sorted(kraus(theta).items()):
        effect = operator.conj().T @ operator
        print(syndrome, np.round(effect, 10).tolist())


def print_fidelity_checks(theta: float) -> None:
    """Compare incumbent and +theta-tuned logical correction tables."""
    identity_table = {s: np.eye(2, dtype=complex) for s in kraus(theta)}
    tuned_table = polar_tables(theta)
    for sign in (1, -1):
        true_theta = sign * theta
        print(
            "true theta",
            true_theta,
            "Fe incumbent",
            round(entanglement_fidelity(true_theta, identity_table), 6),
            "Fe with +theta table",
            round(entanglement_fidelity(true_theta, tuned_table), 6),
        )


def print_sentinel_check(theta: float) -> None:
    """Show that an added Y measurement distinguishes +theta from -theta."""
    y_operator = np.array([[0, -1j], [1j, 0]], dtype=complex)
    ket0 = np.array([1, 0], dtype=complex)
    for sign in (1, -1):
        psi = rx(sign * theta) @ ket0
        expectation = np.real(psi.conj() @ y_operator @ psi)
        print("sentinel <Y>", sign * theta, expectation)


def history_level_check() -> None:
    """Verify equal four-round history laws under global sign reversal.

    The check uses unequal physical angles and complex logical inputs, matching
    the history-level numerical check cited in §6.2 of the paper.
    """
    thetas = np.array([0.13, 0.27, 0.05])
    plus = kraus_unequal(thetas)
    minus = kraus_unequal(-thetas)

    inputs = (
        np.array([1, 1j], dtype=complex) / np.sqrt(2),
        np.array([0.6, 0.8j], dtype=complex),
    )

    for psi in inputs:
        rho = np.outer(psi, psi.conj())
        worst = 0.0
        for history in itertools.product(sorted(plus), repeat=HISTORY_ROUNDS):
            rho_plus = rho.copy()
            rho_minus = rho.copy()
            for syndrome in history:
                rho_plus = plus[syndrome] @ rho_plus @ plus[syndrome].conj().T
                rho_minus = minus[syndrome] @ rho_minus @ minus[syndrome].conj().T
            worst = max(worst, abs(np.trace(rho_plus) - np.trace(rho_minus)))
        print(
            "max |P(+theta) - P(-theta)| over all "
            f"{HISTORY_ROUNDS}-round histories:",
            worst,
        )


def main() -> None:
    """Run all supplementary checks and print their numerical results."""
    print(f"Python {platform.python_version()} | numpy {np.__version__}")
    print_single_round_checks(THETA)
    print("\nFidelity checks:")
    print_fidelity_checks(THETA)
    print("\nExternal sentinel check:")
    print_sentinel_check(THETA)
    print("\nHistory-level check:")
    history_level_check()


if __name__ == "__main__":
    main()
