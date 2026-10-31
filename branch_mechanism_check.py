#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
branch_mechanism_check.py
=========================

Version:
    0.9.0 (2026-10-31)

Supplementary research code for:

    Peter Kahl, 'The Cost of Sign-Blind Recovery under Coherent Noise: Worked Numerical Cases'
    (2026), §4 (mechanism) and §3 (branch undo in diamond norm).

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
infidelity_and_calibration_check.py finds that, by entanglement infidelity,
knowing the sign of the rotation barely helps without the syndrome history,
the history barely helps without the sign, and the two together remove most
of the one-round infidelity. This script shows why, branch by branch, for
the L = 3 toric instrument of Özgüler (2026).

For each syndrome s, K_s is the Kraus operator of that branch on the code
space (dimension 4): rotation, syndrome projection and minimum-weight
correction. The script checks and reports:

1. whether K_s is proportional to a unitary (relative spread of its
   singular values);
2. whether the branch at -theta is the complex conjugate of the branch at
   +theta (so the two signs rotate the logical state in opposite
   directions);
3. the logical rotation angle of the branch, taken as the spread of the
   eigenphases of its nearest unitary (for exp(-i phi P / 2) with P a
   Pauli operator, the spread is phi); and
4. the branch's exact contribution to one-round entanglement infidelity.

The decomposition is exact: for a channel with Kraus operators K_s on
dimension d, entanglement infidelity is sum_s (p_s - |tr K_s|^2 / d^2), with
p_s = tr(K_s^dagger K_s) / d the branch probability. Undoing branch s with
the inverse of its nearest unitary leaves (tr of the singular values)^2 / d^2
in place of |tr K_s|^2 / d^2, which equals p_s when K_s is proportional to
a unitary.

Branches are grouped into the syndrome classes of prop52_numerical_check.py
(syndromes whose branch maps agree up to weight); members of a class share
probability, angle and contribution.

Three repairs are then compared by one-round infidelity: no repair;
undoing each branch's rotation with the sign known (a feasible
sign-known, history-used policy, 'branch undo'); and the same repairs applied
at the wrong sign.

Branch undo is the correction Özgüler (2026) constructs by polar
decomposition of each branch, his 'inverse-polar table' (Supplemental
§IV.A). As a cross-check, the script reproduces his three 300-round figures
at |theta| = 0.1: no correction, the correctly signed table, and the table
applied at the wrong sign (his §II).

Finally, the diamond-norm defect of branch undo is computed for one round at
theta = 0.05, 0.1, 0.2 and 0.3, and for two rounds (one repair after both,
chosen from the two-round history) at theta = 0.05 and 0.1. Each defect is
reported as [entangled-input lower bound, upper bound], with the upper bound
from diamond_norm_of_difference. Branch undo is a feasible policy, so its
defect is an upper bound on the sign-known, history-used optimum; at these
angles it is tighter than the solver's optimum in prop52_numerical_check.py
and two_round_repair_check.py.

Run
---
    python branch_mechanism_check.py

Imports prop52_numerical_check.py and ozguler_prop1_check.py from the same
folder (and two_round_repair_check.py, for the two-round classes). About
a minute.
"""

from __future__ import annotations

import platform

import numpy as np

from ozguler_prop1_check import build_instrument, kraus_operators
from prop52_numerical_check import (choi, diamond_norm_of_difference,
                                    entangled_input_lower_bound,
                                    identity_choi, syndrome_classes)
from two_round_repair_check import two_round_classes

SCRIPT_VERSION = "0.9.0"

THETAS = (0.05, 0.1)
DIAMOND_THETAS = (0.05, 0.1, 0.2, 0.3)
TWO_ROUND_THETAS = (0.05, 0.1)
DIM = 4
TOLERANCE = 1e-12
ROWS_SHOWN = 6


def nearest_unitary(K):
    u, _, vh = np.linalg.svd(K)
    return u @ vh


def rotation_angle(W):
    """Spread of the eigenphases of W after removing a global phase."""
    w = np.linalg.eigvals(W)
    phases = np.angle(w * np.conj(w[0]))
    return float(np.ptp(phases))


def branch_data(K_plus, K_minus):
    p = float(np.real(np.trace(K_plus.conj().T @ K_plus))) / DIM
    sv = np.linalg.svd(K_plus, compute_uv=False)
    spread = float((sv.max() - sv.min()) / sv.max())
    conj_error = float(np.abs(K_minus - K_plus.conj()).max())
    W = nearest_unitary(K_plus)
    angle = rotation_angle(W)
    contribution = p - abs(np.trace(K_plus)) ** 2 / DIM ** 2
    undone = p - sv.sum() ** 2 / DIM ** 2
    wrong = p - abs(np.trace(W.conj().T @ K_minus)) ** 2 / DIM ** 2
    return p, spread, conj_error, angle, contribution, undone, wrong


def report(theta, inst, passed):
    kp, km = kraus_operators(inst, theta), kraus_operators(inst, -theta)
    classes = syndrome_classes(kp, km)
    rows = []
    totals = np.zeros(3)
    worst_conj, worst_spread_rare = 0.0, 0.0
    for c in classes:
        data = [branch_data(kp[s], km[s]) for s in c]
        p = sum(d[0] for d in data)
        contribution = sum(d[4] for d in data)
        undone = sum(d[5] for d in data)
        wrong = sum(d[6] for d in data)
        totals += (contribution, undone, wrong)
        worst_conj = max(worst_conj, max(d[2] for d in data))
        rows.append((p, len(c), data[0][3], max(d[1] for d in data), contribution))
    rows.sort(reverse=True)
    infidelity = totals[0]

    print(f"\n=== theta = {theta:.2f} ===")
    print(f" {len(classes)} syndrome classes; the {ROWS_SHOWN} most probable:")
    print(" class prob.   syndromes  rotation angle   singular-value spread   share of infidelity")
    for p, n, angle, spread, contribution in rows[:ROWS_SHOWN]:
        print(f" {p:.3e}      {n:3d}      {angle:.3e}        {spread:.1e}                 "
              f"{100 * contribution / infidelity:5.1f} per cent")
    trivial = rows[0]
    print(f" most probable class: probability {trivial[0]:.4f}, rotation angle {trivial[2]:.2e} "
          f"(theta^3 = {theta ** 3:.2e}), share of infidelity {100 * trivial[4] / infidelity:.1f} per cent")
    rare = [r for r in rows[1:] if r[0] > 0]
    angles = [r[2] for r in rare[:ROWS_SHOWN]]
    print(f" next classes: rotation angles {min(angles):.4f} to {max(angles):.4f} (theta = {theta})")
    print(f" largest |K(-theta) - conj K(+theta)| over all syndromes: {worst_conj:.1e}")
    print(" One-round entanglement infidelity:")
    print(f"  no repair                                         {infidelity:.4e}")
    print(f"  each branch's rotation undone, sign known         {totals[1]:.4e}   "
          f"({100 * (1 - totals[1] / infidelity):.1f} per cent removed)")
    print(f"  the same repairs applied at the wrong sign        {totals[2]:.4e}   "
          f"({totals[2] / infidelity:.1f} times no repair)")

    passed.append((f"theta = {theta:.2f}: branch at -theta is the complex conjugate of the branch at +theta",
                   worst_conj < TOLERANCE))
    passed.append((f"theta = {theta:.2f}: undoing each branch with the sign known removes more than 90 per cent "
                   "of one-round infidelity", totals[1] < 0.1 * infidelity))
    passed.append((f"theta = {theta:.2f}: the same repairs at the wrong sign are worse than no repair",
                   totals[2] > infidelity))


def undo(K):
    u, _, vh = np.linalg.svd(K)
    return (u @ vh).conj().T @ K


def diamond_of(kraus_list):
    D = choi(kraus_list, DIM) - identity_choi(DIM)
    return entangled_input_lower_bound(D, DIM), diamond_norm_of_difference(D, DIM)


def report_diamond(inst, passed):
    print("\n=== Diamond-norm defect of branch undo, sign known ===")
    print(" [entangled-input lower bound, upper bound]")
    print(" theta  one round                    two rounds (one repair after both)")
    ok = True
    for theta in DIAMOND_THETAS:
        kp = kraus_operators(inst, theta)
        one = diamond_of([undo(K) for K in kp.values()])
        line = f" {theta:4.2f}   [{one[0]:.3e}, {one[1]:.3e}]"
        if theta in TWO_ROUND_THETAS:
            _, reps, residual, covered = two_round_classes(inst, theta)
            two = diamond_of([undo(a) for a, _ in reps])
            line += f"     [{two[0]:.3e}, {two[1]:.3e}]   ({len(reps)} classes)"
            ok &= two[0] <= two[1] * (1 + 1e-9) and covered == 256 ** 2
        ok &= one[0] <= one[1] * (1 + 1e-9)
        print(line)
    passed.append(("branch undo: diamond bounds ordered (lower <= upper)", ok))


OZGULER_300 = {"no correction": 0.1323553, "correctly signed table": 0.0007456,
               "table at the wrong sign": 0.4140468}   # Özgüler (2026), §II


def superop(kraus_list):
    return sum(np.kron(K, K.conj()) for K in kraus_list)


def infidelity_after(S, rounds):
    J = np.zeros((DIM * DIM, DIM * DIM), dtype=complex)
    SN = np.linalg.matrix_power(S, rounds)
    for i in range(DIM):
        for j in range(DIM):
            E = np.zeros((DIM, DIM), dtype=complex)
            E[i, j] = 1
            J[i * DIM:(i + 1) * DIM, j * DIM:(j + 1) * DIM] = (SN @ E.reshape(-1)).reshape(DIM, DIM)
    return 1 - float(np.real(np.trace(identity_choi(DIM) @ J))) / DIM ** 2


def report_ozguler(inst, passed):
    print("\n=== Cross-check: Özgüler (2026, §II), |theta| = 0.10, 300 rounds ===")
    kp, km = kraus_operators(inst, 0.1), kraus_operators(inst, -0.1)
    tables = {s: nearest_unitary(K).conj().T for s, K in kp.items()}
    ours = {
        "no correction": infidelity_after(superop(list(kp.values())), 300),
        "correctly signed table": infidelity_after(superop([tables[s] @ K for s, K in kp.items()]), 300),
        "table at the wrong sign": infidelity_after(superop([tables[s] @ K for s, K in km.items()]), 300),
    }
    ok = True
    for name, value in ours.items():
        print(f" {name:<25} {value:.7f}   (Özgüler: {OZGULER_300[name]:.7f})")
        ok &= round(value, 7) == OZGULER_300[name]
    passed.append(("cross-check: Özgüler's three 300-round infidelities reproduced to seven decimal places", ok))


def main() -> None:
    print("=== Run configuration ===")
    print(f"script version = {SCRIPT_VERSION}")
    print(f"Python version = {platform.python_version()}")
    print(f"NumPy version = {np.__version__}")
    print(f"thetas = {THETAS}")
    print("infidelity = 1 - entanglement fidelity; decomposition sum_s (p_s - |tr K_s|^2 / 16) is exact")

    passed = []
    inst = build_instrument(3)
    for theta in THETAS:
        report(theta, inst, passed)
    report_ozguler(inst, passed)
    report_diamond(inst, passed)

    print("\n=== Summary ===")
    for label, ok in passed:
        print(f"{'PASS' if ok else 'FAIL'}  {label}")


if __name__ == "__main__":
    main()
