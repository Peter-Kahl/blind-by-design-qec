#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
prop52_numerical_check.py
=========================

Version:
    0.2.0 (2026-10-07)

Supplementary research code for:

    Peter Kahl, 'Blind by Design: Quantum Error Correction, Functional
    Incompatibility and the Value of Evidence' (2026), §5.3.1
    (Proposition 5.2, approximate classical-output blindness).

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
Check numerically the inequality of Proposition 5.2. The proposition states
that if the logical channel Lambda of an inquiry can be repaired to within
epsilon of its prescribed operation in diamond norm, then the complete
classical record distinguishes any two encoded states with total-variation
distance at most 2 sqrt(epsilon); over N rounds, chosen adaptively or not, at
most 2 sqrt(epsilon_1 + ... + epsilon_N). The proof rests on Kretschmann,
Schlingemann and Werner (2008), Theorem 3.

Conventions: the diamond norm is unnormalised, with range [0, 2]; the leak is
a total-variation distance, with range [0, 1].

The script uses two instruments:

1. the L = 3 periodic toric code of Özgüler (2026, Proposition 1) under
   uniform coherent X rotations by theta, rebuilt from the preprint's
   definitions by ozguler_prop1_check.py in this repository; and
2. the three-qubit bit-flip repetition code under coherent X rotations
   (§6.2 of the paper).

For the toric code it computes, at each angle:

* leak_1   the largest total-variation distance between the one-round
           syndrome laws of any two encoded states, computed exactly;
* eps_id   the diamond distance ||Lambda - id|| of the one-round logical
           channel from the identity, by semidefinite programme (SDP);
* eps_dec  the same after a unitary decoder fitted to Lambda (the polar part
           of the leading Choi eigenvector); this bounds from above the
           infimum over decoders that appears in the proposition;
* LB       the entangled-input lower bound on each diamond distance; and
* the bounds 2 sqrt(eps) and their ratio to the leak.

For N = 1, 2 and 3 rounds of the same instrument it compares the exact
N-round leak with 2 sqrt(||Lambda^N - id||) and checks that
||Lambda^N - id|| <= N eps_id, as part (ii) of the proposition requires.

For the repetition code it checks that every three-round history law is
independent of the encoded state, although the logical channel is not the
identity: the converse of the proposition fails.

For contrast, it reports the sign-discrimination profile of a single
sentinel qubit measured in the Y basis (per-shot total variation sin theta),
which the passive syndrome record lacks entirely (Özgüler 2026,
Proposition 1; checked by ozguler_prop1_check.py).

Numerical safeguards
--------------------
Every reported diamond distance is a certified upper bound: the SDP's dual
point is shifted until it is exactly feasible before the objective is
evaluated. Every check of the inequality is made at the entangled-input
lower bound on epsilon. Neither verdict therefore depends on solver
tolerance. The SDP is validated against the closed form 2 sin(phi/2) for the
unitary channel diag(1, e^{i phi}).

Scope
-----
The toric-code computations are limited to L = 3 (18 qubits) by the size of
the state vector, and to three rounds by the number of histories (256^N).
These are finite numerical checks of an inequality whose proof is given in
the paper; they do not establish it. The leak is computed exactly using the
parity form of the history effects, F_h = q_h I + g_h P with P the product of
the two logical X operators; the script reports the largest deviation from
that form, which must be at rounding level for the computed leak to be exact.

Requirements
------------
Python 3.9 or later; NumPy 1.20 or later; CVXPY 1.4 or later with the
Clarabel solver (installed with CVXPY by default).

Run
---
    python prop52_numerical_check.py

The script is deterministic and needs no input files. It imports
ozguler_prop1_check.py from the same folder.
"""

from __future__ import annotations

import platform
import warnings
from math import comb

import cvxpy as cp
import numpy as np

from ozguler_prop1_check import build_instrument, kraus_operators

SCRIPT_VERSION = "0.2.0"

THETAS = (0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5)   # one-round angle grid
ROUNDS_THETAS = (0.1, 0.3)                       # angles for the N-round check
MAX_ROUNDS = 3                                   # N = 1, ..., MAX_ROUNDS
REPETITION_ROUNDS = 3                            # history length, repetition code
SENTINEL_SHOTS = (1, 10, 100)
SOLVER = "CLARABEL"
TOLERANCE = 1e-12        # rounding tolerance for the pass/fail summary

X = np.array([[0, 1], [1, 0]], dtype=complex)
PARITY = np.kron(X, X)   # product of the two logical X operators

SOLVER_FLAGS = []        # status of every SDP solve not reported 'optimal'


# ---------------------------------------------------------------------------
# Channels and the diamond norm
# ---------------------------------------------------------------------------

def choi(kraus_list, dim):
    """Choi matrix sum_ij |i><j| (x) Lambda(|i><j|), input factor first."""
    J = np.zeros((dim * dim, dim * dim), dtype=complex)
    for i in range(dim):
        for j in range(dim):
            E = np.zeros((dim, dim), dtype=complex)
            E[i, j] = 1
            out = sum(K @ E @ K.conj().T for K in kraus_list)
            J[i * dim:(i + 1) * dim, j * dim:(j + 1) * dim] = out
    return J


def choi_from_superoperator(S, dim):
    """Choi matrix of the channel whose action on row-major vec(rho) is S."""
    J = np.zeros((dim * dim, dim * dim), dtype=complex)
    for i in range(dim):
        for j in range(dim):
            E = np.zeros((dim, dim), dtype=complex)
            E[i, j] = 1
            J[i * dim:(i + 1) * dim, j * dim:(j + 1) * dim] = (S @ E.reshape(-1)).reshape(dim, dim)
    return J


def identity_choi(dim):
    omega = np.eye(dim).reshape(-1)
    return np.outer(omega, omega).astype(complex)


def diamond_norm_of_difference(J, dim):
    """Certified upper bound on the unnormalised diamond norm of a difference
    of channels with Choi matrix J.

    Watrous (2013) dual: ||Phi||_diamond = 2 min t subject to Z >= 0,
    Z >= J, tr_out Z <= t I. The returned Z is shifted until it is exactly
    feasible, and the objective is evaluated there, so the value cannot fall
    below the true norm through solver tolerance.
    """
    Z = cp.Variable((dim * dim, dim * dim), hermitian=True)
    t = cp.Variable()
    constraints = [Z >> 0, Z - J >> 0,
                   t * np.eye(dim) - cp.partial_trace(Z, [dim, dim], axis=1) >> 0]
    problem = cp.Problem(cp.Minimize(t), constraints)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        problem.solve(solver=SOLVER)
    if problem.status != "optimal":
        SOLVER_FLAGS.append(problem.status)
    Zh = (Z.value + Z.value.conj().T) / 2
    shift = max(0.0, -np.linalg.eigvalsh(Zh).min(), -np.linalg.eigvalsh(Zh - J).min())
    Zf = Zh + shift * np.eye(dim * dim)
    traced = np.einsum('iaja->ij', Zf.reshape(dim, dim, dim, dim))
    return 2 * float(np.linalg.eigvalsh((traced + traced.conj().T) / 2).max())


def entangled_input_lower_bound(J, dim):
    """|| (Phi (x) id)(|Omega><Omega|) ||_1, |Omega> maximally entangled."""
    return float(np.abs(np.linalg.eigvalsh(J / dim)).sum())


def fitted_unitary_decoder(kraus_list, dim):
    """V^dagger, where V is the polar part of the leading Choi eigenvector."""
    _, v = np.linalg.eigh(choi(kraus_list, dim))
    # Choi vector sum_i |i> (x) K|i>  ->  K[out, in] = vec[in * dim + out]
    K = v[:, -1].reshape(dim, dim).T
    U, _, Vh = np.linalg.svd(K)
    return (U @ Vh).conj().T


# ---------------------------------------------------------------------------
# Record leakage
# ---------------------------------------------------------------------------

def parity_form_leak(kraus, rounds):
    """Exact largest total-variation distance between the N-round history laws
    of two encoded states, for an instrument whose history effects have the
    parity form F_h = q_h I + g_h P.

    p(h | rho) = q_h + g_h tr(P rho), and tr(P rho) ranges over [-1, 1], so
    the largest distance over pairs of states is sum_h |g_h|. Returns that
    sum and the largest deviation of any F_h from the parity form.
    """
    A = np.array([kraus[s] for s in sorted(kraus)])
    dim = A.shape[1]
    total, residual = 0.0, 0.0
    prefixes = np.eye(dim, dtype=complex)[None]
    for _ in range(rounds - 1):
        prefixes = np.einsum('aij,bjk->abik', A, prefixes).reshape(-1, dim, dim)
    for M0 in prefixes:                                  # earlier rounds
        M = np.einsum('aij,jk->aik', A, M0)              # last round applied last
        F = np.einsum('hji,hjk->hik', M.conj(), M)
        q = np.einsum('hii->h', F).real / dim
        G = F - q[:, None, None] * np.eye(dim)
        g = np.einsum('hij,ji->h', G, PARITY).real / dim
        residual = max(residual, float(np.abs(G - g[:, None, None] * PARITY).max()))
        total += float(np.abs(g).sum())
    return total, residual


def largest_input_dependence(kraus, rounds):
    """Largest entry of the traceless part of any N-round history effect.

    Zero means every history law is the same for every encoded state.
    """
    A = [kraus[s] for s in sorted(kraus)]
    dim = A[0].shape[0]
    products = [np.eye(dim, dtype=complex)]
    for _ in range(rounds):
        products = [K @ M for M in products for K in A]
    worst = 0.0
    for M in products:
        F = M.conj().T @ M
        G = F - np.trace(F).real / dim * np.eye(dim)
        worst = max(worst, float(np.abs(G).max()))
    return worst


# ---------------------------------------------------------------------------
# Repetition-code instrument
# ---------------------------------------------------------------------------

def repetition_code_kraus(theta):
    """Logical Kraus operators of the three-qubit bit-flip code under coherent
    X rotations, Z1Z2 / Z2Z3 syndrome and majority recovery."""
    idx = np.arange(8)
    code = np.zeros((2, 8), dtype=complex)
    code[0, 0b000] = 1
    code[1, 0b111] = 1
    evolved = code.copy()
    for q in range(3):
        evolved = np.cos(theta / 2) * evolved - 1j * np.sin(theta / 2) * evolved[:, idx ^ (1 << q)]

    def bit(v, q):
        return (v >> q) & 1

    syndrome = np.array([(bit(v, 0) ^ bit(v, 1)) | ((bit(v, 1) ^ bit(v, 2)) << 1) for v in idx])
    recovery = {0b00: 0, 0b01: 1 << 0, 0b11: 1 << 1, 0b10: 1 << 2}
    kraus = {}
    for s, r in recovery.items():
        projected = evolved * (syndrome == s)
        kraus[s] = code.conj() @ projected[:, idx ^ r].T
    return kraus


# ---------------------------------------------------------------------------
# Reports
# ---------------------------------------------------------------------------

def validate_sdp():
    """Largest |SDP - 2 sin(phi/2)| for Ad(diag(1, e^{i phi})) against the identity."""
    worst = 0.0
    for phi in (0.1, 0.7, 2.0):
        V = np.diag([1, np.exp(1j * phi)])
        eps = diamond_norm_of_difference(choi([V], 2) - identity_choi(2), 2)
        worst = max(worst, abs(eps - 2 * np.sin(phi / 2)))
    return worst


def report_one_round(inst, passed):
    print("\n=== 1. L = 3 toric code (Özgüler 2026 instrument), one round ===")
    print(" theta  leak_1     eps_id     eps_dec    2sqrt(eps_id)  2sqrt(eps_dec)  "
          "bound/leak  LB(id)     LB(dec)    parity res.")
    rows, holds, worst_residual = [], True, 0.0
    for theta in THETAS:
        kraus = kraus_operators(inst, theta)
        klist = list(kraus.values())
        leak, residual = parity_form_leak(kraus, 1)
        D = choi(klist, 4) - identity_choi(4)
        eps_id, lb_id = diamond_norm_of_difference(D, 4), entangled_input_lower_bound(D, 4)
        W = fitted_unitary_decoder(klist, 4)
        Dd = choi([W @ K for K in klist], 4) - identity_choi(4)
        eps_dec, lb_dec = diamond_norm_of_difference(Dd, 4), entangled_input_lower_bound(Dd, 4)
        ok = leak <= 2 * np.sqrt(lb_dec) and leak <= 2 * np.sqrt(lb_id)
        holds &= ok
        worst_residual = max(worst_residual, residual)
        rows.append((theta, leak, eps_id, eps_dec))
        print(f" {theta:4.2f}   {leak:.3e}  {eps_id:.3e}  {eps_dec:.3e}  {2 * np.sqrt(eps_id):.3e}      "
              f"{2 * np.sqrt(eps_dec):.3e}       {2 * np.sqrt(eps_dec) / leak:9.1f}  "
              f"{lb_id:.3e}  {lb_dec:.3e}  {residual:.1e}  {'bound holds' if ok else 'BOUND FAILS'}")
    (t0, l0, e0, d0), (t1, l1, e1, d1) = rows[0], rows[1]

    def slope(a, b):
        return np.log(b / a) / np.log(t1 / t0)

    print(f"local log-log slopes, theta {t0} -> {t1}: leak {slope(l0, l1):.2f}, "
          f"eps_id {slope(e0, e1):.2f}, eps_dec {slope(d0, d1):.2f}")
    passed.append(("one round: leak <= 2 sqrt(eps) at every angle", holds))
    passed.append(("one round: effects of parity form (leak exact)", worst_residual < TOLERANCE))


def report_rounds(inst, passed):
    print("\n=== 2. L = 3 toric code, N rounds of the same instrument ===")
    print(" theta  N  leak_N     ||Lam^N-id||  2sqrt(||Lam^N-id||)  N*eps_id   2sqrt(N*eps_id)  parity res.")
    holds, additive, worst_residual = True, True, 0.0
    for theta in ROUNDS_THETAS:
        kraus = kraus_operators(inst, theta)
        klist = list(kraus.values())
        eps1 = diamond_norm_of_difference(choi(klist, 4) - identity_choi(4), 4)
        # vec(K rho K^dagger) = (K (x) K*) vec(rho) for row-major vec
        S1 = sum(np.kron(K, K.conj()) for K in klist)
        SN = np.eye(16, dtype=complex)
        for N in range(1, MAX_ROUNDS + 1):
            SN = S1 @ SN
            D = choi_from_superoperator(SN, 4) - identity_choi(4)
            epsN, lbN = diamond_norm_of_difference(D, 4), entangled_input_lower_bound(D, 4)
            leak, residual = parity_form_leak(kraus, N)
            ok_bound = leak <= 2 * np.sqrt(lbN)
            # Subadditivity, tested conservatively: a lower bound on
            # ||Lam^N - id|| against N times the certified upper bound on eps_1.
            ok_add = lbN <= N * eps1
            holds &= ok_bound
            additive &= ok_add
            worst_residual = max(worst_residual, residual)
            print(f" {theta:4.2f}   {N}  {leak:.3e}  {epsN:.3e}     {2 * np.sqrt(epsN):.3e}            "
                  f"{N * eps1:.3e}  {2 * np.sqrt(N * eps1):.3e}        {residual:.1e}  "
                  f"{'bounds hold' if ok_bound and ok_add else 'BOUND FAILS'}")
    passed.append(("N rounds: leak_N <= 2 sqrt(||Lam^N - id||)", holds))
    passed.append(("N rounds: ||Lam^N - id|| <= N eps_id", additive))
    passed.append(("N rounds: effects of parity form (leak exact)", worst_residual < TOLERANCE))


def report_repetition(passed):
    print("\n=== 3. Three-qubit repetition code: exact blindness without exact preservation ===")
    print(f" theta  input dependence, all {REPETITION_ROUNDS}-round histories   eps_id     "
          "LB(id)     completeness")
    blind, damaged = True, True
    for theta in THETAS:
        kraus = repetition_code_kraus(theta)
        klist = list(kraus.values())
        completeness = float(np.abs(sum(K.conj().T @ K for K in klist) - np.eye(2)).max())
        dependence = largest_input_dependence(kraus, REPETITION_ROUNDS)
        D = choi(klist, 2) - identity_choi(2)
        eps, lb = diamond_norm_of_difference(D, 2), entangled_input_lower_bound(D, 2)
        blind &= dependence < TOLERANCE
        damaged &= lb > 0
        print(f" {theta:4.2f}   {dependence:.1e}                               "
              f"{eps:.3e}  {lb:.3e}  {completeness:.0e}")
    passed.append(("repetition code: every history law independent of the encoded state", blind))
    passed.append(("repetition code: logical channel differs from the identity", damaged))


def report_sentinel():
    print("\n=== 4. Sentinel qubit: sign discrimination (+theta vs -theta), Y-basis readout ===")
    print(" theta  per-shot TV = sin(theta)   TV after N shots, N = " +
          ", ".join(str(n) for n in SENTINEL_SHOTS))
    for theta in THETAS:
        p, q = (1 + np.sin(theta)) / 2, (1 - np.sin(theta)) / 2
        tvs = [0.5 * sum(comb(n, k) * abs(p ** k * q ** (n - k) - q ** k * p ** (n - k))
                         for k in range(n + 1)) for n in SENTINEL_SHOTS]
        print(f" {theta:4.2f}   {np.sin(theta):.4f}                     " +
              "   ".join(f"{t:.4f}" for t in tvs))
    print(" The passive syndrome record's sign profile is identically 0 "
          "(see ozguler_prop1_check.py).")


def main() -> None:
    print("=== Run configuration ===")
    print(f"script version = {SCRIPT_VERSION}")
    print(f"Python version = {platform.python_version()}")
    print(f"NumPy version = {np.__version__}")
    try:
        import clarabel
        solver_version = clarabel.__version__
    except (ImportError, AttributeError):
        solver_version = "unknown"
    print(f"CVXPY version = {cp.__version__}")
    print(f"solver = {SOLVER} {solver_version}")
    print(f"thetas = {THETAS}; N-round thetas = {ROUNDS_THETAS}; max rounds = {MAX_ROUNDS}")
    print("diamond norm unnormalised, range [0, 2]; leak = total-variation distance, range [0, 1]")

    passed = []
    validation = validate_sdp()
    print(f"\nSDP validation: max |SDP - 2 sin(phi/2)| over three unitaries = {validation:.1e}")
    passed.append(("SDP agrees with closed form (to 1e-6)", validation < 1e-6))

    inst = build_instrument(3)
    report_one_round(inst, passed)
    report_rounds(inst, passed)
    report_repetition(passed)
    report_sentinel()

    print("\n=== Summary ===")
    for label, ok in passed:
        print(f"{'PASS' if ok else 'FAIL'}  {label}")
    statuses = ", ".join(sorted(set(SOLVER_FLAGS))) or "none"
    print(f"SDP solves not reported 'optimal' by the solver: {len(SOLVER_FLAGS)} ({statuses}). "
          "Reported distances are certified upper bounds and verdicts use lower bounds, "
          "so solver status does not affect the results.")


if __name__ == "__main__":
    main()
