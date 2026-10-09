#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
prop52_numerical_check.py
=========================

Version:
    0.6.1 (2026-10-08)

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

It then compares the record with the complementary output it is read from.
The complementary channel keeps what the process passes to its environment:
here, the syndrome register with the coherences between syndrome outcomes
intact. The classical record is that output with the coherences discarded.
For four angles the script maximises, over pairs of pure encoded states, the
trace distance between the two complementary outputs, and compares it with
the record's leak and with the bound 2 sqrt(eps_dec). The bound of
Kretschmann, Schlingemann and Werner (2008) applies to the complementary
output, so record leak <= complement distance <= 2 sqrt(eps) must hold.
The maximisation is a local search from fixed random starts, so the reported
complement distance is a lower bound on the true maximum.

It then separates pointwise from uniform repair. The recoverability defect
in Proposition 5.2 allows a different repair for each noise channel. When
the sign of the rotation is unknown, and the syndrome record cannot reveal it
(Özgüler 2026, Proposition 1), a repair must be chosen without it. For four
angles the script reports the fitted unitary repair at the right and the
wrong sign and no repair, and then brackets four optimal repairs over all
channels: for a known sign and for both signs, each with the syndrome
history discarded (one repair after the standard correction) and with the
history used (one repair per class of syndromes whose branch maps agree up to
weight). Each optimum is bracketed by an upper bound, the defect achieved by
the solver's repairs made exact channels, and a lower bound, the dual
objective at a dual point built from the solver's dual variables and made
exactly feasible. One round only. For the history-used repairs, a second
lower bound is computed from the same dual point with the dual constraint
made feasible separately for every one of the 256 syndromes, so it bounds the
optimum over all history-conditioned policies without assuming that the
grouping into classes is exact.

For the repetition code it checks that every three-round history law is
independent of the encoded state, although the logical channel is not the
identity: the converse of the proposition fails.

For contrast, it reports the sign-discrimination profile of a single
sentinel qubit measured in the Y basis (per-shot total variation sin theta),
which the passive syndrome record lacks entirely (Özgüler 2026,
Proposition 1; checked by ozguler_prop1_check.py).

Numerical safeguards
--------------------
Every reported diamond distance is an upper bound constructed so as not to
depend on the solver's tolerance: the SDP's dual point is shifted until it satisfies the
dual constraints, as judged by a floating-point eigenvalue computation, and
the dual objective is evaluated there. Every
check of the inequality is made at a lower bound that does not come from the
solver at all: the trace norm of the output for one explicit input, the
maximally entangled state. Both are numerically constructed estimates,
computed in double-precision floating point; the margins by which the checks
pass far exceed the solver residuals, but the bounds have not been certified
by interval arithmetic or validated eigenvalue bounds. The SDP is validated against the closed form
2 sin(phi/2) for the unitary channel diag(1, e^{i phi}).

What the checks test. Proposition 5.2 bounds the leak by 2 sqrt(eps(Lambda)),
where eps(Lambda) is an infimum over all repair channels. The script cannot
compute that infimum. It uses eps_dec, the defect after one fitted unitary
repair, which bounds eps(Lambda) from above. The check leak <= 2 sqrt(eps_dec)
is therefore a necessary consequence of the proposition, and a violation
would refute it; but passing it does not confirm the bound at the infimum.

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
from scipy.optimize import minimize
import warnings
from math import comb

import cvxpy as cp
import numpy as np

from ozguler_prop1_check import build_instrument, kraus_operators

SCRIPT_VERSION = "0.6.1"

THETAS = (0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5)   # one-round angle grid
ROUNDS_THETAS = (0.1, 0.3)                       # angles for the N-round check
MAX_ROUNDS = 3                                   # N = 1, ..., MAX_ROUNDS
REPETITION_ROUNDS = 3                            # history length, repetition code
SENTINEL_SHOTS = (1, 10, 100)
SOLVER = "CLARABEL"
TOLERANCE = 1e-12        # rounding tolerance for the pass/fail summary

X = np.array([[0, 1], [1, 0]], dtype=complex)
PARITY = np.kron(X, X)   # product of the two logical X operators

COMPLEMENT_THETAS = (0.05, 0.1, 0.2, 0.3)
COMPLEMENT_STARTS = 4
COMPLEMENT_SEED = 0

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
    """Upper bound, independent of solver tolerance, on the unnormalised diamond norm of a difference
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
            # ||Lam^N - id|| against N times the upper bound on eps_1.
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
    print("\n=== 5. Three-qubit repetition code: exact blindness without exact preservation ===")
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


def complement_distance(kraus, starts=COMPLEMENT_STARTS, seed=COMPLEMENT_SEED):
    """Largest trace distance found between complementary outputs of two pure
    encoded states. The complementary output of rho is the matrix
    [tr(K_s rho K_t^dagger)]_{s,t} on the syndrome register; the record is its
    diagonal. Local search from fixed random starts: a lower bound on the maximum."""
    A = np.array([kraus[s] for s in sorted(kraus)])
    dim = A.shape[1]

    def unit(x):
        v = x[:dim] + 1j * x[dim:]
        return v / np.linalg.norm(v)

    def negative_distance(x):
        a, b = unit(x[:2 * dim]), unit(x[2 * dim:])
        X = np.einsum('sij,jk,tik->st', A, np.outer(a, a.conj()) - np.outer(b, b.conj()), A.conj())
        return -0.5 * float(np.abs(np.linalg.eigvalsh((X + X.conj().T) / 2)).sum())

    rng = np.random.default_rng(seed)
    best = 0.0
    for _ in range(starts):
        result = minimize(negative_distance, rng.normal(size=4 * dim), method="L-BFGS-B")
        best = max(best, -float(result.fun))
    return best


def report_complement(inst, passed):
    print("\n=== 3. L = 3 toric code: what the environment receives versus what the record keeps ===")
    print(" theta  complement distance (>=)  record leak   2sqrt(eps_dec)  complement/record")
    ordered, rows = True, []
    for theta in COMPLEMENT_THETAS:
        kraus = kraus_operators(inst, theta)
        klist = list(kraus.values())
        leak, _ = parity_form_leak(kraus, 1)
        comp = complement_distance(kraus)
        W = fitted_unitary_decoder(klist, 4)
        eps_dec = diamond_norm_of_difference(choi([W @ K for K in klist], 4) - identity_choi(4), 4)
        ok = leak <= comp + 1e-12 and comp <= 2 * np.sqrt(eps_dec)
        ordered &= ok
        rows.append((theta, comp))
        print(f" {theta:4.2f}   {comp:.3e}                 {leak:.3e}   {2 * np.sqrt(eps_dec):.3e}       "
              f"{comp / leak:9.1f}  {'ordering holds' if ok else 'ORDERING FAILS'}")
    (t0, c0), (t1, c1) = rows[0], rows[1]
    print(f"local log-log slope of the complement distance, theta {t0} -> {t1}: "
          f"{np.log(c1 / c0) / np.log(t1 / t0):.2f}")
    passed.append(("record leak <= complement distance <= 2 sqrt(eps_dec)", ordered))


def compose(a_list, b_list):
    """Kraus operators of a after b."""
    return [A @ B for A in a_list for B in b_list]


def kraus_from_choi(J, dim):
    """Kraus operators of the channel with Choi matrix J (input factor first)."""
    w, v = np.linalg.eigh((J + J.conj().T) / 2)
    return [np.sqrt(max(x, 0.0)) * v[:, i].reshape(dim, dim).T for i, x in enumerate(w) if x > 1e-14]


def choi_map(kraus_list, dim):
    """Matrix T with vec(Choi(D o Lambda)) = T vec(Choi(D)), row-major vec,
    for the (possibly trace-decreasing) map Lambda with these Kraus operators."""
    Lt = np.zeros((dim, dim, dim, dim), dtype=complex)
    for i in range(dim):
        for j in range(dim):
            E = np.zeros((dim, dim), dtype=complex)
            E[i, j] = 1
            Lt[i, j] = sum(K @ E @ K.conj().T for K in kraus_list)
    I = np.eye(dim)
    return np.einsum('ijkl,pr,qs->ipjqkrls', Lt, I, I).reshape(dim ** 4, dim ** 4)


def adjoint_map(T, W, dim):
    """H with tr(J H) = tr(Choi-image(J) W), for the map vec(J) -> T vec(J)."""
    n = dim * dim
    return (T.T @ W.T.reshape(-1)).reshape(n, n).T


def syndrome_classes(k_plus, k_minus):
    """Group syndromes whose Kraus operators, at both signs, are proportional to
    those of a representative with coefficients of equal modulus: their branch
    maps are then the same up to a common weight, so one repair serves them all."""
    reps = []
    for s in sorted(k_plus):
        a, b = k_plus[s], k_minus[s]
        norm = np.linalg.norm(a)
        if norm < 1e-14:
            continue
        for ra, rb, members in reps:
            c = np.vdot(ra, a) / np.vdot(ra, ra)
            c2 = np.vdot(rb, b) / np.vdot(rb, rb)
            if (np.allclose(a, c * ra, atol=1e-12 * norm) and np.allclose(b, c2 * rb, atol=1e-12 * norm)
                    and abs(abs(c) - abs(c2)) < 1e-9):
                members.append(s)
                break
        else:
            reps.append((a, b, [s]))
    return [members for _, _, members in reps]


def _largest_trace_Y(H, dim):
    """Largest tr Y with H - Y (x) I >= 0, made exactly feasible by a shift."""
    Y = cp.Variable((dim, dim), hermitian=True)
    problem = cp.Problem(cp.Maximize(cp.real(cp.trace(Y))), [H - cp.kron(Y, np.eye(dim)) >> 0])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        problem.solve(solver=SOLVER)
    Yv = (Y.value + Y.value.conj().T) / 2
    shift = max(0.0, -np.linalg.eigvalsh(H - np.kron(Yv, np.eye(dim))).min())
    return Yv - shift * np.eye(dim)


def optimal_repair(groups, dim, all_branches=None):
    """Bracket min over repair policies of max over signs of ||sum_r D_r o I_r^s - id||.

    groups: one entry per class of syndromes sharing a repair, each a list over
    signs of Kraus lists for that class's branch map. A single group holding all
    branches is a repair applied after the history is discarded; one group per
    syndrome class is a repair chosen according to the history.

    all_branches: optional list, one entry per individual syndrome, each a list
    over signs of Kraus lists. If given, a third value is returned: a lower
    bound over policies with one repair per individual syndrome, from the same
    dual point, with the dual constraint made feasible syndrome by syndrome. It
    assumes nothing about proportionality within classes.

    Returns (lower, upper) or (lower, upper, lower_all). upper is the defect achieved by the solver's repairs,
    made exact channels and evaluated with diamond_norm_of_difference. lower is
    the dual objective at a dual point built from the solver's dual variables and
    made exactly feasible (Watrous's primal form of the diamond norm, with the
    channel constraint dualised); by weak duality it bounds the optimum below.
    Both are floating-point computations."""
    n = dim * dim
    n_signs = len(groups[0])
    maps = [[choi_map(group[s], dim) for s in range(n_signs)] for group in groups]
    Om = identity_choi(dim)
    JD = [cp.Variable((n, n), hermitian=True) for _ in groups]
    t = cp.Variable()
    constraints = []
    for J in JD:
        constraints += [J >> 0, cp.partial_trace(J, [dim, dim], axis=1) == np.eye(dim)]
    w_cons, r_cons = [], []
    for s in range(n_signs):
        out = sum(cp.reshape(T[s] @ cp.vec(J, order='C'), (n, n), order='C') for T, J in zip(maps, JD))
        Z = cp.Variable((n, n), hermitian=True)
        c_w = Z - (out - Om) >> 0
        c_r = t * np.eye(dim) - cp.partial_trace(Z, [dim, dim], axis=1) >> 0
        constraints += [Z >> 0, c_w, c_r]
        w_cons.append(c_w)
        r_cons.append(c_r)
    problem = cp.Problem(cp.Minimize(t), constraints)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        problem.solve(solver=SOLVER)
    if problem.status != "optimal":
        SOLVER_FLAGS.append(problem.status)
    # Upper bound: the solver's repairs, made exact channels.
    repairs = []
    for J in JD:
        Jv = (J.value + J.value.conj().T) / 2
        w, v = np.linalg.eigh(Jv)
        Jv = (v * np.clip(w, 0, None)) @ v.conj().T
        A = np.einsum('iaja->ij', Jv.reshape(dim, dim, dim, dim))
        aw, av = np.linalg.eigh((A + A.conj().T) / 2)
        T = np.kron((av / np.sqrt(aw)) @ av.conj().T, np.eye(dim))
        repairs.append(kraus_from_choi(T @ Jv @ T.conj().T, dim))
    upper = 0.0
    for s in range(n_signs):
        total = []
        for repair, group in zip(repairs, groups):
            total += compose(repair, group[s])
        upper = max(upper, diamond_norm_of_difference(choi(total, dim) - Om, dim))
    # Lower bound: a dual point made exactly feasible.
    W_list, rho_list = [], []
    for c_w, c_r in zip(w_cons, r_cons):
        W = (c_w.dual_value + c_w.dual_value.conj().T) / 2
        w, v = np.linalg.eigh(W)
        W = (v * np.clip(w, 0, None)) @ v.conj().T
        rho = (c_r.dual_value + c_r.dual_value.conj().T) / 2
        w, v = np.linalg.eigh(rho)
        rho = (v * np.clip(w, 0, None)) @ v.conj().T
        rho = rho + max(0.0, -np.linalg.eigvalsh(np.kron(rho, np.eye(dim)) - W).min()) * np.eye(dim)
        W_list.append(W)
        rho_list.append(rho)
    scale = 1 / sum(np.trace(r).real for r in rho_list)
    W_list = [scale * W for W in W_list]
    lower = -sum(np.trace(Om @ W).real for W in W_list)
    base = lower
    for T in maps:
        H = sum(adjoint_map(T[s], W_list[s], dim) for s in range(n_signs))
        lower += np.trace(_largest_trace_Y((H + H.conj().T) / 2, dim)).real
    if all_branches is None:
        return 2 * float(lower), float(upper)
    lower_all = base
    for branch in all_branches:
        H = sum(adjoint_map(choi_map(branch[s], dim), W_list[s], dim) for s in range(n_signs))
        lower_all += np.trace(_largest_trace_Y((H + H.conj().T) / 2, dim)).real
    return 2 * float(lower), float(upper), 2 * float(lower_all)


def report_sign_repair(inst, passed):
    print("\n=== 4. L = 3 toric code: repair chosen without knowing the sign ===")
    print(" Fitted unitary repair at the right and the wrong sign, and no repair (upper bounds):")
    print(" theta  right sign   wrong sign   no repair")
    rows = {}
    for theta in COMPLEMENT_THETAS:
        k_plus = list(kraus_operators(inst, theta).values())
        k_minus = list(kraus_operators(inst, -theta).values())
        W = fitted_unitary_decoder(k_plus, 4)
        I = identity_choi(4)
        right = diamond_norm_of_difference(choi([W @ K for K in k_plus], 4) - I, 4)
        wrong = diamond_norm_of_difference(choi([W @ K for K in k_minus], 4) - I, 4)
        none = diamond_norm_of_difference(choi(k_minus, 4) - I, 4)
        rows[theta] = (right, wrong, none)
        print(f" {theta:4.2f}   {right:.3e}    {wrong:.3e}    {none:.3e}")
    print(" Best repair over all channels, as [lower bound, achieved upper bound]:")
    print(" theta  sign known,            both signs,            sign known,            both signs,")
    print("        history discarded      history discarded      history used           history used")
    gap, fitted_optimal, nothing_better, history_widens = True, True, True, True
    all_rows, grouping_safe = {}, True
    for theta in COMPLEMENT_THETAS:
        kp = kraus_operators(inst, theta)
        km = kraus_operators(inst, -theta)
        lp, lm = list(kp.values()), list(km.values())
        classes = syndrome_classes(kp, km)
        p_lo, p_up = optimal_repair([[lp]], 4)
        u_lo, u_up = optimal_repair([[lp, lm]], 4)
        keys = sorted(kp)
        ph_lo, ph_up, ph_all = optimal_repair([[[kp[s] for s in c]] for c in classes], 4,
                                              [[[kp[s]]] for s in keys])
        uh_lo, uh_up, uh_all = optimal_repair([[[kp[s] for s in c], [km[s] for s in c]] for c in classes], 4,
                                              [[[kp[s]], [km[s]]] for s in keys])
        all_rows[theta] = (ph_all, ph_up, uh_all, uh_up, len(keys))
        right, wrong, none = rows[theta]
        gap &= min(uh_lo, uh_all) > 2 * ph_up
        grouping_safe &= uh_all > 0.99 * uh_lo
        fitted_optimal &= right <= 1.05 * p_lo
        nothing_better &= none <= 1.03 * u_lo
        history_widens &= uh_up <= u_up * (1 + 1e-6) and uh_lo / ph_up >= u_up / p_lo
        print(f" {theta:4.2f}   [{p_lo:.3e}, {p_up:.3e}]  [{u_lo:.3e}, {u_up:.3e}]  "
              f"[{ph_lo:.3e}, {ph_up:.3e}]  [{uh_lo:.3e}, {uh_up:.3e}]   ({len(classes)} classes)")
    print(" (history used: one repair per class of syndromes whose branches agree up to weight;")
    print("  one round only; bounds are floating-point computations, not interval-certified)")
    print(" History used, lower bounds over every individual syndrome, without grouping")
    print(" (same dual point, dual constraint made feasible syndrome by syndrome; upper bounds as above,")
    print("  since one repair per class is a policy over individual syndromes):")
    print(" theta  sign known, history used     both signs, history used")
    for theta in COMPLEMENT_THETAS:
        ph_all, ph_up, uh_all, uh_up, n_syn = all_rows[theta]
        print(f" {theta:4.2f}   [{ph_all:.4e}, {ph_up:.4e}]   [{uh_all:.4e}, {uh_up:.4e}]   ({n_syn} syndromes)")
    print(" (a negative lower bound is uninformative)")
    passed.append(("sign-blind repair, history used, more than twice the known-sign optimum", gap))
    passed.append(("sign-blind lower bound, history used, holds over all syndromes without grouping (within 1 per cent)", grouping_safe))
    passed.append(("fitted unitary repair within 5 per cent of the known-sign lower bound", fitted_optimal))
    passed.append(("no repair within 3 per cent of the sign-blind lower bound, history discarded", nothing_better))
    passed.append(("using the history helps the sign-blind repair less than the known-sign one", history_widens))


def report_sentinel():
    print("\n=== 6. Sentinel qubit: sign discrimination (+theta vs -theta), Y-basis readout ===")
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
    report_complement(inst, passed)
    report_sign_repair(inst, passed)
    report_repetition(passed)
    report_sentinel()

    print("\n=== Summary ===")
    for label, ok in passed:
        print(f"{'PASS' if ok else 'FAIL'}  {label}")
    statuses = ", ".join(sorted(set(SOLVER_FLAGS))) or "none"
    print(f"SDP solves not reported 'optimal' by the solver: {len(SOLVER_FLAGS)} ({statuses}). "
          "Reported distances are upper bounds constructed independently of solver tolerance, "
          "and verdicts use explicit lower bounds; both are floating-point estimates, not "
          "interval-certified.")


if __name__ == "__main__":
    main()
