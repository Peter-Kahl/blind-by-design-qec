#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
infidelity_and_calibration_check.py
===================================

Version:
    0.10.0 (2026-10-31)

Supplementary research code for:

    Peter Kahl, 'The Cost of Sign-Blind Recovery under Coherent Noise: Worked Numerical Cases'
    (2026), §4.

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
Section 4 of prop52_numerical_check.py measures the cost of not knowing the
sign of a coherent rotation in diamond norm, a worst-case measure that
weighs coherent errors heavily. This script asks three further questions of
the same L = 3 toric instrument (Özgüler 2026):

1. Infidelity. What are the best one-round repairs, for the same four
   states of knowledge, when the objective is entanglement infidelity
   (1 - entanglement fidelity, an average-case measure) instead of diamond
   norm?
2. Accumulation. Coherent errors add in amplitude, so their infidelity grows
   roughly as the square of the number of rounds. How does the infidelity
   after N rounds compare when a repair is applied every round: no repair; a
   unitary repair fitted at the right sign and at the wrong sign; the
   diamond-optimal history-used repairs with and without the sign; and, for
   reference, an idealised logical Pauli twirl of the unrepaired channel?
3. Calibration error. The sign-known repairs assume theta is known exactly.
   How large is the one-round defect, and the infidelity after up to 300
   rounds, when the repair is fitted at theta (1 + eta) but the channel is
   at theta?

4. Sign error. Branch undo needs the sign, and applied at the wrong sign
   it is worse than no repair. If the sign is estimated once and used for
   all rounds, what is the expected infidelity after 300 rounds for a given
   probability of a wrong sign, and how many independent sentinel shots
   does a majority vote need to reach that probability?

Version 0.9.0 adds the branch-undo repair throughout: for each syndrome,
the inverse of the nearest unitary to that branch's Kraus operator, fitted
with the sign known. Each branch is a coherent logical rotation
(branch_mechanism_check.py), so this simple repair removes almost all of
the one-round infidelity; it is a feasible sign-known, history-used policy,
and it tightens the upper end of the infidelity brackets.

Methods
-------
Entanglement fidelity of a channel Phi on dimension d is
<Omega| (Phi x id)(|Omega><Omega|) |Omega> with |Omega> maximally entangled,
computed as tr(Omega J) / d^2 for the unnormalised Choi matrix J.

Infidelity optima (part 1) are reported as brackets. The upper end is the
infidelity actually achieved by the solver's repairs, made exact channels
and evaluated directly, or by no repair if that is better (no repair is a
member of every class of policies). The lower end is a dual bound: for any
weights lambda_s >= 0 summing to 1, the best worst-sign fidelity is at most
sum over classes of max over repair channels of tr(J H_r), with
H_r = sum_s lambda_s G_r^s; and max over channels of tr(J H) <= tr Y for any
Y with Y (x) I >= H. Y is found by semidefinite programme and shifted until
it is exactly feasible. lambda = (1/2, 1/2) for the sign-blind problem, by
the symmetry between the signs. Both ends are floating-point constructions,
not interval-certified.

Repairs applied every round (part 2) use only that round's syndrome, so the
N-round channel is the N-th power of the repaired one-round channel. The
twirl is applied to the logical channel, averaging over the 16 two-qubit
Pauli operators; it is an idealisation of randomised compiling, which in
practice twirls physical operations, and it leaves one-round entanglement
fidelity unchanged.

Diamond-optimal repairs (parts 2 and 3) are computed with the same
semidefinite programme as optimal_repair in prop52_numerical_check.py, and
returned as exact channels. Defects are evaluated with
diamond_norm_of_difference (an upper bound) and, as a check, the
entangled-input lower bound.

Sign error (part 4): with probability q of a wrong sign, fixed for the
whole run, the expected N-round infidelity is (1 - q) I_right + q I_wrong,
where I_right and I_wrong are the infidelities of branch undo applied at
the right and the wrong sign. A sentinel qubit prepared in |0>, exposed to
the rotation and measured in the Y basis, gives the outcome indicating the
true sign with probability (1 + sin theta) / 2 per shot. With n shots (n
odd) and a majority vote, the probability of a wrong sign is the binomial
tail P(X <= (n - 1) / 2), X ~ Binomial(n, (1 + sin theta) / 2). Shots are
assumed independent, the sign is assumed the same on the sentinel and the
memory, and drift is not modelled.

Run
---
    python infidelity_and_calibration_check.py

Imports prop52_numerical_check.py and ozguler_prop1_check.py from the same
folder. About 10 minutes on the reference machine.
"""

from __future__ import annotations

import itertools
import platform
import time
import warnings

import cvxpy as cp
import numpy as np

from ozguler_prop1_check import build_instrument, kraus_operators
from prop52_numerical_check import (SOLVER, SOLVER_FLAGS, adjoint_map, choi,
                                    choi_from_superoperator, choi_map,
                                    diamond_norm_of_difference,
                                    entangled_input_lower_bound,
                                    fitted_unitary_decoder, identity_choi,
                                    kraus_from_choi, syndrome_classes)

SCRIPT_VERSION = "0.10.0"

DIM = 4
N2 = DIM * DIM
OMEGA = identity_choi(DIM)

FIDELITY_THETAS = (0.05, 0.1, 0.2, 0.3)
ACCUMULATION_THETAS = (0.05, 0.1)
ROUNDS = (1, 10, 30, 100, 300)
CALIBRATION_THETAS = (0.05, 0.1)
ETAS = (0.0, 0.001, 0.01, 0.1, -0.01)
CALIBRATION_ROUNDS = (1, 100, 300)

PAULI_1 = [np.eye(2), np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]),
           np.diag([1, -1])]
PAULIS = [np.kron(a, b).astype(complex) for a, b in itertools.product(PAULI_1, PAULI_1)]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def solve(problem):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        problem.solve(solver=SOLVER)
    if problem.status != "optimal":
        SOLVER_FLAGS.append(problem.status)


def exact_channel(Jvalue):
    """Kraus operators of the solver's Choi matrix made an exact channel."""
    Jv = (Jvalue + Jvalue.conj().T) / 2
    w, v = np.linalg.eigh(Jv)
    Jv = (v * np.clip(w, 0, None)) @ v.conj().T
    A = np.einsum('iaja->ij', Jv.reshape(DIM, DIM, DIM, DIM))
    aw, av = np.linalg.eigh((A + A.conj().T) / 2)
    T = np.kron((av / np.sqrt(aw)) @ av.conj().T, np.eye(DIM))
    return kraus_from_choi(T @ Jv @ T.conj().T, DIM)


def superop(kraus_list):
    """Row-major superoperator: vec(K rho K^dagger) = (K x K*) vec(rho)."""
    return sum(np.kron(K, K.conj()) for K in kraus_list)


def ent_fidelity_superop(S):
    return float(np.real(np.trace(OMEGA @ choi_from_superoperator(S, DIM)))) / N2


def ent_fidelity(kraus_list):
    return float(np.real(np.trace(OMEGA @ choi(kraus_list, DIM)))) / N2


def repaired_kraus(policy, branches):
    """Kraus operators of the channel when branch r is followed by repair r."""
    return [R @ K for repair, group in zip(policy, branches) for R in repair for K in group]


def class_data(inst, theta):
    kp, km = kraus_operators(inst, theta), kraus_operators(inst, -theta)
    classes = syndrome_classes(kp, km)
    return kp, km, classes


# ---------------------------------------------------------------------------
# Part 1: infidelity-optimal repairs, with brackets
# ---------------------------------------------------------------------------

def fidelity_bracket(groups):
    """groups: list over repair slots, each a list over signs of Kraus lists.
    Returns (lower, upper) on the best worst-sign entanglement infidelity."""
    n_signs = len(groups[0])
    maps = [[choi_map(g[s], DIM) for s in range(n_signs)] for g in groups]
    JD = [cp.Variable((N2, N2), hermitian=True) for _ in groups]
    t = cp.Variable()
    cons = []
    for J in JD:
        cons += [J >> 0, cp.partial_trace(J, [DIM, DIM], axis=1) == np.eye(DIM)]
    for s in range(n_signs):
        out = sum(cp.reshape(T[s] @ cp.vec(J, order='C'), (N2, N2), order='C')
                  for T, J in zip(maps, JD))
        cons.append(t <= cp.real(cp.trace(OMEGA @ out)) / N2)
    solve(cp.Problem(cp.Maximize(t), cons))

    # Upper end: achieved by exact channels, or by no repair if better.
    policy = [exact_channel(J.value) for J in JD]
    achieved = min(ent_fidelity(repaired_kraus(policy, [g[s] for g in groups]))
                   for s in range(n_signs))
    none = min(ent_fidelity([K for g in groups for K in g[s]]) for s in range(n_signs))
    upper = 1 - max(achieved, none)

    # Lower end: dual bound with equal weights over signs.
    lam = [1 / n_signs] * n_signs
    fmax = 0.0
    for T in maps:
        H = sum(lam[s] * adjoint_map(T[s], OMEGA / N2, DIM) for s in range(n_signs))
        H = (H + H.conj().T) / 2
        Y = cp.Variable((DIM, DIM), hermitian=True)
        solve(cp.Problem(cp.Minimize(cp.real(cp.trace(Y))), [cp.kron(Y, np.eye(DIM)) - H >> 0]))
        Yv = (Y.value + Y.value.conj().T) / 2
        shift = max(0.0, -np.linalg.eigvalsh(np.kron(Yv, np.eye(DIM)) - H).min())
        fmax += float(np.real(np.trace(Yv))) + DIM * shift
    lower = 1 - fmax
    return lower, upper


def report_infidelity(inst, passed):
    print("\n=== 1. One round: best repairs by entanglement infidelity ===")
    print(" [lower bound, achieved upper bound]; 'no repair' is exact")
    print(" theta  no repair    sign known,              both signs,              "
          "sign known,              both signs,")
    print("                     history discarded        history discarded        "
          "history used             history used")
    rows, undo_rows = {}, {}
    for theta in FIDELITY_THETAS:
        kp, km, classes = class_data(inst, theta)
        lp, lm = list(kp.values()), list(km.values())
        none = 1 - ent_fidelity(lp)
        kd = fidelity_bracket([[lp]])
        bd = fidelity_bracket([[lp, lm]])
        kh = fidelity_bracket([[[kp[s] for s in c]] for c in classes])
        undo = 1 - ent_fidelity(apply_policy(branch_undo_policy(kp), kp))
        kh = (kh[0], min(kh[1], undo))      # branch undo is a feasible sign-known policy
        bh = fidelity_bracket([[[kp[s] for s in c], [km[s] for s in c]] for c in classes])
        rows[theta] = (none, kd, bd, kh, bh)
        undo_rows[theta] = undo
        print(f" {theta:4.2f}   {none:.3e}    [{kd[0]:.3e}, {kd[1]:.3e}]   [{bd[0]:.3e}, {bd[1]:.3e}]   "
              f"[{kh[0]:.3e}, {kh[1]:.3e}]   [{bh[0]:.3e}, {bh[1]:.3e}]")
    print(" (sign known, history used: upper end is the better of the solver's repairs and the")
    print("  branch-undo repair; branch undo alone gives " +
          ", ".join(f"{undo_rows[t]:.3e} at {t}" for t in FIDELITY_THETAS) + ")")
    print(" Ratio, sign-blind lower bound / sign-known upper bound (history used):")
    for theta, (none, kd, bd, kh, bh) in rows.items():
        print(f"  theta {theta:4.2f}: {bh[0] / kh[1]:.1f}   "
              f"(sign known, history discarded, improves on no repair by {100 * (1 - kd[0] / none):.1f} "
              "per cent at most)")
    passed.append(("infidelity: sign-blind, history used, more than twice the known-sign optimum",
                   all(bh[0] > 2 * kh[1] for _, _, _, kh, bh in rows.values())))
    passed.append(("infidelity: brackets ordered (lower <= upper)",
                   all(b[0] <= b[1] + 1e-12 for r in rows.values() for b in r[1:])))


# ---------------------------------------------------------------------------
# Diamond-optimal history-used repairs, returned as policies
# ---------------------------------------------------------------------------

def diamond_policy(groups):
    """Repairs minimising the worst-sign diamond defect (as optimal_repair),
    returned as exact channels, one per group."""
    n_signs = len(groups[0])
    maps = [[choi_map(g[s], DIM) for s in range(n_signs)] for g in groups]
    JD = [cp.Variable((N2, N2), hermitian=True) for _ in groups]
    t = cp.Variable()
    cons = []
    for J in JD:
        cons += [J >> 0, cp.partial_trace(J, [DIM, DIM], axis=1) == np.eye(DIM)]
    for s in range(n_signs):
        out = sum(cp.reshape(T[s] @ cp.vec(J, order='C'), (N2, N2), order='C')
                  for T, J in zip(maps, JD))
        Z = cp.Variable((N2, N2), hermitian=True)
        cons += [Z >> 0, Z - (out - OMEGA) >> 0,
                 t * np.eye(DIM) - cp.partial_trace(Z, [DIM, DIM], axis=1) >> 0]
    solve(cp.Problem(cp.Minimize(t), cons))
    return [exact_channel(J.value) for J in JD]


def per_syndrome(policy, classes):
    out = {}
    for repair, c in zip(policy, classes):
        for s in c:
            out[s] = repair
    return out


def apply_policy(syndrome_repairs, kraus):
    return [R @ kraus[s] for s in kraus for R in syndrome_repairs[s]]


def branch_undo_policy(kraus):
    """Sign-known, history-used repair: for each syndrome, the inverse of the
    nearest unitary to that branch's Kraus operator."""
    out = {}
    for s, K in kraus.items():
        u, _, vh = np.linalg.svd(K)
        out[s] = [(u @ vh).conj().T]
    return out


def twirl_superop(S):
    return sum(np.kron(P.conj().T, P.T) @ S @ np.kron(P, P.conj()) for P in PAULIS) / len(PAULIS)


# ---------------------------------------------------------------------------
# Part 2: accumulation over rounds
# ---------------------------------------------------------------------------

def report_accumulation(inst, passed):
    print("\n=== 2. Infidelity after N rounds, repair applied every round ===")
    print(" Policies: none; fitted unitary at the right / wrong sign; diamond-optimal history-used repair")
    print(" with the sign known; branch undo with the sign known; history-used repair minimising the")
    print(" worst-sign diamond defect for both signs (worse sign shown); logical Pauli twirl of the")
    print(" unrepaired channel (idealised). Exact powers of the one-round channel.")
    ok_growth, ok_order = True, True
    for theta in ACCUMULATION_THETAS:
        kp, km, classes = class_data(inst, theta)
        lp, lm = list(kp.values()), list(km.values())
        W = fitted_unitary_decoder(lp, DIM)
        known = per_syndrome(diamond_policy([[[kp[s] for s in c]] for c in classes]), classes)
        blind = per_syndrome(diamond_policy([[[kp[s] for s in c], [km[s] for s in c]]
                                             for c in classes]), classes)
        channels = {
            "none": superop(lp),
            "fitted, right sign": superop([W @ K for K in lp]),
            "fitted, wrong sign": superop([W @ K for K in lm]),
            "known sign, history": superop(apply_policy(known, kp)),
            "known sign, branch undo": superop(apply_policy(branch_undo_policy(kp), kp)),
            "both signs, diamond-opt": [superop(apply_policy(blind, kp)), superop(apply_policy(blind, km))],
            "twirled, no repair": twirl_superop(superop(lp)),
        }
        print(f"\n theta = {theta:.2f}")
        print("  N     " + "".join(f"{name:<25}" for name in channels))
        values = {name: [] for name in channels}
        for N in ROUNDS:
            line = f"  {N:<5} "
            for name, S in channels.items():
                if isinstance(S, list):
                    v = max(1 - ent_fidelity_superop(np.linalg.matrix_power(x, N)) for x in S)
                else:
                    v = 1 - ent_fidelity_superop(np.linalg.matrix_power(S, N))
                values[name].append(v)
                line += f"{v:<25.3e}"
            print(line)
        i10, i100 = ROUNDS.index(10), ROUNDS.index(100)
        growth = values["none"][i100] / values["none"][i10]
        growth_t = values["twirled, no repair"][i100] / values["twirled, no repair"][i10]
        print(f"  growth from N = 10 to 100: no repair x{growth:.1f}; twirled x{growth_t:.1f}; "
              f"known sign, branch undo "
              f"x{values['known sign, branch undo'][i100] / values['known sign, branch undo'][i10]:.1f}")
        # Sign-blind policies computed here: no repair, the twirl, and the
        # diamond-optimal sign-blind repair. ('Fitted, wrong sign' is a guess,
        # not a policy that serves both signs.)
        blind_names = ("none", "twirled, no repair", "both signs, diamond-opt")
        known_names = ("known sign, history", "known sign, branch undo")
        print("  ratio, best sign-blind policy computed / best known-sign policy computed:")
        for i, N in enumerate(ROUNDS):
            best_blind = min(values[n][i] for n in blind_names)
            best_known = min(values[n][i] for n in known_names)
            print(f"    N = {N:<4} {best_blind / best_known:.0f}")
            ok_order &= best_known < best_blind
        ok_growth &= growth > growth_t and abs(growth_t - 10) < 0.5
    print("\n Note: 'both signs, diamond-opt' minimises the worst-sign diamond defect, which is nearly")
    print(" degenerate (no repair is within 3 per cent of it); the solver's choice among near-optimal")
    print(" repairs is not a good controller by infidelity. Part 1 shows that the best sign-blind repair")
    print(" by one-round infidelity is no repair, to within the bracket.")
    passed.append(("accumulation: unrepaired infidelity grows faster than the twirled channel's, which "
                   "grows linearly (10 -> 100 rounds, within 5 per cent)", ok_growth))
    passed.append(("accumulation: best known-sign policy below every sign-blind policy computed, "
                   "at every N", ok_order))


# ---------------------------------------------------------------------------
# Part 3: calibration error
# ---------------------------------------------------------------------------

def report_calibration(inst, passed):
    print("\n=== 3. Calibration error: repairs fitted at theta (1 + eta), channel at theta ===")
    print(" One-round diamond defect [entangled-input lower bound, upper bound]")
    ok, ok_multi = True, True
    for theta in CALIBRATION_THETAS:
        kp, km, classes = class_data(inst, theta)
        lp = list(kp.values())
        none = diamond_norm_of_difference(choi(lp, DIM) - OMEGA, DIM)
        print(f"\n theta = {theta:.2f}; no repair {none:.3e}")
        print("  eta       fitted unitary, history discarded     known-sign optimal, history used     "
              "branch undo")
        results, multi = {}, {}
        for eta in ETAS:
            kq, kqm, cq = class_data(inst, theta * (1 + eta))
            Wq = fitted_unitary_decoder(list(kq.values()), DIM)
            Dw = choi([Wq @ K for K in lp], DIM) - OMEGA
            policy = per_syndrome(diamond_policy([[[kq[s] for s in c]] for c in cq]), cq)
            Dh = choi(apply_policy(policy, kp), DIM) - OMEGA
            undo = branch_undo_policy(kq)
            Du = choi(apply_policy(undo, kp), DIM) - OMEGA
            w = (entangled_input_lower_bound(Dw, DIM), diamond_norm_of_difference(Dw, DIM))
            h = (entangled_input_lower_bound(Dh, DIM), diamond_norm_of_difference(Dh, DIM))
            u = (entangled_input_lower_bound(Du, DIM), diamond_norm_of_difference(Du, DIM))
            results[eta] = h
            S_h, S_u = superop(apply_policy(policy, kp)), superop(apply_policy(undo, kp))
            multi[eta] = [(1 - ent_fidelity_superop(np.linalg.matrix_power(S_h, N)),
                           1 - ent_fidelity_superop(np.linalg.matrix_power(S_u, N)))
                          for N in CALIBRATION_ROUNDS]
            print(f"  {eta:+.3f}    [{w[0]:.3e}, {w[1]:.3e}]              [{h[0]:.3e}, {h[1]:.3e}]"
                  f"             [{u[0]:.3e}, {u[1]:.3e}]")
        for eta in (0.001, 0.01, 0.1):
            print(f"  gap to no repair at eta = {eta:+.3f}: at least {none / results[eta][1]:.0f}-fold")
        ok &= results[0.01][1] < none
        lp_none = superop(lp)
        twirled = twirl_superop(lp_none)
        print("  Entanglement infidelity after N rounds, repair every round "
              "(diamond-optimal history-used / branch undo):")
        print("  eta      " + "".join(f"N = {N:<24}" for N in CALIBRATION_ROUNDS))
        for eta in ETAS:
            print(f"  {eta:+.3f}   " + "".join(f"{a:.3e} / {b:.3e}      " for a, b in multi[eta]))
        tw = [1 - ent_fidelity_superop(np.linalg.matrix_power(twirled, N)) for N in CALIBRATION_ROUNDS]
        print("  twirled, no repair (best sign-blind policy computed): " +
              ", ".join(f"{v:.3e}" for v in tw))
        i300 = CALIBRATION_ROUNDS.index(300)
        change = max(multi[eta][i300][1] / multi[0.0][i300][1] for eta in ETAS)
        print(f"  largest change in 300-round branch-undo infidelity over the eta grid: x{change:.2f}")
        ok_multi &= all(multi[eta][i300][1] < tw[i300] / 10 for eta in ETAS)
    passed.append(("calibration: known-sign repair at 1 per cent error still beats no repair", ok))
    passed.append(("calibration: over 300 rounds, with a calibration error up to 10 per cent, branch-undo "
                   "infidelity stays more than 10 times below the twirl", ok_multi))


# ---------------------------------------------------------------------------
# Part 4: sign error
# ---------------------------------------------------------------------------

SIGN_ERROR_THETAS = (0.05, 0.1)
SIGN_ERROR_ROUNDS = 300
SIGN_ERROR_PROBS = (0.05, 0.01, 0.001)


def sentinel_shots(theta, q_wrong):
    """Smallest odd n for which a majority vote over n sentinel shots gives the
    wrong sign with probability at most q_wrong."""
    from scipy.stats import binom
    p = (1 + np.sin(theta)) / 2
    n = 1
    while binom.cdf((n - 1) // 2, n, p) > q_wrong:
        n += 2
    return n


def report_sign_error(inst, passed):
    print(f"\n=== 4. Sign error: branch undo with the sign estimated once, {SIGN_ERROR_ROUNDS} rounds ===")
    ok = True
    N = SIGN_ERROR_ROUNDS
    for theta in SIGN_ERROR_THETAS:
        kp, km, _ = class_data(inst, theta)
        undo = branch_undo_policy(kp)
        right = 1 - ent_fidelity_superop(np.linalg.matrix_power(superop(apply_policy(undo, kp)), N))
        wrong = 1 - ent_fidelity_superop(np.linalg.matrix_power(superop(apply_policy(undo, km)), N))
        twirl = 1 - ent_fidelity_superop(np.linalg.matrix_power(twirl_superop(superop(list(kp.values()))), N))
        print(f"\n theta = {theta:.2f}: right sign {right:.4e}; wrong sign {wrong:.4e}; twirl {twirl:.4e}")
        print("  P(wrong sign)   expected infidelity   times ideal   twirl / expected   sentinel shots")
        for q in SIGN_ERROR_PROBS:
            e = (1 - q) * right + q * wrong
            print(f"  {q:<14g}  {e:.3e}             {e / right:6.1f}        {twirl / e:7.1f}            "
                  f"{sentinel_shots(theta, q)}")
        q_even = (twirl - right) / (wrong - right)
        q_double = right / (wrong - right)
        print(f"  break-even with the twirl: P(wrong sign) = {q_even:.3f} ({sentinel_shots(theta, q_even)} shots)")
        print(f"  within a factor of 2 of the ideal: P(wrong sign) <= {q_double:.2e} "
              f"({sentinel_shots(theta, q_double)} shots)")
        ok &= 0 < q_even < 0.5 and q_double < q_even
    passed.append(("sign error: break-even with the twirl at a wrong-sign probability below one half, "
                   "and the factor-of-2 threshold stricter", ok))


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
    print(f"fidelity thetas = {FIDELITY_THETAS}; accumulation thetas = {ACCUMULATION_THETAS}; "
          f"rounds = {ROUNDS}")
    print(f"calibration thetas = {CALIBRATION_THETAS}; eta = {ETAS}")
    print("diamond norm unnormalised, range [0, 2]; infidelity = 1 - entanglement fidelity")

    passed = []
    inst = build_instrument(3)
    t0 = time.time()
    report_infidelity(inst, passed)
    report_accumulation(inst, passed)
    report_calibration(inst, passed)
    report_sign_error(inst, passed)

    print("\n=== Summary ===")
    for label, ok in passed:
        print(f"{'PASS' if ok else 'FAIL'}  {label}")
    statuses = ", ".join(sorted(set(SOLVER_FLAGS))) or "none"
    print(f"SDP solves not reported 'optimal' by the solver: {len(SOLVER_FLAGS)} ({statuses}). "
          "Brackets and defects are constructed independently of solver tolerance; "
          "floating-point estimates, not interval-certified.")
    print(f"(elapsed {time.time() - t0:.0f} s)")


if __name__ == "__main__":
    main()
