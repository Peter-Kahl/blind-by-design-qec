#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
two_round_repair_check.py
=========================

Version:
    0.7.0 (2026-10-31)

Supplementary research code for:

    Peter Kahl, 'Sign-Blind Recovery, Record Leakage and Model
    Misspecification in Quantum Error Correction: Worked Numerical Cases'
    (2026), §3 (two rounds); and

    Peter Kahl, 'Blind by Design: Quantum Error Correction, Functional
    Incompatibility and the Value of Evidence' (2026), §5.3.1, extending
    Table 7 from one round to two.

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
Section 4 of prop52_numerical_check.py compares, for one round of the L = 3
toric instrument of Özgüler (2026), the best repair over all channels for a
known rotation sign and for both signs, with the syndrome history discarded
and with it used. The complete syndrome history has the same law at +theta
and -theta (Özgüler 2026, Proposition 1), so a longer history cannot reveal
the sign; but it might still help a repair that must serve both signs. This
script repeats the comparison over two rounds.

One repair channel is applied after both rounds, chosen either without the
syndrome history or according to the complete two-round history. Repairs
applied between rounds are not considered.

Method
------
* One-round syndromes are grouped into classes whose branch maps agree up to
  weight at both signs (syndrome_classes in prop52_numerical_check.py). Each
  class is replaced by one Kraus operator per sign, scaled so that its map
  equals the sum of the members' maps. The largest deviation is reported.
* Two-round histories are formed from ordered pairs of one-round classes (the
  second round applied after the first) and grouped again in the same way.
  Every one of the 256^2 two-round histories belongs to exactly one pair.
* Sign known, history used: one repair per two-round class.
* Both signs, history used: classes whose probability is below POOL share a
  single repair; the others have one each. Without this pooling the solver
  fails numerically on the full problem (probabilities span about 1 to
  1e-16). The pooled policy is a valid policy, so its achieved defect is a
  valid upper bound on the optimum over policies with one repair per class.
  The lower bound is computed from the same dual point with the dual
  constraint made feasible for every class individually (the all_branches
  argument of optimal_repair), so it holds for that optimum too and does not
  depend on the pooling.
* For the sign-known repair, the lower bound reported is the smaller of the
  grouped and the class-by-class bound.
* Upper bounds are the defects achieved by the solver's repairs, made exact
  channels and evaluated with diamond_norm_of_difference; lower bounds are
  dual objectives at dual points made exactly feasible. Both are
  floating-point constructions, not interval-certified, as in
  prop52_numerical_check.py.

Known limitation
----------------
At theta = 0.05 the sign-known, history-used bracket is wide (its lower end
is about 40 times below its upper end): the defects are of order 1e-6 and
the solver reports its solutions as 'optimal_inaccurate'. Only the upper end
of that bracket is informative there, and the summary checks use it.

Requirements
------------
Python 3.9 or later; NumPy 1.20 or later; CVXPY 1.4 or later with the
Clarabel solver. Imports prop52_numerical_check.py and ozguler_prop1_check.py
from the same folder.

Run
---
    python two_round_repair_check.py             # theta = 0.05 and 0.1
    python two_round_repair_check.py 0.1         # chosen angles

About 4 to 8 minutes per angle on the reference machine.
"""

from __future__ import annotations

import platform
import sys
import time

import cvxpy as cp
import numpy as np

from ozguler_prop1_check import build_instrument, kraus_operators
from prop52_numerical_check import (SOLVER, SOLVER_FLAGS, choi,
                                    diamond_norm_of_difference,
                                    fitted_unitary_decoder, identity_choi,
                                    optimal_repair, syndrome_classes)

SCRIPT_VERSION = "0.7.0"

DEFAULT_THETAS = (0.05, 0.1)
POOL = 1e-8          # classes below this probability share one sign-blind repair
TOLERANCE = 1e-12    # rounding tolerance for the structural checks


# ---------------------------------------------------------------------------
# Class representatives and two-round histories
# ---------------------------------------------------------------------------

def class_representatives(kp, km, classes):
    """One Kraus operator per class and sign, scaled so that its map equals the
    sum of the class members' maps."""
    reps = []
    for c in classes:
        rp, rm = kp[c[0]], km[c[0]]
        wp = sum(np.linalg.norm(kp[s]) ** 2 for s in c) / np.linalg.norm(rp) ** 2
        wm = sum(np.linalg.norm(km[s]) ** 2 for s in c) / np.linalg.norm(rm) ** 2
        reps.append((np.sqrt(wp) * rp, np.sqrt(wm) * rm))
    return reps


def class_map_residual(kp, km, classes, reps):
    """Largest deviation of a representative's map from its class's summed map."""
    worst = 0.0
    for c, (ap, am) in zip(classes, reps):
        for K, a in ((kp, ap), (km, am)):
            S = sum(np.kron(K[s], K[s].conj()) for s in c)
            worst = max(worst, float(np.abs(S - np.kron(a, a.conj())).max()))
    return worst


def two_round_classes(inst, theta):
    kp, km = kraus_operators(inst, theta), kraus_operators(inst, -theta)
    cls1 = syndrome_classes(kp, km)
    reps1 = class_representatives(kp, km, cls1)
    res1 = class_map_residual(kp, km, cls1, reps1)
    n1 = len(cls1)
    kp2 = {(a, b): reps1[b][0] @ reps1[a][0] for a in range(n1) for b in range(n1)}
    km2 = {(a, b): reps1[b][1] @ reps1[a][1] for a in range(n1) for b in range(n1)}
    cls2 = syndrome_classes(kp2, km2)
    reps2 = class_representatives(kp2, km2, cls2)
    res2 = class_map_residual(kp2, km2, cls2, reps2)
    covered = sum(len(cls1[a]) * len(cls1[b]) for c in cls2 for a, b in c)
    return n1, reps2, max(res1, res2), covered


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------

def run(theta, inst, passed):
    n1, reps, residual, covered = two_round_classes(inst, theta)
    lp, lm = [r[0] for r in reps], [r[1] for r in reps]
    completeness = max(float(np.abs(sum(K.conj().T @ K for K in L) - np.eye(4)).max())
                       for L in (lp, lm))
    weights = [float(np.linalg.norm(a) ** 2) / 4 for a in lp]

    print(f"\n=== theta = {theta:.2f}, two rounds, one repair after both rounds ===")
    print(f" one-round classes {n1}; two-round classes {len(reps)}, covering {covered} "
          f"of {256 ** 2} histories")
    print(f" largest class-map residual {residual:.1e}; completeness residual {completeness:.1e}")

    I = identity_choi(4)
    W = fitted_unitary_decoder(lp, 4)
    right = diamond_norm_of_difference(choi([W @ K for K in lp], 4) - I, 4)
    wrong = diamond_norm_of_difference(choi([W @ K for K in lm], 4) - I, 4)
    none_p = diamond_norm_of_difference(choi(lp, 4) - I, 4)
    none_m = diamond_norm_of_difference(choi(lm, 4) - I, 4)
    print(" Fitted unitary repair at the right and the wrong sign, and no repair (upper bounds):")
    print(f"  right sign {right:.3e}   wrong sign {wrong:.3e}   no repair {none_p:.3e} (+theta), "
          f"{none_m:.3e} (-theta)")

    t0 = time.time()
    p_lo, p_up = optimal_repair([[lp]], 4)
    u_lo, u_up = optimal_repair([[lp, lm]], 4)
    ph = optimal_repair([[[a]] for a in lp], 4, [[[a]] for a in lp])
    ph_lo, ph_up = min(ph[0], ph[2]), ph[1]
    print(f" (sign-known optima done, {time.time() - t0:.0f} s)", flush=True)

    kept = [r for r, w in zip(reps, weights) if w >= POOL]
    pooled = [r for r, w in zip(reps, weights) if w < POOL]
    groups = [[[a], [b]] for a, b in kept]
    if pooled:
        groups.append([[a for a, _ in pooled], [b for _, b in pooled]])
    _, uh_up, uh_lo = optimal_repair(groups, 4, [[[a], [b]] for a, b in reps])

    print(" Best repair over all channels, as [lower bound, achieved upper bound]:")
    print("  sign known,  history discarded  "
          f"[{p_lo:.3e}, {p_up:.3e}]")
    print("  both signs,  history discarded  "
          f"[{u_lo:.3e}, {u_up:.3e}]")
    print("  sign known,  history used       "
          f"[{ph_lo:.3e}, {ph_up:.3e}]   ({len(reps)} repairs)")
    print("  both signs,  history used       "
          f"[{uh_lo:.3e}, {uh_up:.3e}]   ({len(kept)} repairs, plus one shared by "
          f"{len(pooled)} classes of total probability {sum(w for w in weights if w < POOL):.1e})")
    print(f"  ({time.time() - t0:.0f} s)")
    gain = 1 - uh_lo / none_p
    print(f" Sign-blind repair: using the history lowers the defect by at most {100 * gain:.1f} per cent "
          "relative to no repair")
    print(f" Conservative ratio, sign-blind lower bound / sign-known upper bound (history used): "
          f"{uh_lo / ph_up:.0f}")

    tag = f"theta = {theta:.2f}"
    passed.append((f"{tag}: class maps and completeness exact to rounding",
                   residual < TOLERANCE and completeness < TOLERANCE and covered == 256 ** 2))
    passed.append((f"{tag}: sign-blind repair, history used, more than twice the known-sign optimum",
                   uh_lo > 2 * ph_up))
    passed.append((f"{tag}: using the history helps the sign-blind repair less than the known-sign one",
                   uh_up <= u_up * (1 + 1e-6) and uh_lo / ph_up >= u_up / p_lo))
    passed.append((f"{tag}: no repair within 3 per cent of the sign-blind lower bound, history discarded",
                   none_p <= 1.03 * u_lo))


def main() -> None:
    thetas = tuple(float(x) for x in sys.argv[1:]) or DEFAULT_THETAS
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
    print(f"thetas = {thetas}; pooling threshold for sign-blind repairs = {POOL:g}")
    print("diamond norm unnormalised, range [0, 2]")

    passed = []
    inst = build_instrument(3)
    for theta in thetas:
        run(theta, inst, passed)

    print("\n=== Summary ===")
    for label, ok in passed:
        print(f"{'PASS' if ok else 'FAIL'}  {label}")
    statuses = ", ".join(sorted(set(SOLVER_FLAGS))) or "none"
    print(f"SDP solves not reported 'optimal' by the solver: {len(SOLVER_FLAGS)} ({statuses}). "
          "Reported brackets are constructed independently of solver tolerance; both ends are "
          "floating-point estimates, not interval-certified.")


if __name__ == "__main__":
    main()
