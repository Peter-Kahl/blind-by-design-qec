#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
interval_certification_check.py
===============================

Version:
    0.12.0 (2026-10-31)

Supplementary research code for:

    Peter Kahl, 'The Cost of Sign-Blind Recovery under Coherent Noise: Worked Numerical Cases'
    (2026), §§2-4 (certified bounds behind the headline ratios).

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
The other scripts report primal-dual brackets constructed in double
precision. This script certifies, by interval arithmetic, the four bounds
behind the companion's headline ratios at theta = 0.05 and 0.1, for one
round of the L = 3 toric instrument of Özgüler (2026):

1. an upper bound on the diamond-norm defect of branch undo (the
   sign-known, history-used repair);
2. a lower bound on the diamond-norm defect of every sign-blind repair
   that may depend on the syndrome (one repair per syndrome, worst sign);
3. the entanglement infidelity of branch undo, enclosed exactly; and
4. a lower bound on the worst-sign entanglement infidelity of every
   sign-blind, syndrome-dependent repair.

Method
------
Exact branch operators. With c = cos(theta/2) and s = sin(theta/2), one
round applies prod_q (c I - i s X_q) to the code state. Each X-error
pattern E with syndrome sigma is followed by the deterministic recovery
r(sigma); E + r(sigma) has trivial syndrome, so on the code space it acts as
one of the four logical X operators X_L, L in {00, 01, 10, 11} (the star
stabilisers act as the identity there). Hence

    K_sigma = sum_L a_{sigma,L} X_L,
    a_{sigma,L} = sum_w N[sigma, L, w] c^(18 - w) (-i s)^w,

where N[sigma, L, w] counts the error patterns of weight w with syndrome
sigma whose corrected pattern lies in class L. The integer counts are
computed once by enumerating all 2^18 patterns. The script checks that this
formula reproduces kraus_operators() of ozguler_prop1_check.py in double
precision.

Logical-X basis. All X_L are diagonal in the product Hadamard basis, so
every K_sigma is diagonal there, with eigenvalues
lambda_{sigma,k} = sum_L a_{sigma,L} (-1)^(k.L). The repair problems are
unitarily equivalent in this basis (the same unitary on input and output
leaves every defect and fidelity unchanged), and all certificates are built
in it. The branch at -theta has eigenvalues conj(lambda).

Intervals. c and s are enclosed with mpmath interval arithmetic
(iv.dps = 40) for theta given as the exact decimals 0.05 and 0.1; every
lambda, matrix entry and objective value is an interval enclosing the true
value. Positive definiteness is proved by interval Cholesky factorisation
of the real symmetric 2n x 2n form [[Re, -Im], [Im, Re]] of each Hermitian
matrix: if the factorisation completes with every pivot interval strictly
positive, every symmetric matrix in the interval box is positive definite.

Certificates (D_sigma ranges over channels, one per syndrome).

1. Branch undo maps each branch to diag(|lambda_{sigma,k}|), so the repaired
   channel is the Schur multiplier rho_kl -> M_kl rho_kl with
   M_kl = sum_sigma |lambda_{sigma,k}| |lambda_{sigma,l}|, and the Choi
   matrix of (repaired - identity) is sum_kl (M_kl - 1)|kk><ll|. By Watrous
   (2013), for a difference of channels, ||Phi||_diamond = 2 min t over
   Z >= 0, Z >= J, Tr_out Z <= t I. Any Z supported on span{|kk>} with
   4 x 4 block Zt >= 0 and Zt >= M - 1 is feasible, with Tr_out Z = diag(Zt);
   so 2 max_k Zt_kk is an upper bound. Zt comes from a floating-point
   solve, is shifted by a small multiple of I, and its feasibility is
   proved by interval Cholesky.
2. The dual of min over policies of max over signs of the defect (as in
   optimal_repair() of prop52_numerical_check.py): for W_pm >= 0 and
   rho_pm >= 0 with rho_pm (x) I >= W_pm, and Y_sigma with
   H_sigma - Y_sigma (x) I >= 0, where
   H_sigma = sum_pm (conj(K_sigma^pm) (x) I) W_pm (conj(K_sigma^pm) (x) I)^dagger,
   the optimum is at least
   2 (sum_sigma tr Y_sigma - sum_pm tr(W_pm Omega)) / sum_pm tr rho_pm.
   W and rho come from a floating-point solve over syndrome classes, are
   shifted by small multiples of I, and Y_sigma is solved for each of the
   256 syndromes and shifted down; every constraint is proved by interval
   Cholesky, syndrome by syndrome.
3. Branch-undo fidelity: F = sum_sigma (sum_k |lambda_{sigma,k}|)^2 / 16,
   evaluated in intervals.
4. For any policy, the worse sign's fidelity is at most the average, and
   the average is at most sum_sigma max over channels of tr(J H_sigma) with
   H_sigma = sum_pm (1/2) (conj(K^pm) (x) I) (Omega / 16) (conj(K^pm) (x) I)^dagger,
   which is at most tr Y_sigma for any Y_sigma with Y_sigma (x) I >= H_sigma.
   Y_sigma is solved for each syndrome, shifted up, and proved feasible.

What is certified. The bounds hold for the exact instrument at the exact
angles, for repairs that may depend on the syndrome of the round (one
round). They do not depend on the solver's accuracy. The instrument's
definition (recovery table, code states, logical operators) is taken from
ozguler_prop1_check.py.

Two stages
----------
Propose (needs CVXPY and Clarabel): a floating-point solver proposes the
certificates, which are shifted until interval arithmetic verifies them,
and saved to output/interval_certificates.json.

Verify (needs only NumPy and mpmath): the saved certificates are loaded and
every bound is proved from them in interval arithmetic. The derivations are
given in Appendix B of the companion note.

Run
---
    python interval_certification_check.py            # propose, save, verify
    python interval_certification_check.py --verify   # verify saved certificates only

Requires mpmath (1.3 or later). The propose stage also needs the packages of
prop52_numerical_check.py. About 3 minutes for both stages; the verify stage
alone takes about 2 minutes.
"""

from __future__ import annotations

import argparse
import json
import os
import platform
import time
import warnings

import mpmath
import numpy as np
from mpmath import iv

from ozguler_prop1_check import build_instrument, kraus_operators, mask_of

SCRIPT_VERSION = "0.12.0"
THETAS = ("0.05", "0.1")
DIM = 4
N2 = DIM * DIM
iv.dps = 40
SHIFT = 1e-12          # relative shift used to make float certificates strictly feasible
CERT_PATH = os.path.join("output", "interval_certificates.json")
SOLVER = "CLARABEL"


# ---------------------------------------------------------------------------
# Exact branch structure
# ---------------------------------------------------------------------------

def branch_counts(inst):
    """N[sigma, L, w]: number of X-error patterns of weight w with syndrome
    sigma whose corrected pattern acts on the code space as X_L."""
    n = inst["n"]
    syn = inst["syndrome"]
    rec = np.zeros(int(syn.max()) + 1, dtype=np.int64)
    for s, r in inst["recovery"].items():
        rec[s] = r
    star = {0}
    for row in inst["hx"]:
        m = mask_of(row)
        star |= {v ^ m for v in star}
    lx1, lx2 = inst["logical_x1"], inst["logical_x2"]
    class_of = np.full(1 << n, -1, dtype=np.int64)
    for a in (0, 1):
        for b in (0, 1):
            shift = (lx1 if a else 0) ^ (lx2 if b else 0)
            for g in star:
                class_of[g ^ shift] = 2 * a + b
    E = np.arange(1 << n, dtype=np.int64)
    M = E ^ rec[syn]
    L = class_of[M]
    assert (L >= 0).all(), "a corrected pattern is not syndrome-free"
    w = np.bitwise_count(E.astype(np.uint64)).astype(np.int64)
    N = np.zeros((len(rec), 4, n + 1), dtype=np.int64)
    np.add.at(N, (syn, L, w), 1)
    return N


def float_branches(N, theta, n=18):
    c, s = np.cos(theta / 2), np.sin(theta / 2)
    w = np.arange(n + 1)
    coeff = c ** (n - w) * (-1j * s) ** w
    a = N @ coeff                                    # (sigma, L)
    X = np.array([[0, 1], [1, 0]])
    XL = [np.kron(np.linalg.matrix_power(X, L >> 1), np.linalg.matrix_power(X, L & 1)) for L in range(4)]
    return {sig: sum(a[sig, L] * XL[L] for L in range(4)) for sig in range(a.shape[0])}


def characters():
    return [[(-1) ** (((k >> 1) & (L >> 1)) ^ ((k & 1) & (L & 1))) for L in range(4)] for k in range(4)]


def float_eigs(N, theta, n=18):
    c, s = np.cos(theta / 2), np.sin(theta / 2)
    w = np.arange(n + 1)
    a = N @ (c ** (n - w) * (-1j * s) ** w)
    chi = np.array(characters())
    return a @ chi.T                                 # (sigma, k)


def interval_eigs(N, theta_str, n=18):
    """Interval enclosures (re, im) of lambda_{sigma,k}."""
    th = iv.mpf(theta_str)
    c, s = iv.cos(th / 2), iv.sin(th / 2)
    cp_ = [c ** (n - w) for w in range(n + 1)]
    sp_ = [s ** w for w in range(n + 1)]
    chi = characters()
    out = []
    for sig in range(N.shape[0]):
        a_re, a_im = [], []
        for L in range(4):
            re, im = iv.mpf(0), iv.mpf(0)
            for w in range(n + 1):
                cnt = int(N[sig, L, w])
                if cnt == 0:
                    continue
                term = cnt * cp_[w] * sp_[w]
                r = w % 4                            # (-i)^w = 1, -i, -1, i
                if r == 0:
                    re += term
                elif r == 1:
                    im -= term
                elif r == 2:
                    re -= term
                else:
                    im += term
            a_re.append(re)
            a_im.append(im)
        row = []
        for k in range(4):
            re = sum((chi[k][L] * a_re[L] for L in range(4)), iv.mpf(0))
            im = sum((chi[k][L] * a_im[L] for L in range(4)), iv.mpf(0))
            row.append((re, im))
        out.append(row)
    return out


# ---------------------------------------------------------------------------
# Interval linear algebra
# ---------------------------------------------------------------------------

def lo(x, digits=8):
    return mpmath.nstr(mpmath.mpf(x._mpi_[0]), digits)


def hi(x, digits=8):
    return mpmath.nstr(mpmath.mpf(x._mpi_[1]), digits)


def iv_abs(z):
    return iv.sqrt(z[0] ** 2 + z[1] ** 2)


def to_iv(x):
    return iv.mpf(float(x))


def pd_certified(re, im):
    """re, im: n x n lists of intervals (re symmetric, im antisymmetric), the
    real and imaginary parts of a Hermitian interval matrix. True if interval
    Cholesky of [[re, -im], [im, re]] completes with positive pivots."""
    n = len(re)
    m = 2 * n

    def entry(i, j):
        bi, ii = divmod(i, n)
        bj, jj = divmod(j, n)
        if bi == bj:
            return re[ii][jj]
        return -im[ii][jj] if bi == 0 else im[ii][jj]

    Lm = [[None] * m for _ in range(m)]
    for j in range(m):
        acc = entry(j, j)
        for k in range(j):
            acc -= Lm[j][k] ** 2
        if not acc.a > 0:
            return False
        Lm[j][j] = iv.sqrt(acc)
        for i in range(j + 1, m):
            acc = entry(i, j)
            for k in range(j):
                acc -= Lm[i][k] * Lm[j][k]
            Lm[i][j] = acc / Lm[j][j]
    return True


def hermitian_parts(A):
    """Split an n x n matrix of complex pairs (re, im) into symmetric re and
    antisymmetric im, imposing exact Hermitian structure on the diagonal."""
    n = len(A)
    re = [[A[i][j][0] if i <= j else A[j][i][0] for j in range(n)] for i in range(n)]
    im = [[(A[i][j][1] if i < j else (-A[j][i][1] if i > j else iv.mpf(0))) for j in range(n)] for i in range(n)]
    return re, im


def point_matrix(X):
    return [[(to_iv(X[i, j].real), to_iv(X[i, j].imag)) for j in range(X.shape[1])] for i in range(X.shape[0])]


def sub(A, B):
    return [[(A[i][j][0] - B[i][j][0], A[i][j][1] - B[i][j][1]) for j in range(len(A))] for i in range(len(A))]


def kron_eye(Y):
    """Y (x) I_4 for a 4 x 4 interval matrix Y."""
    out = [[(iv.mpf(0), iv.mpf(0)) for _ in range(N2)] for _ in range(N2)]
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                out[i * DIM + k][j * DIM + k] = Y[i][j]
    return out


def conj_scale(W, lam_list, weights):
    """sum_pm weight * (diag(conj lam) (x) I) W_pm (diag(conj lam) (x) I)^dagger,
    entrywise: H[i, j] = sum_pm w conj(lam[k(i)]) W[i, j] lam[k(j)], k(i) = i // 4."""
    H = [[(iv.mpf(0), iv.mpf(0)) for _ in range(N2)] for _ in range(N2)]
    for Wpm, lam, wt in zip(W, lam_list, weights):
        for i in range(N2):
            ar, ai = lam[i // DIM]
            for j in range(N2):
                br, bi = lam[j // DIM]
                # conj(a) * b
                pr, pi = ar * br + ai * bi, ar * bi - ai * br
                wr, wi = Wpm[i][j]
                H[i][j] = (H[i][j][0] + wt * (pr * wr - pi * wi), H[i][j][1] + wt * (pr * wi + pi * wr))
    return H


# ---------------------------------------------------------------------------
# Serialisation of certificates (exact binary floats, round-tripped by JSON)
# ---------------------------------------------------------------------------

def enc(X):
    X = np.asarray(X, dtype=complex)
    return {"re": X.real.tolist(), "im": X.imag.tolist()}


def dec(d):
    return np.array(d["re"], dtype=float) + 1j * np.array(d["im"], dtype=float)


# ---------------------------------------------------------------------------
# Interval building blocks shared by both stages
# ---------------------------------------------------------------------------

def interval_data(theta_str, N, sigmas):
    lam = interval_eigs(N, theta_str)
    absl = [[iv_abs(lam[s][k]) for k in range(DIM)] for s in range(len(lam))]
    M = [[sum((absl[s][k] * absl[s][l] for s in sigmas), iv.mpf(0)) for l in range(DIM)] for k in range(DIM)]
    Wm = [[M[k][l] - 1 for l in range(DIM)] for k in range(DIM)]
    return lam, absl, Wm


def lam_pair(lam, s):
    return [lam[s], [(re, -im) for re, im in lam[s]]]


def omega16():
    return [[((iv.mpf(1) / 16 if (i % (DIM + 1) == 0 and j % (DIM + 1) == 0) else iv.mpf(0)), iv.mpf(0))
             for j in range(N2)] for i in range(N2)]


def check_upper(Zt, Wm):
    Zi = [[to_iv(Zt[k, l]) for l in range(DIM)] for k in range(DIM)]
    zero = [[iv.mpf(0)] * DIM for _ in range(DIM)]
    return pd_certified(Zi, zero) and pd_certified([[Zi[k][l] - Wm[k][l] for l in range(DIM)]
                                                    for k in range(DIM)], zero)


def check_W_rho(Wc, Rc):
    Wi = [point_matrix(W) for W in Wc]
    Ri = [point_matrix(R) for R in Rc]
    ok = all(pd_certified(*hermitian_parts(W)) for W in Wi)
    ok &= all(pd_certified(*hermitian_parts(sub(kron_eye(R), W))) for R, W in zip(Ri, Wi))
    return ok, Wi, Ri


def check_Y_below(Hi, Y):
    return pd_certified(*hermitian_parts(sub(Hi, kron_eye(point_matrix(Y)))))


def check_Y_above(Hi, Y):
    return pd_certified(*hermitian_parts(sub(kron_eye(point_matrix(Y)), Hi)))


def trace_iv(Y):
    return sum((to_iv(np.real(Y[k, k])) for k in range(DIM)), iv.mpf(0))


# ---------------------------------------------------------------------------
# Stage 1: propose certificates (needs CVXPY)
# ---------------------------------------------------------------------------

def propose(theta_str, N, inst):
    import cvxpy as cp
    from prop52_numerical_check import choi_map, syndrome_classes

    def solve(problem):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            problem.solve(solver=SOLVER)

    def float_H(W, lams, weights):
        H = np.zeros((N2, N2), dtype=complex)
        for Wpm, lam, wt in zip(W, lams, weights):
            A = np.kron(np.diag(np.conj(lam)), np.eye(DIM))
            H += wt * A @ Wpm @ A.conj().T
        return (H + H.conj().T) / 2

    def y_below(H):
        Y = cp.Variable((DIM, DIM), hermitian=True)
        solve(cp.Problem(cp.Maximize(cp.real(cp.trace(Y))), [H - cp.kron(Y, np.eye(DIM)) >> 0]))
        return (Y.value + Y.value.conj().T) / 2

    def y_above(H):
        Y = cp.Variable((DIM, DIM), hermitian=True)
        solve(cp.Problem(cp.Minimize(cp.real(cp.trace(Y))), [cp.kron(Y, np.eye(DIM)) - H >> 0]))
        return (Y.value + Y.value.conj().T) / 2

    theta = float(theta_str)
    sigmas = sorted(inst["recovery"])
    lam_f = float_eigs(N, theta)
    lam, absl, Wm = interval_data(theta_str, N, sigmas)
    cert = {"theta": theta_str}

    # B.1: upper bound for branch undo
    Wm_f = np.array([[float(mpmath.mpf(Wm[k][l].mid)) for l in range(DIM)] for k in range(DIM)])
    Zt = cp.Variable((DIM, DIM), symmetric=True)
    tt = cp.Variable()
    solve(cp.Problem(cp.Minimize(tt), [Zt >> 0, Zt - Wm_f >> 0, cp.diag(Zt) <= tt]))
    Ztv = (Zt.value + Zt.value.T) / 2
    scale = max(1e-30, float(np.abs(Wm_f).max()))
    for attempt in range(10):
        shift = scale * SHIFT * 10 ** attempt + max(0.0, -np.linalg.eigvalsh(Ztv - Wm_f).min(),
                                                    -np.linalg.eigvalsh(Ztv).min())
        Zc = Ztv + shift * np.eye(DIM)
        if check_upper(Zc, Wm):
            break
    cert["Z"] = enc(Zc)

    # B.2: dual point for the sign-blind diamond problem, over syndrome classes
    kp = {s: np.diag(lam_f[s]) for s in sigmas}
    km = {s: np.diag(np.conj(lam_f[s])) for s in sigmas}
    classes = syndrome_classes(kp, km)
    Om = np.outer(np.eye(DIM).reshape(-1), np.eye(DIM).reshape(-1)).astype(complex)
    groups = [[[kp[s] for s in c], [km[s] for s in c]] for c in classes]
    maps = [[choi_map(g[sgn], DIM) for sgn in range(2)] for g in groups]
    JD = [cp.Variable((N2, N2), hermitian=True) for _ in groups]
    t = cp.Variable()
    cons, w_cons, r_cons = [], [], []
    for J in JD:
        cons += [J >> 0, cp.partial_trace(J, [DIM, DIM], axis=1) == np.eye(DIM)]
    for sgn in range(2):
        out = sum(cp.reshape(T[sgn] @ cp.vec(J, order='C'), (N2, N2), order='C') for T, J in zip(maps, JD))
        Z = cp.Variable((N2, N2), hermitian=True)
        c_w = Z - (out - Om) >> 0
        c_r = t * np.eye(DIM) - cp.partial_trace(Z, [DIM, DIM], axis=1) >> 0
        cons += [Z >> 0, c_w, c_r]
        w_cons.append(c_w)
        r_cons.append(c_r)
    solve(cp.Problem(cp.Minimize(t), cons))
    cert["float_primal_sign_blind"] = 2 * float(t.value)
    W, R = [], []
    for c_w, c_r in zip(w_cons, r_cons):
        Wv = (c_w.dual_value + c_w.dual_value.conj().T) / 2
        e, v = np.linalg.eigh(Wv)
        W.append((v * np.clip(e, 0, None)) @ v.conj().T)
        Rv = (c_r.dual_value + c_r.dual_value.conj().T) / 2
        e, v = np.linalg.eigh(Rv)
        R.append((v * np.clip(e, 0, None)) @ v.conj().T)
    a = SHIFT * max(np.abs(W[0]).max(), np.abs(W[1]).max())
    Wc = [Wp + a * np.eye(N2) for Wp in W]
    Rc = [Rp + (max(0.0, -np.linalg.eigvalsh(np.kron(Rp, np.eye(DIM)) - Wp).min()) + SHIFT * np.abs(Rp).max() + a)
          * np.eye(DIM) for Wp, Rp in zip(Wc, R)]
    cert["W"] = [enc(Wp) for Wp in Wc]
    cert["rho"] = [enc(Rp) for Rp in Rc]
    _, Wi, _ = check_W_rho(Wc, Rc)
    Yd = {}
    for s in sigmas:
        Hi = conj_scale(Wi, lam_pair(lam, s), [iv.mpf(1), iv.mpf(1)])
        Hf = float_H(Wc, [lam_f[s], np.conj(lam_f[s])], [1.0, 1.0])
        Yf = y_below(Hf)
        for attempt in range(10):
            c = max(0.0, -np.linalg.eigvalsh(Hf - np.kron(Yf, np.eye(DIM))).min()) + \
                SHIFT * 10 ** attempt * max(1e-30, np.abs(Hf).max())
            Yc = Yf - c * np.eye(DIM)
            if check_Y_below(Hi, Yc):
                break
        Yd[str(s)] = enc(Yc)
    cert["Y_diamond"] = Yd

    # B.3: infidelity dual
    Om16i = omega16()
    Om16f = Om / 16
    Yi_ = {}
    for s in sigmas:
        Hi = conj_scale([Om16i, Om16i], lam_pair(lam, s), [iv.mpf(1) / 2, iv.mpf(1) / 2])
        Hf = float_H([Om16f, Om16f], [lam_f[s], np.conj(lam_f[s])], [0.5, 0.5])
        Yf = y_above(Hf)
        for attempt in range(10):
            c = max(0.0, -np.linalg.eigvalsh(np.kron(Yf, np.eye(DIM)) - Hf).min()) + \
                SHIFT * 10 ** attempt * max(1e-30, np.abs(Hf).max())
            Yc = Yf + c * np.eye(DIM)
            if check_Y_above(Hi, Yc):
                break
        Yi_[str(s)] = enc(Yc)
    cert["Y_infidelity"] = Yi_
    return cert


# ---------------------------------------------------------------------------
# Stage 2: verify certificates (NumPy and mpmath only)
# ---------------------------------------------------------------------------

def verify(cert, N, inst, passed):
    theta_str = cert["theta"]
    theta = float(theta_str)
    t0 = time.time()
    print(f"\n=== theta = {theta_str} ===")
    Kf = kraus_operators(inst, theta)
    Kx = float_branches(N, theta)
    gap = max(np.abs(Kf[s] - Kx[s]).max() for s in Kf)
    print(f" exact branch formula vs state-vector instrument: max difference {gap:.1e}")
    passed.append((f"theta = {theta_str}: branch formula reproduces the instrument (to 1e-12)", gap < 1e-12))

    sigmas = sorted(Kf)
    lam, absl, Wm = interval_data(theta_str, N, sigmas)

    # B.4 and no repair
    F_undo = sum((sum(absl[s], iv.mpf(0)) ** 2 for s in sigmas), iv.mpf(0)) / 16
    F_none = sum(((sum((lam[s][k][0] for k in range(DIM)), iv.mpf(0))) ** 2 +
                  (sum((lam[s][k][1] for k in range(DIM)), iv.mpf(0))) ** 2 for s in sigmas), iv.mpf(0)) / 16
    inf_undo, inf_none = 1 - F_undo, 1 - F_none
    print(f" infidelity, no repair     in [{lo(inf_none, 12)}, {hi(inf_none, 12)}]")
    print(f" infidelity, branch undo   in [{lo(inf_undo, 12)}, {hi(inf_undo, 12)}]")

    # B.1
    Zc = dec(cert["Z"]).real
    ok_up = check_upper(Zc, Wm)
    d_undo_upper = iv.mpf(2) * max(to_iv(Zc[k, k]) for k in range(DIM))
    print(f" diamond defect, branch undo        <= {hi(d_undo_upper)}   (certificate "
          f"{'verified' if ok_up else 'NOT verified'})")
    passed.append((f"theta = {theta_str}: branch-undo diamond upper bound certified", ok_up))

    # B.2
    Wc = [dec(w) for w in cert["W"]]
    Rc = [dec(r) for r in cert["rho"]]
    ok, Wi, Ri = check_W_rho(Wc, Rc)
    total_Y = iv.mpf(0)
    for s in sigmas:
        Y = dec(cert["Y_diamond"][str(s)])
        Hi = conj_scale(Wi, lam_pair(lam, s), [iv.mpf(1), iv.mpf(1)])
        ok &= check_Y_below(Hi, Y)
        total_Y += trace_iv(Y)
    Om_tr = [sum((Wp[k * DIM + k][l * DIM + l][0] for k in range(DIM) for l in range(DIM)), iv.mpf(0)) for Wp in Wi]
    norm = sum((sum((Rp[k][k][0] for k in range(DIM)), iv.mpf(0)) for Rp in Ri), iv.mpf(0))
    ok &= norm.a > 0
    d_blind_lower = 2 * (total_Y - Om_tr[0] - Om_tr[1]) / norm
    print(f" diamond defect, sign-blind (any syndrome-dependent repair) >= {lo(d_blind_lower)}"
          f"   (float primal {cert.get('float_primal_sign_blind', float('nan')):.6e}; certificate "
          f"{'verified' if ok else 'NOT verified'})")
    passed.append((f"theta = {theta_str}: sign-blind diamond lower bound certified for all 256 syndromes", ok))

    # B.3
    Om16i = omega16()
    okf = True
    total_Yf = iv.mpf(0)
    for s in sigmas:
        Y = dec(cert["Y_infidelity"][str(s)])
        Hi = conj_scale([Om16i, Om16i], lam_pair(lam, s), [iv.mpf(1) / 2, iv.mpf(1) / 2])
        okf &= check_Y_above(Hi, Y)
        total_Yf += trace_iv(Y)
    inf_blind_lower = 1 - total_Yf
    print(f" infidelity, sign-blind (any syndrome-dependent repair, worse sign) >= "
          f"{lo(inf_blind_lower)}   (certificate {'verified' if okf else 'NOT verified'})")
    passed.append((f"theta = {theta_str}: sign-blind infidelity lower bound certified for all 256 syndromes", okf))

    r_d = d_blind_lower / d_undo_upper
    r_f = inf_blind_lower / inf_undo
    print(f" certified ratios: diamond >= {lo(r_d, 6)}; infidelity >= {lo(r_f, 6)}   ({time.time() - t0:.0f} s)")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--verify", action="store_true",
                        help="verify the saved certificates only (needs NumPy and mpmath, not CVXPY)")
    parser.add_argument("--certificates", default=CERT_PATH, help="certificate file (default: %(default)s)")
    args = parser.parse_args()

    print("=== Run configuration ===")
    print(f"script version = {SCRIPT_VERSION}")
    print(f"Python version = {platform.python_version()}")
    print(f"NumPy version = {np.__version__}")
    print(f"mpmath version = {mpmath.__version__}; interval precision = {iv.dps} decimal digits")
    if args.verify:
        print(f"mode = verify only; certificates from {args.certificates}")
    else:
        import cvxpy as cp
        print(f"mode = propose, save and verify; CVXPY version = {cp.__version__}; solver = {SOLVER} "
              "(used only to propose certificates)")
    print(f"thetas = {THETAS} (exact decimals); one round; L = 3")

    inst = build_instrument(3)
    N = branch_counts(inst)
    print(f"\nbranch counts: {int(N.sum())} error patterns, {int((N.sum(axis=(1, 2)) > 0).sum())} syndromes")

    if not args.verify:
        certs = {"script_version": SCRIPT_VERSION, "certificates": [propose(th, N, inst) for th in THETAS]}
        os.makedirs(os.path.dirname(args.certificates) or ".", exist_ok=True)
        with open(args.certificates, "w") as f:
            json.dump(certs, f)
        print(f"certificates saved to {args.certificates}")
    with open(args.certificates) as f:
        certs = json.load(f)

    passed = []
    for cert in certs["certificates"]:
        verify(cert, N, inst, passed)

    print("\n=== Summary ===")
    for label, ok in passed:
        print(f"{'PASS' if ok else 'FAIL'}  {label}")
    print("Certified bounds are proved by interval arithmetic from the saved certificates; "
          "the floating-point solver only proposes them.")


if __name__ == "__main__":
    main()
