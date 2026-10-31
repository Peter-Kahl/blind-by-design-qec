# Blind by Design: Computational Supplement

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

<p align="left">
  <img
    src="assets/blind-by-design-qec-3.jpg"
    alt="Blind by Design: schematic of a constrained measurement interface, quantum error-correction structure, and observable records."
  >
</p>

Computational supplement to:

**Peter Kahl, ‘Blind by Design: Quantum Error Correction, Functional Incompatibility and the Value of Evidence’ (2026). Lex et Ratio Working Paper LXR-2026-PHI-BBDQEC-WP, Version 1.0. doi:10.5281/zenodo.23278643.**

and to its technical companion:

**Peter Kahl, ‘The Cost of Sign-Blind Recovery under Coherent Noise: Worked Numerical Cases’ (2026). In preparation.**

The paper asks when evidence a system lacks marks a defect in how it inquires, and when it is required by what the system must do. Quantum error correction (QEC) is its test case: syndrome measurements are built to reveal errors while revealing nothing about the encoded information. The paper separates three kinds of blindness:

- **model-relative blindness**: evidence is present in the record but unused by the model through which the record is read;
- **observation-relative blindness**: a distinction the current procedures do not reveal, but some further procedure that preserves the system's function would; and
- **function-relative blindness**: a distinction no function-preserving procedure can reveal, as with an arbitrary unknown encoded quantum state.

Two propositions carry the function-relative case. Proposition 5.1 shows that an inquiry which performs its required operation exactly produces a classical record that is the same for every encoded state. Proposition 5.2 extends this to inquiries that perform the operation only approximately: if the inquiry falls short of its function by ε in diamond norm, after the best possible repair, its complete record distinguishes encoded states with total-variation distance at most 2√ε, and the shortfalls add over rounds.

This repository contains the five scripts behind the paper's computations, three further scripts extending one of them for the technical companion, the complete output of each reference run, and the reliance note on the preprint result the paper uses.

## Contents and correspondence with the paper

| Script | Paper | What it shows |
|---|---|---|
| `prop52_numerical_check.py` | §5.3.1; Appendix B | Numerical check of Proposition 5.2 (approximate classical-output blindness) on the toric and repetition codes; comparison of the complementary output with the record; pointwise versus uniform repair over all channels, with and without the syndrome history |
| `blind_by_design_toy_repetition_code.py` | §§6.2–6.3 | Observation-relative and function-relative blindness in one minimal system |
| `ozguler_prop1_check.py` | §6.1; reliance note | Independent numerical check of the toric-code sign-symmetry result used in §6.1 |
| `audit_checks.py` | §6.4.1; Appendix A | Construction audit of the surface-code simulations |
| `blind_by_design_continuous_misspecification.py` | §6.4.2 | Model-relative blindness, detected and repaired from the records already held |
| `two_round_repair_check.py` | Companion note, §3 (extends §5.3.1, Table 7) | Pointwise versus uniform repair over two rounds of the toric instrument, with the two-round syndrome history discarded or used |
| `branch_mechanism_check.py` | Companion note, §4 (and §3 for branch undo in diamond norm) | Why the sign matters: each syndrome branch is a coherent logical rotation whose direction reverses with the sign; branch-undo repair in infidelity and diamond norm |
| `interval_certification_check.py` | Companion note, §2 | Certified bounds, by interval arithmetic, behind the headline sign-blind to sign-aware ratios at θ = 0.05 and 0.1 |
| `infidelity_and_calibration_check.py` | Companion note, §4 | The same comparison by entanglement infidelity; infidelity accumulated over up to 300 rounds, against an idealised logical twirl; sensitivity of the known-sign repairs to calibration error, over one and many rounds |

## The scripts

### `prop52_numerical_check.py`

A numerical check of the inequality in Proposition 5.2. Two measures are used. The *diamond norm* measures the distance between two quantum channels as the largest trace distance between their outputs over all inputs, including inputs entangled with an untouched reference system; it is unnormalised here, with range 0 to 2. The *leak* is the largest total-variation distance between the syndrome-history laws of any two encoded states, with range 0 to 1; it bounds the advantage of any test, however computed, that tries to tell from the record which state was stored.

On the L = 3 toric-code instrument of Özgüler (2026), rebuilt by `ozguler_prop1_check.py` (two logical qubits, uniform coherent X rotations by θ, ideal syndrome extraction, minimum-weight recovery), the script computes for seven angles from θ = 0.05 to 0.5:

1. the exact one-round leak;
2. the diamond distance of the one-round logical channel from the identity, ε_id, by semidefinite programme (SDP);
3. the same distance after a unitary decoder fitted to the logical channel, ε_dec, which bounds from above the 'best possible repair' in the proposition; and
4. the bounds 2√ε_id and 2√ε_dec, and their ratio to the leak.

It then compares the record with the *complementary output* it is read from. The complementary channel keeps what the process passes to its environment: here, the syndrome register with the coherences between syndrome outcomes intact. The classical record is that output with the coherences discarded. For four angles the script finds, by local search over pairs of pure encoded states, how far apart the two complementary outputs can be, which is a lower bound on the true maximum. In the reference run the complementary output distinguishes encoded states about 14,000 times better than the record at θ = 0.05, growing as θ² rather than θ⁶, and comes within a factor of about 1.4 of the bound 2√ε_dec in the tested regime. Most of the observed gap between bound and record therefore lies in what the record discards. The complementary output is a mathematical construction from a dilation of the channel; how much of it is physically present in, or accessible from, a real apparatus depends on the implementation.

It then separates pointwise from uniform repair. The recoverability defect allows a different repair for each noise channel, but the record cannot reveal the rotation's sign, so a repair chosen from the record alone must serve both signs. By semidefinite programming over all repair channels, the script brackets four optima: the best repair for a known sign and for both signs, each with the syndrome history discarded and with the history used (one repair per class of syndromes whose branch maps agree up to weight). Each optimum is bracketed by the defect achieved by the solver's repairs, made exact channels (an upper bound), and by the dual objective at a dual point made exactly feasible (a lower bound). For the history-used repairs, a second lower bound is computed from the same dual point with the dual constraint made feasible for each of the 256 syndromes separately; it bounds the optimum over all history-conditioned policies without assuming that the grouping into classes is exact, and it agrees with the grouped bound to within 0.1 per cent for the sign-blind repair. With the history discarded, the fitted unitary repair is within about 3 per cent of the optimum, and no repair is within 3 per cent of the best sign-blind repair. With the history used, sign-blind repair improves by 1 to 29 per cent across the angles, but known-sign repair improves far more, so the cost of not knowing the sign grows: at θ = 0.05 the sign-blind optimum is at least 3.7 × 10⁻⁴ and the known-sign one at most 1.0 × 10⁻⁶. One round only.

It also repeats the instrument for N = 1, 2 and 3 rounds at θ = 0.1 and 0.3, compares the exact N-round leak with its bound, and checks that the N-round logical channel's distance from the identity is at most N times the one-round distance, as part (ii) of the proposition requires.

Two contrasts follow. In the three-qubit repetition code, every three-round history law is exactly independent of the encoded state, although the logical channel is not the identity: a blind record does not certify that the state is preserved, so the converse of the proposition fails. And a single sentinel qubit measured in the Y basis distinguishes the rotation's sign with per-shot total variation sin θ, whereas the passive syndrome record's sign discrimination is identically zero (checked by `ozguler_prop1_check.py`).

In the reference run the bound holds at every angle and round count. It is loose by a factor of about 20,000 at θ = 0.05 and about 7 at θ = 0.5. At small angles the leak grows as θ⁶ (that is, θ^2L), ε_id as θ³ (the coherent logical rotation) and ε_dec as θ⁴ (what remains once the fitted decoder removes that rotation).

Two safeguards keep the verdicts from depending on solver tolerance. Every reported diamond distance is an upper bound, because the SDP's dual solution is shifted until it satisfies the dual constraints, as judged by a floating-point eigenvalue computation, before the objective is evaluated. Every check of the inequality is made at a lower bound that does not come from the solver at all: the trace distance produced by one explicit input, the maximally entangled state. Both are numerically constructed estimates in double-precision floating point. The margins by which the checks pass far exceed the solver residuals, but the bounds have not been certified by interval arithmetic or validated eigenvalue bounds; the computations illustrate the proposition and test the implementation, and the proposition itself rests on its proof. One limit should be stated plainly: Proposition 5.2 is stated with the best repair over all channels, which the script cannot compute. It uses one fitted unitary repair, which can only overstate the defect, so its checks test a necessary consequence of the proposition. A violation would refute it; passing does not confirm the bound at its strongest. The SDP is validated against the closed form 2 sin(φ/2) for a single-qubit phase rotation. The leak is computed exactly from the parity form of the history effects, F_h = q_h I + g_h P with P the product of the two logical X operators; the script reports the largest deviation from that form, which is at rounding level throughout.

### `two_round_repair_check.py`

Section 4 of `prop52_numerical_check.py` computes the cost of not knowing the rotation's sign for one round only. The complete syndrome history has the same law at +θ and −θ (Özgüler 2026, Proposition 1), so a longer history cannot reveal the sign; but it might still help a repair that must serve both signs. This script repeats the comparison over two rounds of the same L = 3 toric instrument, at θ = 0.05 and 0.1, with one repair applied after both rounds. Repairs applied between rounds are not considered.

The 256² two-round histories are formed from ordered pairs of the 20 one-round syndrome classes and regrouped into classes whose branch maps agree up to weight: 201 classes at θ = 0.05 and 206 at θ = 0.1. Each class is represented by one Kraus operator per sign, scaled so that its map equals the class's summed map; the script reports the largest deviation, which is at rounding level. The four optima are bracketed as in `prop52_numerical_check.py`. For the sign-blind, history-used repair, classes with probability below 10⁻⁸ share one repair, since without this pooling the solver fails numerically (class probabilities span about 1 to 10⁻¹⁶). The pooled policy is a valid policy, so its achieved defect is a valid upper bound; the lower bound is computed for every class individually from the same dual point, so it does not depend on the pooling.

In the reference run at θ = 0.1, using the two-round history lowers the sign-blind defect by about 5 per cent (from 6.15 × 10⁻³ with no repair to at least 5.84 × 10⁻³), about the same as in one round. Every defect roughly doubles from one round to two. The solver's best known-sign repair reaches 1.8 × 10⁻⁵ at θ = 0.1 and 4.0 × 10⁻⁶ at θ = 0.05, so the script prints conservative ratios of 330 and 185. These are limited by the solver's precision, not by the problem. `branch_mechanism_check.py` constructs a repair of exactly the same form (one repair after both rounds, chosen from the complete history), branch undo, which reaches 9.95 × 10⁻⁶ and 1.69 × 10⁻⁷; it is a valid upper bound on the same optimum, and gives ratios of about 590 and 4,400, as in one round. The technical companion reports these tighter figures.

### `branch_mechanism_check.py`

Explains the infidelity results below, branch by branch, at θ = 0.05 and 0.1. For each syndrome s, K_s is the operator the round applies to the code space (rotation, syndrome projection and minimum-weight correction). The script reports, by syndrome class:

- how close K_s is to a multiple of a unitary;
- whether the branch at −θ is the complex conjugate of the branch at +θ;
- the logical rotation angle (the spread of the eigenphases of the nearest unitary); and
- the branch's share of one-round entanglement infidelity, from the exact decomposition 1 − F = Σ_s (p_s − |tr K_s|²/16).

In the reference run every branch is a coherent logical rotation (singular values equal to within 10⁻⁸ for the trivial syndrome and 10⁻⁴ for single-error syndromes at θ = 0.05). The branch at −θ is exactly the complex conjugate of the branch at +θ, and two classes of single-error syndromes, with total probability about 0.011 at θ = 0.05, rotate the logical state by θ and carry 95 per cent of the one-round infidelity (82 per cent at θ = 0.1).

*Branch undo*, the inverse of each branch's nearest unitary, applied with the sign known, removes 99.4 and 98.0 per cent of one-round infidelity at the two angles. At the wrong sign it is about four times worse than no repair. Branch undo is the correction Özgüler (2026) constructs by polar decomposition, his 'inverse-polar table'. As a cross-check the script reproduces his three 300-round infidelities at |θ| = 0.1 (0.1323553 with no correction, 0.0007456 with the correctly signed table, 0.4140468 at the wrong sign) to seven decimal places. The script also reports the diamond-norm defect of branch undo for one round at θ = 0.05 to 0.3 and for two rounds at θ = 0.05 and 0.1:

| θ | One round | Two rounds |
|---|---|---|
| 0.05 | 8.45 × 10⁻⁸ | 1.69 × 10⁻⁷ |
| 0.1 | 4.97 × 10⁻⁶ | 9.95 × 10⁻⁶ |

These are tighter than the solver's sign-known, history-used upper bounds in `prop52_numerical_check.py` (9.9 × 10⁻⁷ and 7.2 × 10⁻⁶) and `two_round_repair_check.py` (4.0 × 10⁻⁶ and 1.8 × 10⁻⁵). The conservative sign-blind to sign-known ratio is therefore about 4,400 at θ = 0.05 and 590 at θ = 0.1, the same over one and two rounds. For the rarest branches the eigenphase spread exceeds π, where it no longer measures the rotation angle reliably; the infidelity shares do not depend on it.

### `interval_certification_check.py`

The other scripts report primal–dual brackets constructed in double precision. This script certifies the four bounds behind the companion's headline ratios at θ = 0.05 and 0.1 (one round, L = 3) by interval arithmetic:

1. an upper bound on the diamond-norm defect of branch undo;
2. a lower bound on the diamond-norm defect of every sign-blind repair that may depend on the syndrome (one repair per syndrome, worst sign);
3. the entanglement infidelity of branch undo, enclosed exactly; and
4. a lower bound on the worst-sign infidelity of every sign-blind, syndrome-dependent repair.

Each branch operator is written exactly as K_σ = Σ_L a_{σ,L} X_L, where X_L ranges over the four logical X operators and a_{σ,L} = Σ_w N[σ, L, w] cos(θ/2)^(18−w) (−i sin(θ/2))^w, with integer counts N obtained by enumerating all 2¹⁸ X-error patterns. The script checks that the formula reproduces the state-vector instrument (to 3 × 10⁻¹⁵). In the logical-X basis every branch is diagonal. A floating-point solver proposes each certificate (a Watrous dual point for the upper bound; dual points of the repair problems for the lower bounds); mpmath interval arithmetic at 40 significant digits then proves feasibility by interval Cholesky factorisation for each of the 256 syndromes and encloses the bound. The docstring gives the derivations.

Reference run: branch undo has diamond defect at most 8.450 × 10⁻⁸ and 4.975 × 10⁻⁶, and infidelity 4.2210 × 10⁻⁸ and 2.4872 × 10⁻⁶; every sign-blind, syndrome-dependent repair has diamond defect at least 3.608 × 10⁻⁴ and 2.919 × 10⁻³, and worst-sign infidelity at least 7.178 × 10⁻⁶ and 1.2629 × 10⁻⁴, at θ = 0.05 and 0.1. The certified ratios are at least 4,269 and 586 in diamond norm and 170 and 50.7 in infidelity. The certified sign-blind diamond bounds are about 3 per cent weaker than the floating-point brackets of `prop52_numerical_check.py`, because of the shifts that make the certificates strictly feasible.

### `infidelity_and_calibration_check.py`

`prop52_numerical_check.py` and `two_round_repair_check.py` measure repairs in diamond norm, a worst-case measure that weighs coherent errors heavily. This script asks four further questions of the same L = 3 toric instrument.

1. **Infidelity.** The four optimal repairs are recomputed with entanglement infidelity (1 − entanglement fidelity) as the objective, at θ = 0.05, 0.1, 0.2 and 0.3, and bracketed. The upper end is the best of the solver's repairs made exact channels, branch undo, and no repair; the lower end is a dual bound made exactly feasible. In the reference run the sign-blind optimum, with the history used, is at least 171, 51, 19 and 13 times the known-sign one at the four angles. At θ = 0.1, branch undo meets the dual bound to within 0.1 per cent. With the history discarded, knowing the sign improves on no repair by at most 0.3 to 5.5 per cent.
2. **Accumulation.** Infidelity after N = 1 to 300 rounds, with a repair applied every round from that round's syndrome alone, so that the N-round channel is the N-th power of the repaired one-round channel. Policies: no repair; a unitary fitted at the right and at the wrong sign; the diamond-optimal history-used repair with the sign known; branch undo with the sign known; the diamond-optimal history-used repair for both signs; and, for reference, a logical Pauli twirl of the unrepaired channel, an idealisation of randomised compiling. Unrepaired infidelity grows faster than linearly (12 to 17 times from 10 to 100 rounds), the twirled channel's linearly. The best sign-blind policy computed is the twirl; the best known-sign policy, branch undo, beats it by a factor of 172 at θ = 0.05 and 51 at θ = 0.1, the same at every N. At θ = 0.1 the unrepaired infidelity after 300 rounds is 0.132; the diamond-optimal sign-blind repair is a poor controller by infidelity, because the diamond optimum is nearly degenerate.
3. **Calibration error.** Repairs fitted at θ(1 + η) are applied to the channel at θ, for η = 0, ±0.001, ±0.01 and 0.1. One-round diamond defects are reported for the fitted unitary, the diamond-optimal history-used repair and branch undo, and the infidelity after 1, 100 and 300 rounds for the last two. A 1 per cent error changes almost nothing. A 10 per cent error raises the one-round diamond defect of branch undo 60-fold at θ = 0.05, but after 300 rounds raises its infidelity only by a factor of 2.7 (θ = 0.05) and 1.5 (θ = 0.1), and it remains more than 60 and 30 times better than the twirl. Only fixed offsets are tested, not drift, correlated errors or non-uniform rotations.
4. **Sign error.** With the sign estimated once, wrong with probability q, and used for all 300 rounds, the expected infidelity of branch undo is (1 − q) I_right + q I_wrong. The number of independent sentinel shots (Y-basis readout of a qubit prepared in |0⟩, majority vote) needed to reach q is the smallest odd n whose binomial tail is at most q. In the reference run the advantage over the twirl falls to 2.9-fold (θ = 0.05) and 1.7-fold (θ = 0.1) at q = 0.05; break-even with the twirl is at q ≈ 0.15 and 0.09 (447 and 183 shots), and staying within a factor of 2 of the ideal needs q below 8.5 × 10⁻⁴ and 1.8 × 10⁻³ (about 3,900 and 850 shots).

### `blind_by_design_toy_repetition_code.py`

The three-qubit repetition code (|0_L⟩ = |000⟩, |1_L⟩ = |111⟩; stabilisers Z₁Z₂ and Z₂Z₃) under coherent rotations R_x(θ) = exp(−iθX/2), with majority-vote recovery. The script builds the syndrome instrument explicitly and checks that:

1. the syndrome probabilities are the same for every logical input (blindness to the stored state) and at +θ and −θ (blindness to the rotation's sign);
2. every syndrome effect is a multiple of the identity on the code space, and the effects sum to the identity;
3. the sign matters: a correction tuned to +θ restores the memory at +θ (entanglement fidelity 1.000 after 50 rounds) and damages it at −θ (0.913), against 0.976 with no extra correction;
4. a sentinel qubit measured in the Y basis reveals the sign (⟨Y⟩ = ∓0.199) without touching the memory; and
5. every four-round syndrome history has the same probability at +θ and −θ, for complex inputs and unequal angles on the three qubits reversed together.

Blindness to the sign is observation-relative: the sentinel is a function-preserving procedure that recovers it. Blindness to the stored state is function-relative: no procedure that reveals it leaves the stored state intact. The agreement in check 5 is a check on the implementation; the equality itself follows from a real-structure symmetry of the instrument (§6.2).

### `ozguler_prop1_check.py`

An independent reconstruction of the periodic toric-code instrument of Özgüler (2026, Proposition 1), built from the preprint's definitions without using its code: edge indexing, plaquette and star check matrices, logical frame and code states, reduced syndrome, minimum-weight X recovery with a lexicographic tie rule, and coherent X rotations on every edge. For L = 2 and L = 3 it checks sign symmetry of every syndrome effect, the real parity form of the effects, completeness, equality of multi-round history probabilities for random complex inputs, and, at L = 3, sign symmetry under unequal per-edge angles reversed together. Sanity checks confirm the code states and logical operators.

The proposition is stated for odd L ≥ 3. L = 2 is included as a contrast outside its scope: there the parity form fails, as expected, although sign symmetry holds. State-vector simulation scales as 2^(2L²), so L = 3 (18 qubits, 2¹⁸ amplitudes) is the largest tractable case. These are finite checks; the general result rests on the proof, re-derived in `docs/reliance_note_ozguler2026.md`. The instrument this script builds is also the one `prop52_numerical_check.py` uses.

### `audit_checks.py`

A construction audit of the Stim and PyMatching simulations of §6.4.1: a rotated surface-code memory, distance 5, five rounds, with two-qubit depolarising noise at p = 0.005. It reports validation, not findings. In pipeline mode it rebuilds the paper's circuits and runs:

1. **Native versus baseline.** The baseline re-expresses each depolarising channel as a PAULI_CHANNEL with individually adjustable components (labelled `B***` in the code). Its matching-graph probabilities differ from the native circuit's on all 502 edges, by at most 6.96 × 10⁻⁵, because Stim converts PAULI_CHANNEL noise only approximately. The script confirms this from Stim's refusal of exact conversion, the sign of the gap, and its scaling when the error rate is halved (0.259, against 0.25 expected for a second-order effect).
2. **The origin of a lone-D6 edge.** The paper's correlated error (a Y error on qubits 6, 7, 19 and 20 together, labelled `A5`) produces a genuine single fault whose signature is detector D6 alone. An earlier comparison family (labelled `C_B`) produces a lone-D6 edge only through Stim's hyperedge decomposition; the paper does not use it as evidence.
3. **The geometry of D6.** D6 is a first-round detector on the outermost row, with matching edges to D3, D9 and D13 only and no edge to the boundary; its only route to the boundary runs through D9.
4. **Implementation checks of the lemma in Appendix A.3** (optional, `--lemma-checks`). Deleting each of the 520 noise-channel occurrences in turn, and reweighting every channel at random within the same support (131 reweightings, fixed seed), never creates the boundary edge e* and never adds an edge.

Check 4 re-implements, in a single script, the deletion and reweighting audits run during development and reported in Appendix A.

### `blind_by_design_continuous_misspecification.py`

Continuous error correction on the three-qubit bit-flip code: each qubit flips at rate γ, and the stabilisers are measured continuously with record increments dQ_k = 2ηκ s_k dt + dW_k. A Wonham filter tracks the syndrome while assuming a detector efficiency of 1 when the true value is 0.5. The script shows that:

- the misspecification degrades syndrome tracking (3.1 against 6.7 per cent error at T = 25);
- an innovation statistic computed from the controller's own records flags it in every run; and
- re-estimating the efficiency by maximum likelihood from the same records, and rerunning the filter retrospectively, restores the tracking error to 3.1 per cent.

No new measurement is needed: the blindness is model-relative. With no feedback the problem is classical (a hidden Markov chain with a misspecified observation strength), and the methods are standard filtering techniques; the case is included to show the remedy at work. The null distribution of the innovation statistic is checked empirically, not derived.

## What the computations do not establish

- Finite numerical agreement does not prove a symmetry for arbitrary code size; the toric-code result rests on its proof.
- Proposition 5.2 rests on its proof, not on `prop52_numerical_check.py`. The script checks the inequality on two instruments, at seven angles and up to three rounds; it does not test adaptive inquiries, other codes or other noise.
- The fitted decoder in `prop52_numerical_check.py` is one unitary repair, not the best possible repair over all channels. Its checks test a consequence of Proposition 5.2, not the bound at its strongest.
- The complement distances in `prop52_numerical_check.py` come from a local search and are lower bounds on the true maxima.
- `infidelity_and_calibration_check.py` twirls the logical channel, not the physical operations as randomised compiling does, and does not combine twirling with syndrome-conditioned repair. Its accumulation results assume the same channel in every round, with no drift.
- `two_round_repair_check.py` considers one repair after two rounds. A controller that repairs between rounds, or that uses longer histories, has a larger class of policies, which has not been computed.
- How loose the bound is on these instruments says nothing about how loose it is in general.
- Small state-vector calculations say nothing about asymptotic behaviour.
- Numerical equality holds only to floating-point precision.
- A decoder's or controller's failure does not by itself show that information is absent from the record.
- Detecting one misspecification from existing records does not show that every blindness is recoverable that way.
- Failure to distinguish two conditions through a given syndrome interface does not show that they are physically indistinguishable under every measurement.
- The audit tests specified circuit constructions and their software representations, not quantum hardware. Its lemma checks support, but do not prove, the statement that Stim's decomposition depends only on which mechanisms are present.
- The surface-code decoding results of §6.4.1 (the logical-error penalty and the detectability of the correlated error) were produced by the author's development scripts; this repository contains the construction audit of those simulations, not the decoding runs themselves.

The analytical and interpretive claims depend on the definitions, assumptions and arguments of the paper.

## Requirements

- Python 3.9 or later for `prop52_numerical_check.py`; Python 3.8 or later for the other NumPy scripts; for `audit_checks.py`, a Python version supported by current Stim releases
- NumPy 1.20 or later for `prop52_numerical_check.py`; NumPy 1.17 or later for the other NumPy scripts
- CVXPY 1.4 or later, with the Clarabel solver installed by default (`prop52_numerical_check.py`, `two_round_repair_check.py`, `branch_mechanism_check.py` and `infidelity_and_calibration_check.py` only)
- Stim and PyMatching (`audit_checks.py` only)
- mpmath 1.3 or later (`interval_certification_check.py` only)

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

On macOS and some Linux systems the interpreter is called `python3`, and a bare `python` gives 'command not found'. Create the environment with `python3` as above; once it is activated, `python` refers to the environment's interpreter and the commands below work as written. The reference runs used Python 3.12; current NumPy releases need Python 3.11 or later.

No external data are needed: every script generates its own data.

## Running

```bash
python prop52_numerical_check.py
python blind_by_design_toy_repetition_code.py
python ozguler_prop1_check.py
python audit_checks.py --pipeline --lemma-checks
python blind_by_design_continuous_misspecification.py
python blind_by_design_continuous_misspecification.py --kappa 1
python two_round_repair_check.py
python branch_mechanism_check.py
python infidelity_and_calibration_check.py
python interval_certification_check.py
```

`prop52_numerical_check.py` imports `ozguler_prop1_check.py`, and the three companion scripts import both (`branch_mechanism_check.py` also imports `two_round_repair_check.py`), so run them from the repository folder. The last command reproduces the parameter-sensitivity figures quoted in §6.4.2. Run times on a laptop are seconds for the repetition-code and continuous scripts, under a minute for the toric-code check and for the audit without `--lemma-checks`, about one to two minutes for the audit with `--lemma-checks`, about 20 to 30 minutes for `prop52_numerical_check.py`, most of it in section 4, the optimal repairs, about 10 to 15 minutes for `two_round_repair_check.py` at its two default angles, about a minute for `branch_mechanism_check.py`, about 10 to 12 minutes for `infidelity_and_calibration_check.py`, and about 4 to 5 minutes for `interval_certification_check.py`. `audit_checks.py` also has a `--mirror` mode, a mode for saved `.stim` circuits, and a self-contained demonstration when run without arguments; the demonstration uses a different circuit and does not reproduce the paper's figures. Use `--help` on either of the two scripts with options.

## Reference output

The `output/` folder holds the complete, unedited console output of each reference run. The first six were produced by version **0.6.1** of the scripts, the seventh by version **0.7.0** of `two_round_repair_check.py`, the eighth by version **0.9.0** of `branch_mechanism_check.py`, the ninth by version **0.10.0** of `infidelity_and_calibration_check.py`, and the last by version **0.11.0** of `interval_certification_check.py` (with mpmath 1.3.0), all with Python 3.12.3, NumPy 2.4.4, Stim 1.16.0, PyMatching 2.4.0, CVXPY 1.9.3 and Clarabel 0.11.1:

| File | Command |
|---|---|
| `prop52_numerical_check_output.txt` | `python prop52_numerical_check.py` |
| `toy_repetition_code_output.txt` | `python blind_by_design_toy_repetition_code.py` |
| `ozguler_prop1_check_output.txt` | `python ozguler_prop1_check.py` |
| `audit_checks_pipeline_output.txt` | `python audit_checks.py --pipeline --lemma-checks` |
| `continuous_misspecification_output.txt` | `python blind_by_design_continuous_misspecification.py` |
| `continuous_misspecification_kappa1_output.txt` | `python blind_by_design_continuous_misspecification.py --kappa 1` |
| `two_round_repair_check_output.txt` | `python two_round_repair_check.py` |
| `branch_mechanism_check_output.txt` | `python branch_mechanism_check.py` |
| `infidelity_and_calibration_check_output.txt` | `python infidelity_and_calibration_check.py` |
| `interval_certification_check_output.txt` | `python interval_certification_check.py` |

The five outputs of the other scripts were regenerated in the same environment and are identical to the 0.1.0 reference files except for the version line. The elapsed-time line printed by the continuous script is omitted from its reference files, since it varies between runs. To compare a new run with a reference file:

```bash
python ozguler_prop1_check.py > my_run.txt
diff -u output/ozguler_prop1_check_output.txt my_run.txt
```

Exact textual identity is not guaranteed across Python, NumPy, Stim, PyMatching, CVXPY, solver or platform versions. Differences in the last digits of quantities at the level of floating-point rounding (around 10⁻¹⁶ or smaller) do not indicate a failure to reproduce. Diamond distances in `prop52_numerical_check.py` are computed by an iterative solver, and may differ in the third or fourth significant figure across solver versions; the solver reports some solves as 'optimal_inaccurate'. Its pass/fail verdicts do not depend on either, because of the bounds described above.

As a cross-platform check, `prop52_numerical_check.py` version 0.6.1 was also run on macOS with Python 3.14.7, NumPy 2.5.3, CVXPY 1.9.3 and Clarabel 0.11.1. All 14 checks passed, sections 1–3 and 5–6 matched the reference output to rounding, and every optimal-repair bracket in section 4 was consistent with the reference, with some lower bounds slightly lower: the fitted repair at θ = 0.05 came within 4 per cent of the optimum rather than 3, and within 0.5 per cent at θ = 0.1.

## Reproducibility

All randomness uses fixed seeds set in the source: 20261002 for the continuous simulation, 1 and 7 for the random inputs and unequal angles of the toric-code check, 20260925 for the audit's random reweightings, and 0 for the starting points of the complement search in `prop52_numerical_check.py`. The repetition-code script is deterministic. Each script prints its version and the package versions it ran with.

## Repository structure

```text
blind-by-design-qec/
├── README.md
├── CHANGELOG.md
├── LICENSE
├── CITATION.cff
├── .zenodo.json
├── .gitignore
├── requirements.txt
├── prop52_numerical_check.py
├── blind_by_design_toy_repetition_code.py
├── ozguler_prop1_check.py
├── audit_checks.py
├── blind_by_design_continuous_misspecification.py
├── two_round_repair_check.py
├── branch_mechanism_check.py
├── infidelity_and_calibration_check.py
├── interval_certification_check.py
├── docs/
│   └── reliance_note_ozguler2026.md
├── output/
│   ├── prop52_numerical_check_output.txt
│   ├── toy_repetition_code_output.txt
│   ├── ozguler_prop1_check_output.txt
│   ├── audit_checks_pipeline_output.txt
│   ├── continuous_misspecification_output.txt
│   ├── continuous_misspecification_kappa1_output.txt
│   ├── two_round_repair_check_output.txt
│   ├── branch_mechanism_check_output.txt
│   ├── infidelity_and_calibration_check_output.txt
│   └── interval_certification_check_output.txt
└── assets/
    └── blind-by-design-qec-3.jpg
```

## Version correspondence

The figures reported in Version 1.0 of the paper, including the numerical checks of §5.3.1, correspond to version **0.6.1** of this software. Version **0.7.0** adds `two_round_repair_check.py`, version **0.8.0** adds `infidelity_and_calibration_check.py`, and version **0.9.0** adds `branch_mechanism_check.py` and extends `infidelity_and_calibration_check.py`; their figures are reported in the technical companion. Version **0.9.1** changes documentation and metadata only. Version **0.10.0** adds the sign-error analysis to `infidelity_and_calibration_check.py`, and version **0.11.0** adds `interval_certification_check.py`. None of these versions changes the five earlier scripts or their reference outputs, so each reproduces every figure of the paper as well. Versions 0.1.0 to 0.6.0 accompanied earlier drafts; the results of the other four scripts are unchanged throughout. Use the archived release of the matching version when reproducing or citing figures. See `CHANGELOG.md` for the history of changes.

## Citation

If you use the theoretical argument, please cite the paper:

> Kahl, P. (2026) *Blind by Design: Quantum Error Correction, Functional Incompatibility and the Value of Evidence*. Lex et Ratio Working Paper LXR-2026-PHI-BBDQEC-WP, Version 1.0. doi:10.5281/zenodo.23278643.

For the two-round, mechanism, infidelity, accumulation and calibration results, please cite the technical companion:

> Kahl, P. (2026) *The Cost of Sign-Blind Recovery under Coherent Noise: Worked Numerical Cases*. In preparation.

If you use or modify the software, please also cite the archived software release. Citation metadata are in `CITATION.cff`. The paper and the software are separate scholarly objects with separate persistent identifiers.

The toric-code result checked by `ozguler_prop1_check.py`, and the instrument used by `prop52_numerical_check.py`, are due to:

> Özgüler, A.B. (2026) ‘Securing quantum error correction against misleading advice from AI agents’. arXiv:2609.19090 [quant-ph]. Preprint.

The proof of Proposition 5.2 rests on, and the diamond-norm computation in `prop52_numerical_check.py` uses, respectively:

> Kretschmann, D., Schlingemann, D. and Werner, R.F. (2008) ‘The information-disturbance tradeoff and the continuity of Stinespring’s representation’, *IEEE Transactions on Information Theory*, 54(4), pp. 1708–1717. doi:10.1109/TIT.2008.917696.

> Watrous, J. (2013) ‘Simpler semidefinite programs for completely bounded norms’, *Chicago Journal of Theoretical Computer Science*, 2013, Article 8. doi:10.4086/cjtcs.2013.008.

## Licence

The source code is released under the **MIT License** (see `LICENSE`). The accompanying paper and the reliance note are separate works, licensed under **CC BY 4.0**.

## Disclaimer

This software supports reproducibility and scrutiny of the computations accompanying *Blind by Design*. Its output should be read together with the assumptions, definitions, analytical results and limitations stated in the paper. It is not, by itself, sufficient to establish the paper's general claims about measurement, evidence or functional incompatibility.

The software is provided ‘as is’, without warranty of any kind, as specified in the MIT License.

## Author

**Peter Kahl**\
Independent researcher, Lex et Ratio\
ORCID: [0009-0003-1616-4843](https://orcid.org/0009-0003-1616-4843)\
[www.lexetratio.com](https://www.lexetratio.com)
