# Reliance note: Özgüler (2026), Proposition 1

*Supplementary material for Kahl, P. (2026) 'Blind by Design: Quantum Error Correction, Functional Incompatibility and the Value of Evidence', Version 1.0 • Peter Kahl • 3 October 2026 (software references updated 4 October 2026)*

---

## 1. Purpose and scope

Section 6.1 of *Blind by Design* relies on a result from a recent, unrefereed preprint: Özgüler, A.B. (2026) 'Securing quantum error correction against misleading advice from AI agents', arXiv:2609.19090v1 [quant-ph]. Because the preprint is recent and the result structurally important to the paper, this note records exactly what is relied on, how it was checked, what the checks establish, and what they do not.

The note has three parts:

1. **Mathematical re-derivation** of the steps of the proposition's proof that the paper relies on.
2. **Numerical reconstruction** of the instrument at small distance.
3. **Credit allocation** between the preprint and earlier work, taken from the preprint's own comparison table.

It does not review the preprint as a whole.

## 2. The result relied on

**Özgüler (2026), Proposition 1**, equation (25) of the supplement (equation (2) of the main text).

In plain terms: in a periodic square toric code of odd distance L ≥ 3, subject to uniform coherent X rotations of angle θ on every data qubit, with ideal syndrome extraction and a fixed minimum-weight recovery, the probability distribution of the complete passive syndrome history is identical at +θ and −θ. This holds for every logical input, every angle and every history length.

A corollary follows from the identity of the distributions alone. No decoder, estimator or learned model, however large, can recover the sign of θ from the passive record.

**Assumptions under which the result holds**, which should accompany any citation:

- odd distance L ≥ 3;
- a periodic square toric code encoding two logical qubits;
- ideal preparation, syndrome extraction, readout and recovery;
- coherent X rotations as the only noise, uniform in magnitude as the theorem is stated;
- complete Z-syndrome extraction, with deterministic minimum-weight X recovery under the paper's fixed tie-break rule;
- 'passive record' meaning the syndrome history alone, with no terminal logical measurement.

The result is a property of this instrument. It is not a claim that inverse unitary operations are indistinguishable under every measurement, and the preprint itself says so.

## 3. Sources read

The full preprint was read, including the supplement. That covers:

- the main text (Sections I–VI);
- the supplement on the instrument (Sec. I), the exact-information results with Proposition 1 and its proof (Sec. II), logical dynamics and history bounds (Sec. III), recovery from the record (Sec. IV), finite-angle validation (Sec. V), and calibration and acceptance (Secs. VI–IX);
- channel proofs (Sec. A), disconnected supports (Sec. G), the character expansion (Sec. H), comparison with prior coefficient results (Sec. I), and the comparison with prior work and Table S5 (Sec. S).

## 4. Mathematical re-derivation

The printed proof of Proposition 1 uses two steps it does not spell out in the main text. Both were checked against the supplement, and both were also derived independently.

### 4.1 Evenness of error supports contained in a cut

The proof requires that every surviving error support f — a cycle contained in a cut — has even weight, including when f is disconnected. The supplement establishes this by decomposing cycles into circuits (Sec. G).

**Independent proof.** Let f be contained in the cut between a set S of plaquettes and its complement. Each edge of f has exactly one endpoint in S, so |f| equals the sum over plaquettes in S of the number of edges of f meeting each. Because f is a cycle, every plaquette meets an even number of its edges. A sum of even numbers is even, so |f| is even. Connectedness plays no role.

This rests on identifying the kernel of the check matrix H_Z with the cycle space and its row space with the cut space, as the preprint asserts in the proof of its Theorem 3.

### 4.2 The parity form and evenness in θ

The parts of Proposition 1 that the paper uses are that each syndrome effect has the form F_s = q_s I + b_s P, with q_s and b_s real (P being the product of the two logical X operators), and that these coefficients are even functions of θ. Both can be derived from the preprint's main-text expression for the corrected Kraus operators as combinations of logical X operators, without the character expansion of Sec. H:

1. **Fixed phases.** The coefficient of F_s = K_s†K_s on a logical class is a sum of terms from pairs of error supports with the same syndrome. Each term carries the phase i^(|a|+|a′|), up to sign. Since the two supports share a syndrome, their symmetric difference is a cycle of that class, whose weight has the parity of |a| + |a′|.
2. **Parity by class.** Trivial-class cycles are sums of star boundaries of weight 4, so they have even weight. Cycles in the two single-logical classes have odd weight, since their representatives have odd weight L. Cycles in the product class have even weight.
3. **Parity form.** Coefficients on the trivial and product classes are therefore real, and those on the single-logical classes purely imaginary. Since F_s and the logical operators are Hermitian and the operators linearly independent, every coefficient must be real, so the imaginary ones vanish.
4. **Evenness in θ.** Reversing θ multiplies each term by (−1)^(|a|+|a′|), which is +1 in the two surviving classes.

The history law then follows from the commutativity of the corrected Kraus operators, as in the printed proof, and sign symmetry follows from evenness.

### 4.3 A second route to sign symmetry

In the computational basis, the code states, stabiliser projectors and X-type recoveries are all real, and reversing θ is the same as taking the complex conjugate of the rotation. Hence K_s(−θ) equals the complex conjugate of K_s(θ), and once F_s is real, F_s(−θ) = F_s(θ). This is the author of this note's own argument, offered as a sufficient condition for the instruments considered, not as a general theorem.

### 4.4 Agreement with the supplement

- **Sec. G** proves the evenness claim by circuit decomposition and 'cycle–cut orthogonality'. That is the same fact as §4.1, and the two arguments agree.
- **Sec. H** derives the character expansion from four per-edge pair sums, with the recovery cancelling in F = K†K. This is consistent with §4.2. The support count reported there is confirmed (for example, 42 at L = 3).

## 5. Numerical reconstruction

### 5.1 Method

The instrument was rebuilt from the preprint's own definitions (supplement Sec. I A–B): the same edge labels, check matrices, logical frame and code states, the reduced syndrome, minimum-weight recovery with the stated tie rule, and X rotations on every edge. The full state was simulated at L = 3 (18 qubits, 256 syndromes) and, outside the proposition's scope, at L = 2.

### 5.2 Results

At θ = 0.1:

| Check | L = 3 | L = 2 (outside scope) |
|---|---|---|
| Effects identical at +θ and −θ (largest difference) | 4.2 × 10⁻²² | 0 |
| Effects real (largest imaginary part) | 2.1 × 10⁻²² | 0 |
| Effects of the form qI + bP (largest deviation) | 3.4 × 10⁻²¹ | 9.8 × 10⁻³ (form fails) |
| Completeness (effects sum to identity) | 1.0 × 10⁻¹⁴ | 5.6 × 10⁻¹⁶ |
| History probabilities at ±θ (largest difference, random complex inputs) | 1.1 × 10⁻¹⁶ (all two-round histories) | 0 (all three-round histories) |

At L = 3, with unequal angles on the 18 edges reversed together, the effects agree to 4.3 × 10⁻¹⁹, and K_s(−θ) equals the complex conjugate of K_s(θ) exactly.

The L = 3 results agree with the preprint's own numerical check, which reports a largest residual of 4.44 × 10⁻¹⁶ over all 256 syndromes. At L = 2 the parity form fails, as §4.2 predicts for even distance. Sign symmetry still holds, consistent with §4.3.

### 5.3 Software and environment

- **Script:** `ozguler_prop1_check.py`, version 0.1.0, in this repository. Besides the checks above, it confirms that the constructed code states are orthonormal and stabilised by every plaquette and star check, and that the logical operators commute with the plaquette checks.
- **Reference output:** `output/ozguler_prop1_check_output.txt`.
- **Environment of the reference run:** Python 3.12.3, NumPy 2.4.4. {Confirm on a second machine and record that environment here.}
- **Run time:** under a minute on a laptop.

## 6. What the checks do and do not establish

**They establish:**

- that the steps of Proposition 1's proof on which *Blind by Design* relies are mathematically sound, by an argument independent of the supplement's character expansion (§§4.1–4.2);
- that the printed definitions, implemented as stated, produce the claimed sign symmetry at L = 3 (§5).

**They do not establish:**

- **The result numerically for arbitrary distance or history length.** The numerical checks are finite: L = 3, two-round histories, θ = 0.1, plus the unequal-angle case. The general claim rests on the proof, re-derived in §4.
- **The preprint's other results.** These include Theorem 4 (the order and coefficient of input dependence), Theorem 5 (logical drift), Proposition 6 and Theorem 7 (history contraction and operational separation), Proposition 8, and the calibration, acceptance and drift-bound framework. *Blind by Design* cites the calibration framework only for what its abstract and Figs 2–3 report, not for any proof.
- **Anything beyond the stated assumptions** of §2: for example, noisy syndrome extraction, non-uniform noise beyond the reversed-together case checked, or other codes.

**An observation the paper does not use.** At L = 2 sign symmetry holds although the proposition does not cover even distance. By §4.3, sign symmetry needs only real syndrome effects, while odd distance is what produces the parity form. This rests on one even case and the note author's own argument. It is not attributed to the preprint and is not used in the paper.

## 7. Credit allocation

From the preprint's own comparison with prior work (Sec. S and Table S5):

- **Input blindness under coherent errors is prior work.** For odd-distance planar surface codes of the standard kind, syndrome probabilities independent of the encoded input are established in earlier work (Bravyi et al. 2018; Behrends and Béri 2025; Cheng et al. 2025).
- **Sign ambiguity in individual stabiliser outcomes is prior work** (Orsucci, Tiersch and Briegel 2016).
- **The preprint's contribution** is that sign symmetry holds for the complete joint distribution of arbitrarily long passive histories, which the single-outcome result does not by itself establish.
- **The basic criterion for input dependence** predates the preprint, which says so in its Sec. I, citing earlier coefficient results; the preprint's addition there is the exact order and coefficient, which *Blind by Design* does not use.

## 8. Citation form

*Blind by Design* cites the result as 'Özgüler (2026, Proposition 1); preprint'. Its citing sentence is:

> In a toric code subject to coherent rotations, the complete passive syndrome record is exactly blind to the sign of the rotation at every history length (Özgüler 2026, Proposition 1; preprint). The corresponding ambiguity in individual stabiliser outcomes was identified earlier (Orsucci, Tiersch and Briegel 2016), and blindness of syndrome statistics to the encoded input under coherent errors is established for odd-distance surface codes (Bravyi et al. 2018; Behrends and Béri 2025; Cheng et al. 2025).

If the preprint is revised or withdrawn, the paper's argument stands on the three-qubit example of its §6.2, which exhibits the same structure with exact, independently verified results. It would lose only the generality of a genuine two-dimensional topological code.

## References

Behrends, J. and Béri, B. (2025) ‘The surface code beyond Pauli channels: logical noise coherence, information-theoretic measures, and errorfield-double phenomenology’, PRX Quantum, 6, 040350.

Bravyi, S., Englbrecht, M., König, R. and Peard, N. (2018) ‘Correcting coherent errors with surface codes’, npj Quantum Information, 4, 55. doi: 10.1038/s41534-018-0106-y.

Cheng, Z., Huang, E., Khemani, V., Gullans, M.J. and Ippoliti, M. (2025) ‘Emergent unitary designs for encoded qubits from coherent errors and syndrome measurements’, PRX Quantum, 6, 030333.

Orsucci, D., Tiersch, M. and Briegel, H.J. (2016) ‘Estimation of coherent error sources from stabilizer measurements’, Physical Review A, 93, 042303.

Özgüler, A.B. (2026) ‘Securing quantum error correction against misleading advice from AI agents’. arXiv:2609.19090 [quant-ph]. Preprint.
