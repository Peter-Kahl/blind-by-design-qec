#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
audit_checks.py
===============

Version:
    0.1.0 (2026-10-04)

Supplementary research code for:

    Peter Kahl, 'Blind by Design: Quantum Error Correction, Functional
    Incompatibility and the Value of Evidence' (2026), §6.4.1 and Appendix A.

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
Construction audit of the Stim/PyMatching surface-code simulations of §6.4.1.
It reports validation, not findings: it checks that the pipeline behaves as
the construction rules of a detector error model (DEM) require.

The circuit is a rotated surface-code memory in the X basis, distance 5,
five rounds, with two-qubit depolarising noise at p = 0.005 after every
two-qubit gate (Stim's ``after_clifford_depolarization``). Labels used below:

- native   Stim's generated circuit with DEPOLARIZE1/DEPOLARIZE2 noise.
- B***     the *baseline* of the paper (Appendix A): the native circuit with
           each depolarising channel re-expressed as a PAULI_CHANNEL with
           individually adjustable components, one per qubit or qubit pair.
- A5       the paper's correlated error: a Y error on each of qubits 6, 7, 19
           and 20 occurring together with probability q, inserted into the
           baseline. It is the positive control and the mechanism used in
           §6.4.1.
- C_B      an earlier comparison family (YY x YY x YY on pairs (2,3), (6,7),
           (19,20)). It is retained only to show that its lone-D6 edge comes
           from Stim's hyperedge decomposition; the paper does not use it as
           evidence.
- e*       the matching edge from detector D6 directly to the boundary,
           flipping no logical observable.

Checks (pipeline mode)
----------------------
1. Native versus baseline: the cause of their small DEM discrepancy (Stim's
   approximate conversion of PAULI_CHANNEL noise, a second-order effect),
   tested by Stim's refusal of exact conversion, the sign of the gap, and its
   scaling when the error rate is halved.
2. The origin of a lone-D6 edge: a physical fault (A5) or Stim's
   decomposition step (C_B).
3. The geometry of D6: its position, its matching-graph edges, and its only
   route to the boundary.
4. (optional, --lemma-checks) Implementation checks of the lemma of
   Appendix A.3: e* is absent from the baseline matching graph after every
   single-channel deletion, and after random reweighting of every channel
   within the same support.

Modes
-----
    python audit_checks.py --pipeline                 # the paper's circuits
    python audit_checks.py --pipeline --lemma-checks  # plus check 4 (a few minutes)
    python audit_checks.py --mirror                   # native circuit, two stand-in baselines
    python audit_checks.py --native n.stim --bstar b.stim [--cb c.stim]
                           [--native-half nh.stim --bstar-half bh.stim]
    python audit_checks.py                            # self-contained demonstration

The self-contained demonstration uses a different Stim-generated circuit and
does not reproduce the paper's figures.

Requirements
------------
Python 3.8 or later; Stim; PyMatching (for matching-graph comparisons, used
in pipeline mode and check 4).
"""

from __future__ import annotations

import argparse
import platform
from collections import defaultdict

import numpy as np

import stim

SCRIPT_VERSION = "0.1.0"

DETECTOR = 6
PARTNER = 9
E_STAR = ((6, None), ())   # D6 to the boundary, no logical flip
NOISE_INSTRUCTIONS = {
    "DEPOLARIZE1", "DEPOLARIZE2", "X_ERROR", "Y_ERROR", "Z_ERROR",
    "PAULI_CHANNEL_1", "PAULI_CHANNEL_2",
}
COMMON_CLEAN_BOUNDARY = 55
PAIR_0 = (2, 3)
PAIR_78 = (6, 7)
PAIR_118 = (19, 20)
PAULIS_1 = ("X", "Y", "Z")
PAULIS_2 = tuple(
    (a, b) for a in "IXYZ" for b in "IXYZ" if (a, b) != ("I", "I")
)


def print_versions() -> None:
    """Print runtime versions needed to reproduce an audit."""
    try:
        import pymatching
        pm_version = pymatching.__version__
    except Exception:
        pm_version = "not installed"
    print("=== Run configuration ===")
    print(f"script version = {SCRIPT_VERSION}")
    print(f"Python version = {platform.python_version()}")
    print(f"Stim version = {stim.__version__}")
    print(f"PyMatching version = {pm_version}")


def xor_combine(p: float, q: float) -> float:
    """Return the probability that exactly one of two independent events occurs."""
    return p * (1 - q) + q * (1 - p)


def dem_table(dem: stim.DetectorErrorModel) -> dict[str, float]:
    """Map full DEM target strings to parity-merged probabilities."""
    table: defaultdict[str, float] = defaultdict(float)
    for instruction in dem.flattened():
        if instruction.type != "error":
            continue
        key = " ".join(str(t) for t in instruction.targets_copy())
        table[key] = xor_combine(table[key], instruction.args_copy()[0])
    return dict(table)


def signature(key: str) -> tuple[str, ...]:
    """Return a detector/observable signature with '^' separators XOR-reduced."""
    parity: defaultdict[str, int] = defaultdict(int)
    for token in key.split():
        if token != "^":
            parity[token] ^= 1
    return tuple(sorted(token for token, value in parity.items() if value))


def pieces(key: str) -> list[tuple[str, ...]]:
    """Return graph-like components of a possibly decomposed DEM error."""
    return [tuple(sorted(part.split())) for part in key.split(" ^ ")]


def build_dem(circuit: stim.Circuit, label: str, decompose: bool = True):
    """Build a DEM exactly where possible, otherwise use Stim's approximation."""
    try:
        dem = circuit.detector_error_model(decompose_errors=decompose)
        print(f"  [{label}] converted EXACTLY (no approximation needed)")
        return dem, False
    except ValueError as error:
        first_line = str(error).splitlines()[0][:120]
        print(f"  [{label}] exact conversion REFUSED by Stim: {first_line}")
        dem = circuit.detector_error_model(
            decompose_errors=decompose,
            approximate_disjoint_errors=True,
        )
        print(f"  [{label}] converted with approximate_disjoint_errors=True")
        return dem, True


def demo_native(p: float) -> stim.Circuit:
    """Build the self-contained distance-5 demonstration circuit."""
    return stim.Circuit.generated(
        "surface_code:rotated_memory_z",
        distance=5,
        rounds=5,
        after_clifford_depolarization=p,
        before_round_data_depolarization=p,
        before_measure_flip_probability=p,
        after_reset_flip_probability=p,
    )


def mirror_native(distance: int, rounds: int, p: float) -> stim.Circuit:
    """Reproduce the paper's native rotated-memory-X construction."""
    return stim.Circuit.generated(
        "surface_code:rotated_memory_x",
        distance=distance,
        rounds=rounds,
        after_clifford_depolarization=p,
    ).flattened()


def _pauli_target(letter: str, qubit: int):
    """Convert X/Y/Z plus a qubit index to the corresponding Stim target."""
    return {"X": stim.target_x, "Y": stim.target_y, "Z": stim.target_z}[letter](qubit)


def to_independent_components(circuit: stim.Circuit) -> stim.Circuit:
    """Rewrite depolarising channels as independent Pauli error components."""
    output = stim.Circuit()
    for instruction in circuit.flattened():
        if instruction.name == "DEPOLARIZE1":
            p = instruction.gate_args_copy()[0]
            for pauli in PAULIS_1:
                output.append(f"{pauli}_ERROR", instruction.targets_copy(), [p / 3])
        elif instruction.name == "DEPOLARIZE2":
            p = instruction.gate_args_copy()[0]
            targets = [target.value for target in instruction.targets_copy()]
            for q1, q2 in zip(targets[0::2], targets[1::2]):
                for a, b in PAULIS_2:
                    pauli_targets = []
                    if a != "I":
                        pauli_targets.append(_pauli_target(a, q1))
                    if b != "I":
                        pauli_targets.append(_pauli_target(b, q2))
                    output.append("CORRELATED_ERROR", pauli_targets, [p / 15])
        else:
            output.append(instruction)
    return output


def to_pauli_channel(circuit: stim.Circuit) -> stim.Circuit:
    """Rewrite depolarising instructions as PAULI_CHANNEL instructions."""
    output = stim.Circuit()
    for instruction in circuit:
        if isinstance(instruction, stim.CircuitRepeatBlock):
            output.append(
                stim.CircuitRepeatBlock(
                    instruction.repeat_count,
                    to_pauli_channel(instruction.body_copy()),
                )
            )
        elif instruction.name == "DEPOLARIZE1":
            p = instruction.gate_args_copy()[0]
            output.append("PAULI_CHANNEL_1", instruction.targets_copy(), [p / 3] * 3)
        elif instruction.name == "DEPOLARIZE2":
            p = instruction.gate_args_copy()[0]
            output.append("PAULI_CHANNEL_2", instruction.targets_copy(), [p / 15] * 15)
        else:
            output.append(instruction)
    return output


def signature_table(dem: stim.DetectorErrorModel) -> dict[str, float]:
    """Merge DEM errors by undecomposed XOR signature."""
    table: defaultdict[str, float] = defaultdict(float)
    for key, probability in dem_table(dem).items():
        sig = " ".join(signature(key))
        table[sig] = xor_combine(table[sig], probability)
    return dict(table)


def compare(dem_a, dem_b, label: str) -> float:
    """Compare probabilities of shared undecomposed DEM signatures."""
    a, b = signature_table(dem_a), signature_table(dem_b)
    shared = sorted(set(a) & set(b))
    differences = [b[key] - a[key] for key in shared]
    max_abs = max((abs(value) for value in differences), default=0.0)
    lower = sum(value < 0 for value in differences)
    higher = sum(value > 0 for value in differences)
    print(
        f"  {label}: {len(shared)} shared signatures | only in native: "
        f"{len(set(a) - set(b))} | only in B***: {len(set(b) - set(a))}"
    )
    print(f"  max |p_B*** - p_native| = {max_abs:.3e}")
    print(
        f"  B*** lower in {lower}, higher in {higher}, equal in "
        f"{len(differences) - lower - higher} of {len(differences)} shared signatures"
    )
    return max_abs


def check1(native, bstar, native_half=None, bstar_half=None) -> None:
    """Audit the source and scaling of native-versus-B*** DEM discrepancies."""
    print("\n=== CHECK 1: native vs B*** ===")
    native_dem, _ = build_dem(native, "native", decompose=False)
    bstar_dem, approximate = build_dem(bstar, "B***", decompose=False)
    full_gap = compare(native_dem, bstar_dem, "full error rate")
    if native_half is not None and bstar_half is not None:
        native_half_dem, _ = build_dem(native_half, "native, half rate", decompose=False)
        bstar_half_dem, _ = build_dem(bstar_half, "B***, half rate", decompose=False)
        half_gap = compare(native_half_dem, bstar_half_dem, "half error rate")
        if half_gap > 0:
            ratio = half_gap / full_gap if full_gap else float("nan")
            print(
                f"  scaling: max gap at half rate / max gap at full rate = {ratio:.3f} "
                "(about 0.25 supports a second-order approximation effect; about "
                "0.5 points to a first-order construction difference)"
            )
    else:
        print("  (scaling test skipped: supply --native-half and --bstar-half)")
    if approximate:
        print("  READING: Stim refused exact B*** conversion; approximation is implicated.")
    else:
        print("  READING: B*** converted exactly; inspect modelling/construction differences.")


def explain(circuit: stim.Circuit, keys: list[str], max_items: int = 6) -> None:
    """Ask Stim for representative circuit locations producing selected DEM errors."""
    if not keys:
        return
    filter_dem = stim.DetectorErrorModel()
    for key in keys[:max_items]:
        targets = []
        for token in key.split():
            if token == "^":
                continue
            if token.startswith("D"):
                targets.append(stim.target_relative_detector_id(int(token[1:])))
            elif token.startswith("L"):
                targets.append(stim.target_logical_observable_id(int(token[1:])))
        filter_dem.append("error", 0.001, targets)
    try:
        explained = circuit.explain_detector_error_model_errors(
            dem_filter=filter_dem,
            reduce_to_one_representative_error=True,
        )
        if not explained:
            print("    (Stim returned no circuit-location explanation.)")
        for item in explained:
            print("    ", str(item).replace("\n", "\n     ")[:1500])
    except Exception as error:
        print(f"    (Stim could not explain these errors: {error})")


def check2(circuit: stim.Circuit, label: str) -> None:
    """Determine whether the apparent lone-D6 edge is physical or decompositional."""
    print(f"\n=== CHECK 2: where does a lone-D{DETECTOR} edge come from? [{label}] ===")
    dem, _ = build_dem(circuit, label)
    table = dem_table(dem)
    touching = {
        key: p for key, p in table.items()
        if f"D{DETECTOR}" in key.replace("^", "").split()
    }
    direct, split = [], []
    accepted = [(f"D{DETECTOR}",), tuple(sorted((f"D{DETECTOR}", "L0")))]
    for key, probability in touching.items():
        for piece in pieces(key):
            if piece in accepted:
                (direct if "^" not in key else split).append((key, probability))
    print(f"  {len(touching)} DEM errors touch D{DETECTOR}. Listing all of them:")
    for key, probability in sorted(touching.items(), key=lambda item: -item[1]):
        flag = "   <-- decomposed ('^')" if "^" in key else ""
        print(f"    p={probability:.3e}  {key}{flag}")
    print(f"  Lone-D{DETECTOR} edge as a PHYSICAL fault (no '^'): {len(direct)}")
    print(f"  Lone-D{DETECTOR} edge created by DECOMPOSITION ('^'): {len(split)}")
    if direct or split:
        print("  Circuit fault locations behind these errors:")
        explain(circuit, [key for key, _ in direct + split])


def check3(circuit: stim.Circuit, label: str) -> None:
    """Locate D6 and audit its boundary/logical-error connectivity."""
    print(f"\n=== CHECK 3: where is D{DETECTOR}? [{label}] ===")
    coordinates = circuit.get_detector_coordinates()
    print(f"  circuit has {circuit.num_detectors} detectors, {circuit.num_observables} observable(s)")
    if not coordinates or DETECTOR not in coordinates:
        print("  no detector coordinates available")
        return
    xs = [c[0] for c in coordinates.values() if len(c) >= 2]
    ys = [c[1] for c in coordinates.values() if len(c) >= 2]
    for detector in (DETECTOR, PARTNER):
        if detector in coordinates:
            print(f"  D{detector}: coordinates {coordinates[detector]}")
    c = coordinates[DETECTOR]
    print(f"  detector x range {min(xs)}..{max(xs)}, y range {min(ys)}..{max(ys)}")
    print(f"  D{DETECTOR} on outermost x column: {c[0] in (min(xs), max(xs))}; on outermost y row: {c[1] in (min(ys), max(ys))}")
    if len(c) >= 3:
        times = sorted({value[2] for value in coordinates.values() if len(value) >= 3})
        print(f"  D{DETECTOR} is in time slice {c[2]} of {times}")
    dem, _ = build_dem(circuit, label)
    table = dem_table(dem)
    alone = [(k, p) for k, p in table.items() if signature(k) == (f"D{DETECTOR}",)]
    logical = [(k, p) for k, p in table.items() if signature(k) == tuple(sorted((f"D{DETECTOR}", "L0")))]
    print(f"  errors whose TOTAL signature is D{DETECTOR} alone: {len(alone)}")
    print(f"  errors whose TOTAL signature is D{DETECTOR} + L0: {len(logical)}")
    try:
        edges = matching_edges(circuit)
    except ImportError:
        print("  (PyMatching not installed: matching-graph edges of D6 not listed)")
        return
    print(f"  matching-graph edges at D{DETECTOR} (endpoints, logical fault ids, probability):")
    for key in sorted((k for k in edges if DETECTOR in k[0]), key=str):
        print(f"    {key}  p={edges[key]:.3e}")
    boundary = [k for k in edges if k[0] == (DETECTOR, None)]
    print(f"  boundary edges at D{DETECTOR}: {len(boundary)} (e* present: {E_STAR in edges})")
    for key in sorted((k for k in edges if k[0] == (PARTNER, None)), key=str):
        print(f"  boundary edge at partner D{PARTNER}: {key}  p={edges[key]:.3e}")


def _qubits(instruction) -> list[int]:
    """Extract qubit indices from a Stim instruction."""
    return [int(t.value) for t in instruction.targets_copy() if t.is_qubit_target]


def _pairs(instruction) -> list[tuple[int, int]]:
    """Extract consecutive two-qubit target pairs from a Stim instruction."""
    qubits = _qubits(instruction)
    return list(zip(qubits[0::2], qubits[1::2]))


def _pauli_pair_targets(pair, label):
    """Create Stim Pauli targets for a two-qubit Pauli label such as 'YY'."""
    return [_pauli_target(p, q) for q, p in zip(pair, label) if p != "I"]


def build_pipeline(distance: int, rounds: int, p: float):
    """Reconstruct native and B*** circuits from Experiment 7H-C3 v6."""
    clean = stim.Circuit.generated("surface_code:rotated_memory_x", distance=distance, rounds=rounds).flattened()
    native = stim.Circuit.generated("surface_code:rotated_memory_x", distance=distance, rounds=rounds, after_clifford_depolarization=p).flattened()
    clean_inst, native_inst = list(clean), list(native)
    mapping, i, j = {0: 0}, 0, 0
    while i < len(clean_inst):
        while j < len(native_inst) and native_inst[j].name in NOISE_INSTRUCTIONS:
            j += 1
        if str(clean_inst[i]) != str(native_inst[j]):
            raise RuntimeError("alignment failure; circuit differs from paper pipeline")
        mapping[i] = j
        i += 1
        j += 1
    common_noise_index = mapping[COMMON_CLEAN_BOUNDARY] - 1
    if native_inst[common_noise_index].name != "DEPOLARIZE2":
        raise RuntimeError("common layer is not DEPOLARIZE2")
    bstar = stim.Circuit()
    cb_insert = None
    expanded = {}
    for k, instruction in enumerate(native_inst):
        if instruction.name == "DEPOLARIZE1":
            pp = instruction.gate_args_copy()[0]
            for q in _qubits(instruction):
                bstar.append("PAULI_CHANNEL_1", [q], [pp / 3] * 3)
        elif instruction.name == "DEPOLARIZE2":
            pp = instruction.gate_args_copy()[0]
            for pair in _pairs(instruction):
                expanded.setdefault((k, pair), len(bstar))
                bstar.append("PAULI_CHANNEL_2", list(pair), [pp / 15] * 15)
            if k == common_noise_index:
                cb_insert = len(bstar)
        else:
            bstar.append(instruction)
    return clean, native, native_inst, bstar, cb_insert, expanded


def _one(name, targets, probability):
    """Return a one-instruction Stim circuit."""
    circuit = stim.Circuit()
    circuit.append(name, targets, probability)
    return circuit


def insert_at(circuit, index, targets, probability):
    """Insert one correlated-error instruction at a flattened circuit index."""
    return circuit[:index] + _one("CORRELATED_ERROR", targets, probability) + circuit[index:]


def pipeline_cb(bstar, index, q=1e-3, reps=("YY", "YY", "YY")):
    """Construct comparison family C_B from the paper pipeline."""
    targets = _pauli_pair_targets(PAIR_0, reps[0]) + _pauli_pair_targets(PAIR_78, reps[1]) + _pauli_pair_targets(PAIR_118, reps[2])
    return insert_at(bstar, index, targets, q)


def pipeline_a5(bstar, expanded, q=1e-3):
    """Construct the A5 correlated-pair surgical arm."""
    by_layer = defaultdict(dict)
    for (layer, pair), index in expanded.items():
        by_layer[layer][pair] = index
    layer = min(k for k, d in by_layer.items() if PAIR_78 in d and PAIR_118 in d)
    index = min(by_layer[layer][PAIR_78], by_layer[layer][PAIR_118])
    targets = _pauli_pair_targets(PAIR_78, "YY") + _pauli_pair_targets(PAIR_118, "YY")
    return insert_at(bstar, index, targets, q), index


def matching_edges(circuit):
    """Return PyMatching graph edges keyed by endpoints and logical fault IDs."""
    import pymatching
    dem = circuit.detector_error_model(decompose_errors=True, flatten_loops=True, approximate_disjoint_errors=True)
    matching = pymatching.Matching.from_detector_error_model(dem)
    output = {}
    for u, v, data in matching.edges():
        ends = (int(u), None) if v is None else tuple(sorted((int(u), int(v))))
        output[(ends, tuple(sorted(data.get("fault_ids", set()))))] = float(data["error_probability"])
    return output


def compare_edges(a, b, label: str) -> float:
    """Compare matching-edge probabilities and return the maximum shared-edge gap."""
    shared = set(a) & set(b)
    if not shared:
        print(f"  {label}: no shared matching edges")
        return 0.0
    differences = [b[key] - a[key] for key in shared]
    worst = max(shared, key=lambda key: abs(b[key] - a[key]))
    gap = abs(b[worst] - a[worst])
    print(f"  {label}: native edges {len(a)}, B*** edges {len(b)}, added {len(set(b)-set(a))}, removed {len(set(a)-set(b))}")
    print(f"  max |p_B*** - p_native| on matching edges = {gap:.3e} at edge {worst}")
    print(f"  B*** lower on {sum(d < 0 for d in differences)}, higher on {sum(d > 0 for d in differences)}, equal on {sum(d == 0 for d in differences)} edges")
    return gap


def run_pipeline(distance: int, rounds: int, p: float, lemma: bool = False) -> None:
    """Rebuild the paper pipeline and execute all three audit checks."""
    print(f"PIPELINE MODE: distance={distance}, rounds={rounds}, p={p}")
    _, native, _, bstar, cb_insert, expanded = build_pipeline(distance, rounds, p)
    _, native_half, _, bstar_half, _, _ = build_pipeline(distance, rounds, p / 2)
    print(f"  detectors {native.num_detectors}; B*** instructions {len(bstar)}; C_B insertion index {cb_insert}")
    check1(native, bstar, native_half, bstar_half)
    print("\n  Matching-edge comparison after decomposition:")
    full = compare_edges(matching_edges(native), matching_edges(bstar), "full error rate")
    half = compare_edges(matching_edges(native_half), matching_edges(bstar_half), "half error rate")
    print(f"  scaling on matching edges: half / full = {half / full:.3f}")
    cb = pipeline_cb(bstar, cb_insert)
    check2(cb, "C_B (YY x YY x YY, q=1e-3)")
    a5, index = pipeline_a5(bstar, expanded)
    print(f"\n  (A5 correlated event inserted at expanded index {index})")
    check2(a5, "A5 correlated pair (q=1e-3)")
    check3(native, "native circuit")
    if lemma:
        lemma_checks(bstar)


def _channel_occurrences(bstar):
    """List (instruction index, targets) for every PAULI_CHANNEL occurrence."""
    occurrences = []
    for k, instruction in enumerate(bstar):
        if instruction.name == "PAULI_CHANNEL_1":
            occurrences += [(k, (q,)) for q in _qubits(instruction)]
        elif instruction.name == "PAULI_CHANNEL_2":
            occurrences += [(k, pair) for pair in _pairs(instruction)]
    return occurrences


def _rebuild(bstar, drop=None, scale=None):
    """Rebuild the baseline, deleting one occurrence or rescaling every channel.

    ``drop`` is an (instruction index, targets) pair to omit; ``scale`` maps
    each occurrence to a vector of positive factors for its probabilities.
    """
    circuit = stim.Circuit()
    for k, instruction in enumerate(bstar):
        if instruction.name not in ("PAULI_CHANNEL_1", "PAULI_CHANNEL_2"):
            circuit.append(instruction)
            continue
        args = instruction.gate_args_copy()
        groups = [(q,) for q in _qubits(instruction)] if instruction.name == "PAULI_CHANNEL_1" else _pairs(instruction)
        for targets in groups:
            if drop == (k, targets):
                continue
            probs = list(args)
            if scale is not None:
                probs = [min(0.5, p * f) for p, f in zip(args, scale[(k, targets)])]
            circuit.append(instruction.name, list(targets), probs)
    return circuit


def lemma_checks(bstar, n_perturbations: int = 131, seed: int = 20260925) -> None:
    """Check 4: e* stays absent under deletions and reweightings of the baseline.

    The lemma of Appendix A.3 says that, under independent mechanisms, the
    matching graph's edge set is fixed by which mechanisms are active, their
    signatures and the decomposition rule, so reweighting cannot create e*,
    and deleting mechanisms can remove edges but (if decomposition depends only
    on the active signatures) cannot create e*. These are implementation
    checks of that statement, not a proof of it.
    """
    print("\n=== CHECK 4: implementation checks of the Appendix A.3 lemma ===")
    base_edges = set(matching_edges(bstar))
    occurrences = _channel_occurrences(bstar)
    print(f"  baseline: {len(base_edges)} matching edges; e* present: {E_STAR in base_edges}; "
          f"{len(occurrences)} channel occurrences")
    created, new_edges = 0, 0
    for occurrence in occurrences:
        edges = set(matching_edges(_rebuild(bstar, drop=occurrence)))
        created += E_STAR in edges
        new_edges += len(edges - base_edges)
    print(f"  single-channel deletions: {len(occurrences)} tested; e* created in {created}; "
          f"edges not in the baseline appeared in {new_edges} cases")
    rng = np.random.default_rng(seed)
    created, changed = 0, 0
    for _ in range(n_perturbations):
        scale = {occ: rng.uniform(0.1, 3.0, 15 if len(occ[1]) == 2 else 3) for occ in occurrences}
        edges = set(matching_edges(_rebuild(bstar, scale=scale)))
        created += E_STAR in edges
        changed += edges != base_edges
    print(f"  random reweightings within the same support: {n_perturbations} tested (seed {seed}); "
          f"e* created in {created}; edge set changed in {changed}")


def load(path: str | None):
    """Load a Stim circuit when a path was supplied."""
    return stim.Circuit.from_file(path) if path else None


def parse_args():
    """Parse command-line options."""
    parser = argparse.ArgumentParser(description="Construction audit for Blind by Design, §6.4.1 and Appendix A.")
    parser.add_argument("--native")
    parser.add_argument("--bstar")
    parser.add_argument("--cb")
    parser.add_argument("--native-half")
    parser.add_argument("--bstar-half")
    parser.add_argument("--p", type=float, default=0.005)
    parser.add_argument("--mirror", action="store_true")
    parser.add_argument("--distance", type=int, default=5)
    parser.add_argument("--rounds", type=int, default=5)
    parser.add_argument("--pipeline", action="store_true",
                        help="rebuild the paper's circuits (native, baseline, C_B, A5) and run checks 1-3")
    parser.add_argument("--lemma-checks", action="store_true",
                        help="with --pipeline, also run check 4 (deletions and reweightings; a few minutes)")
    return parser.parse_args()


def main() -> None:
    """Dispatch to demonstration, mirror, pipeline or saved-circuit mode."""
    args = parse_args()
    print_versions()
    if args.pipeline:
        run_pipeline(args.distance, args.rounds, args.p, lemma=args.lemma_checks)
        return
    if args.mirror:
        print(f"MIRROR MODE: distance={args.distance}, rounds={args.rounds}, p={args.p}")
        native = mirror_native(args.distance, args.rounds, args.p)
        native_half = mirror_native(args.distance, args.rounds, args.p / 2)
        print("\n--- B*** stand-in A: PAULI_CHANNEL rewrite ---")
        check1(native, to_pauli_channel(native), native_half, to_pauli_channel(native_half))
        print("\n--- B*** stand-in B: independent-component rewrite ---")
        check1(native, to_independent_components(native), native_half, to_independent_components(native_half))
        check2(native, "native circuit")
        check3(native, "native circuit")
        return
    if args.native:
        native = load(args.native)
        bstar = load(args.bstar)
        cb = load(args.cb)
        native_half = load(args.native_half)
        bstar_half = load(args.bstar_half)
        if bstar is not None:
            check1(native, bstar, native_half, bstar_half)
        check2(cb if cb is not None else native, "C_B" if cb is not None else "native")
        check3(native, "native")
        return
    print("DEMONSTRATION MODE (Stim-generated circuit; not the paper pipeline)")
    native = demo_native(args.p)
    native_half = demo_native(args.p / 2)
    check1(native, to_pauli_channel(native), native_half, to_pauli_channel(native_half))
    check2(native, "demo native")
    check3(native, "demo native")


if __name__ == "__main__":
    main()
