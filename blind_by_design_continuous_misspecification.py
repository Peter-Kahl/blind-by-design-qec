#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
blind_by_design_continuous_misspecification.py
==============================================

Supplementary simulation for:

    Peter Kahl, "Blind by Design: Quantum Error Correction, Functional
    Incompatibility and the Value of Evidence" (2026), §6.4.

Author
------
Peter Kahl
Independent Researcher, Lex et Ratio
https://www.lexetratio.com
ORCID: 0009-0003-1616-4843

Purpose
-------
Demonstrate model-relative blindness in continuous quantum error correction and
show that it can be remedied using information already present in the
measurement record.

System
------
A three-qubit bit-flip code. Each physical qubit suffers bit flips at rate
``gamma``. The stabilisers S1 = Z1 Z2 and S2 = Z2 Z3 are measured continuously
with strength ``kappa`` and true detector efficiency ``eta_true``. With no
feedback, the syndrome sector is a classical continuous-time Markov chain and
measurement increments are modelled as

    dQ_k = 2 eta kappa s_k dt + dW_k.

See Ahn, Doherty and Landahl (2002) and van Handel and Mabuchi (2005).

Controller and tests
--------------------
A Wonham filter tracks the four syndrome sectors while assuming an efficiency
``eta_hat``. The script asks two questions.

1. Detection: can misspecification be detected from the record already held?
   The diagnostic is a normalised correlation between the innovation and the
   predicted syndrome. Under the correctly specified filter the innovation has
   zero conditional mean, so the statistic Z is approximately standard normal.
2. Remedy: can eta be re-estimated by maximum likelihood from the same records,
   and does using that estimate improve syndrome tracking?

Requirements
------------
Python 3.10+ and NumPy.

Run
---
    python blind_by_design_continuous_misspecification.py

The random seed and all default simulation parameters are fixed below for
reproducibility.

Licence
-------
MIT License. See LICENSE in the repository root.
SPDX-License-Identifier: MIT

Copyright (c) 2026 Peter Kahl.
"""

from __future__ import annotations

import platform
import time

import numpy as np


# Sector order: code (+,+), X1 (-,+), X2 (-,-), X3 (+,-).
SYNDROME_SIGNS = np.array(
    [[1, 1], [-1, 1], [-1, -1], [1, -1]], dtype=float
)


def simulate(
    eta_true: float,
    kappa: float,
    gamma: float,
    dt: float,
    n_steps: int,
    n_runs: int,
    seed: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Simulate true syndrome trajectories and continuous measurement records.

    Returns
    -------
    S
        Array of shape ``(n_steps, n_runs, 2)`` containing the two stabiliser
        signs at each time step.
    dQ
        Array of the same shape containing the corresponding noisy measurement
        increments.
    """
    rng = np.random.default_rng(seed)
    flipped = np.zeros((n_runs, 3), dtype=int)
    syndromes = np.empty((n_steps, n_runs, 2))
    records = np.empty((n_steps, n_runs, 2))

    for step in range(n_steps):
        physical_flips = rng.random((n_runs, 3)) < gamma * dt
        flipped ^= physical_flips

        s1 = 1 - 2 * ((flipped[:, 0] + flipped[:, 1]) % 2)
        s2 = 1 - 2 * ((flipped[:, 1] + flipped[:, 2]) % 2)
        syndromes[step, :, 0] = s1
        syndromes[step, :, 1] = s2

        noise = rng.normal(0, np.sqrt(dt), (n_runs, 2))
        records[step] = 2 * eta_true * kappa * syndromes[step] * dt + noise

    return syndromes, records


def run_filter(
    dQ: np.ndarray,
    eta_hats: np.ndarray,
    kappa: float,
    gamma: float,
    dt: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Run Wonham filters for all records and candidate detector efficiencies.

    Parameters
    ----------
    dQ
        Continuous measurement increments with shape ``(steps, runs, 2)``.
    eta_hats
        Candidate efficiencies used by the controller.
    kappa, gamma, dt
        Measurement strength, physical flip rate and integration time step.

    Returns
    -------
    log_likelihood
        Record log-likelihood for every run and candidate efficiency, with
        candidate-independent Gaussian constants omitted.
    z_statistic
        Normalised innovation/syndrome correlation for model criticism.
    map_path
        Maximum-a-posteriori syndrome-sector estimate at each time step.
    """
    n_steps, n_runs, _ = dQ.shape
    n_models = len(eta_hats)
    eta = np.asarray(eta_hats)[None, :, None]

    probabilities = np.zeros((n_runs, n_models, 4))
    probabilities[..., 0] = 1.0
    decay = np.exp(-4 * gamma * dt)

    log_likelihood = np.zeros((n_runs, n_models))
    numerator = np.zeros((n_runs, n_models))
    denominator = np.zeros((n_runs, n_models))
    map_path = np.empty((n_steps, n_runs, n_models), dtype=int)

    for step in range(n_steps):
        predicted = 0.25 + decay * (probabilities - 0.25)
        predicted_syndrome = predicted @ SYNDROME_SIGNS
        observation = dQ[step][:, None, :]

        innovation = observation - 2 * eta * kappa * predicted_syndrome * dt
        numerator += (innovation * predicted_syndrome).sum(axis=-1)
        denominator += (predicted_syndrome**2).sum(axis=-1) * dt

        raw_exponent = 2 * eta * kappa * (observation @ SYNDROME_SIGNS.T)
        maximum = raw_exponent.max(axis=-1, keepdims=True)
        weights = predicted * np.exp(raw_exponent - maximum)
        normaliser = weights.sum(axis=-1, keepdims=True)

        # Candidate-independent Gaussian terms are omitted. The final term is
        # the quadratic drift contribution for the two measurement channels.
        log_likelihood += (
            np.log(normaliser[..., 0])
            + maximum[..., 0]
            - 4 * (eta[..., 0] ** 2) * kappa**2 * dt
        )

        probabilities = weights / normaliser
        map_path[step] = probabilities.argmax(axis=-1)

    z_statistic = numerator / np.sqrt(denominator)
    return log_likelihood, z_statistic, map_path


def sector_index(syndromes: np.ndarray) -> np.ndarray:
    """Map pairs of stabiliser signs to the four-sector integer convention."""
    lookup = {(1, 1): 0, (-1, 1): 1, (-1, -1): 2, (1, -1): 3}
    indices = np.empty(syndromes.shape[:2], dtype=int)
    for signs, sector in lookup.items():
        mask = (
            (syndromes[..., 0] == signs[0])
            & (syndromes[..., 1] == signs[1])
        )
        indices[mask] = sector
    return indices


def run_experiment(
    duration: float,
    *,
    kappa: float,
    gamma: float,
    dt: float,
    eta_true: float,
    eta_assumed: float,
    n_runs: int,
    seed: int,
) -> None:
    """Run one duration of the misspecification experiment and print results."""
    n_steps = int(duration / dt)
    syndromes, records = simulate(
        eta_true, kappa, gamma, dt, n_steps, n_runs, seed
    )
    grid = np.round(np.arange(0.1, 1.51, 0.05), 3)
    log_likelihood, z_statistic, maps = run_filter(
        records, grid, kappa, gamma, dt
    )
    truth = sector_index(syndromes)

    true_index = int(np.where(np.isclose(grid, eta_true))[0][0])
    assumed_index = int(np.where(np.isclose(grid, eta_assumed))[0][0])
    eta_ml = grid[log_likelihood.argmax(axis=1)]

    tracking_error = (maps != truth[..., None]).mean(axis=0)
    ml_error = np.array(
        [tracking_error[r, log_likelihood[r].argmax()] for r in range(n_runs)]
    )

    z_true = z_statistic[:, true_index]
    z_wrong = z_statistic[:, assumed_index]

    print(
        f"\nT = {duration:g} ({n_steps} steps), {n_runs} runs, "
        f"true eta = {eta_true}, assumed eta = {eta_assumed}"
    )
    print(
        "  innovation statistic Z, correct model : "
        f"mean {z_true.mean():+.2f}, sd {z_true.std(ddof=1):.2f}, "
        f"|Z|>3 in {np.sum(np.abs(z_true) > 3)}/{n_runs} runs"
    )
    print(
        "  innovation statistic Z, assumed model : "
        f"mean {z_wrong.mean():+.2f}, sd {z_wrong.std(ddof=1):.2f}, "
        f"|Z|>3 in {np.sum(np.abs(z_wrong) > 3)}/{n_runs} runs"
    )
    print(
        "  maximum-likelihood eta from the same records: "
        f"mean {eta_ml.mean():.3f}, range {eta_ml.min():.2f}-{eta_ml.max():.2f}"
    )
    print("  syndrome-tracking error (fraction of time MAP sector wrong):")
    print(
        f"    true eta {eta_true}: {tracking_error[:, true_index].mean():.4f}   "
        f"assumed eta {eta_assumed}: {tracking_error[:, assumed_index].mean():.4f}   "
        f"re-estimated: {ml_error.mean():.4f}"
    )


def main() -> None:
    """Run the two published-duration simulations with fixed parameters."""
    start = time.time()
    print(f"Python {platform.python_version()} | numpy {np.__version__}")

    parameters = {
        "kappa": 4.0,
        "gamma": 0.05,
        "dt": 0.005,
        "eta_true": 0.5,
        "eta_assumed": 1.0,
        "n_runs": 20,
        "seed": 20261002,
    }
    for duration in (5.0, 25.0):
        run_experiment(duration, **parameters)

    print(f"\nelapsed {time.time() - start:.1f} s")


if __name__ == "__main__":
    main()
