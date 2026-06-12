"""Grid search over departure dates / flight times and Pareto-front selection."""

from __future__ import annotations

from dataclasses import dataclass

import astropy.units as u
import numpy as np
from astropy.time import Time

from .bodies import CelestialBodyOrbitalInfo, compute_hohmann_time_of_flight_in_earth_days, compute_synodic_period_days
from .ephemeris import MU_SUN, get_state
from .lambert import izzo_v0
from .transfer import TransferResult


@dataclass
class PorkchopGrid:
    departure_times: Time  # shape (n_depart,)
    tofs_days: np.ndarray  # shape (n_tof,)
    total_dv: np.ndarray  # shape (n_depart, n_tof), km/s, NaN where no solution


@dataclass
class SearchResult:
    grid: PorkchopGrid
    options: list[TransferResult]


def find_transfer_options(
    origin: CelestialBodyOrbitalInfo,
    destination: CelestialBodyOrbitalInfo,
    depart_after: Time,
    window_days: float | None = None,
    min_tof_days: float | None = None,
    max_tof_days: float | None = None,
    n_depart: int = 80,
    n_tof: int = 50,
    top_n: int = 5,
) -> SearchResult:
    """Search a grid of (departure date, time of flight) for transfer options.

    Returns the full porkchop grid plus a Pareto-optimal subset of `top_n`
    options trading off time of flight against total delta-v.
    """
    if window_days is None:
        window_days = compute_synodic_period_days(origin, destination)
    if min_tof_days is None or max_tof_days is None:
        hohmann = compute_hohmann_time_of_flight_in_earth_days(origin, destination)
        if min_tof_days is None:
            min_tof_days = max(10.0, 0.4 * hohmann)
        if max_tof_days is None:
            max_tof_days = 1.8 * hohmann

    departure_offsets = np.linspace(0, window_days, n_depart)
    tofs = np.linspace(min_tof_days, max_tof_days, n_tof)

    departure_times = depart_after + departure_offsets * u.day
    r1_all, v1_all = get_state(origin.ephemeris_key, departure_times)

    arrival_offsets = departure_offsets[:, None] + tofs[None, :]
    arrival_times_flat = depart_after + arrival_offsets.ravel() * u.day
    r2_flat, v2_flat = get_state(destination.ephemeris_key, arrival_times_flat)
    r2_all = r2_flat.reshape(n_depart, n_tof, 3)
    v2_all = v2_flat.reshape(n_depart, n_tof, 3)

    total_dv_grid = np.full((n_depart, n_tof), np.nan)
    departure_dv_grid = np.full((n_depart, n_tof), np.nan)
    arrival_dv_grid = np.full((n_depart, n_tof), np.nan)

    with np.errstate(invalid="ignore", divide="ignore"):
        for i in range(n_depart):
            for j in range(n_tof):
                try:
                    v1t, v2t = izzo_v0(MU_SUN, r1_all[i], r2_all[i, j], tofs[j] * 86400)
                except (ValueError, RuntimeError):
                    continue
                dep_dv = np.linalg.norm(v1t - v1_all[i])
                arr_dv = np.linalg.norm(v2t - v2_all[i, j])
                departure_dv_grid[i, j] = dep_dv
                arrival_dv_grid[i, j] = arr_dv
                total_dv_grid[i, j] = dep_dv + arr_dv

    grid = PorkchopGrid(
        departure_times=departure_times,
        tofs_days=tofs,
        total_dv=total_dv_grid,
    )

    # Best (lowest total dv) departure date for each flight time.
    best_per_tof = []
    for j in range(n_tof):
        col = total_dv_grid[:, j]
        if np.all(np.isnan(col)):
            continue
        i_best = int(np.nanargmin(col))
        best_per_tof.append((tofs[j], total_dv_grid[i_best, j], i_best, j))

    # Pareto front: as TOF increases, keep only points that strictly improve dv.
    best_per_tof.sort(key=lambda x: x[0])
    pareto = []
    min_dv_so_far = np.inf
    for tof, dv, i, j in best_per_tof:
        if dv < min_dv_so_far:
            pareto.append((tof, dv, i, j))
            min_dv_so_far = dv

    if len(pareto) <= top_n:
        selected = pareto
    else:
        indices = sorted(set(np.linspace(0, len(pareto) - 1, top_n).round().astype(int)))
        selected = [pareto[k] for k in indices]

    options = []
    for tof, dv, i, j in selected:
        dep_time = depart_after + departure_offsets[i] * u.day
        arr_time = dep_time + tofs[j] * u.day
        options.append(
            TransferResult(
                departure_time=dep_time,
                arrival_time=arr_time,
                tof_days=float(tofs[j]),
                departure_dv=float(departure_dv_grid[i, j]),
                arrival_dv=float(arrival_dv_grid[i, j]),
                total_dv=float(dv),
            )
        )

    return SearchResult(grid=grid, options=options)
