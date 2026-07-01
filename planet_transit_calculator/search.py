"""Grid search over departure dates / flight times and Pareto-front selection."""

from __future__ import annotations

from dataclasses import dataclass

import astropy.units as units
import numpy as np
from astropy.time import Time

from .calculate_search_windows import CalcluateSearchWindows
from .bodies import (
    CelestialBodyOrbit,
)
from .ephemeris import MU_SUN, get_state
from .lambert import izzo_v0
from .transfer import SECONDS_PER_DAY, TransferResult


@dataclass
class PorkchopGrid:
    departure_times: Time  # shape (n_depart,)
    tofs_days: np.ndarray  # shape (n_tof,)
    total_dv: np.ndarray  # shape (n_depart, n_tof), km/s, NaN where no solution


@dataclass
class SearchResult:
    grid: PorkchopGrid
    options: list[TransferResult]


@dataclass
class _CandidateTransfer:
    """A single (departure, time of flight) grid point under consideration."""

    tof_days: float
    total_dv: float
    depart_index: int
    tof_index: int


def find_transfer_options(
    origin: CelestialBodyOrbit,
    destination: CelestialBodyOrbit,
    depart_after: Time,
    window_seconds: float | None = None,
    min_time_of_flight_seconds: float | None = None,
    max_time_of_flight_seconds: float | None = None,
    n_depart: int = 80,
    n_tof: int = 50,
    top_n: int = 5,
) -> SearchResult:
    """Search a grid of (departure date, time of flight) for transfer options.

    Returns the full porkchop grid plus a Pareto-optimal subset of `top_n`
    options trading off time of flight against total delta-v.
    """
    search_windows = CalcluateSearchWindows(
        origin_celestial_body=origin,
        destination_celestial_body=destination,
        mu=MU_SUN,
        min_time_of_flight_seconds=min_time_of_flight_seconds,
        max_time_of_flight_seconds=max_time_of_flight_seconds,
        n_depart=n_depart,
        n_tof=n_tof,
    )
    departure_offsets, times_of_flight = (
        search_windows.calculate_depature_and_times_of_flight(
            window_seconds=window_seconds
        )
    )

    (
        departure_times,
        origin_positions,
        origin_velocities,
        destination_positions,
        destination_velocities,
    ) = _lookup_body_states(
        origin,
        destination,
        depart_after,
        departure_offsets,
        times_of_flight,
        n_depart,
        n_tof,
    )

    total_dv_grid, departure_dv_grid, arrival_dv_grid = _compute_delta_v_grids(
        origin_positions,
        origin_velocities,
        destination_positions,
        destination_velocities,
        times_of_flight,
        n_depart,
        n_tof,
    )

    grid = PorkchopGrid(
        departure_times=departure_times,
        tofs_days=times_of_flight / SECONDS_PER_DAY,
        total_dv=total_dv_grid,
    )

    pareto_front = _find_pareto_front(times_of_flight, total_dv_grid, n_tof)
    selected_candidates = _select_top_n(pareto_front, top_n)
    options = _build_transfer_results(
        selected_candidates,
        depart_after,
        departure_offsets,
        departure_dv_grid,
        arrival_dv_grid,
    )

    return SearchResult(grid=grid, options=options)


def _lookup_body_states(
    origin: CelestialBodyOrbit,
    destination: CelestialBodyOrbit,
    depart_after: Time,
    departure_offsets: np.ndarray,
    times_of_flight: np.ndarray,
    n_depart: int,
    n_tof: int,
) -> tuple[Time, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Look up origin/destination position and velocity for every grid point.

    Returns the concrete departure times, the origin's state at each
    departure time (shape (n_depart, 3)), and the destination's state at
    every (departure, time of flight) pair (shape (n_depart, n_tof, 3)).
    """
    departure_times = depart_after + departure_offsets * units.second
    origin_positions, origin_velocities = get_state(
        origin.ephemeris_key, departure_times
    )

    arrival_offsets = departure_offsets[:, None] + times_of_flight[None, :]
    arrival_times = depart_after + arrival_offsets.ravel() * units.second
    destination_positions_flat, destination_velocities_flat = get_state(
        destination.ephemeris_key, arrival_times
    )
    destination_positions = destination_positions_flat.reshape(n_depart, n_tof, 3)
    destination_velocities = destination_velocities_flat.reshape(n_depart, n_tof, 3)

    return (
        departure_times,
        origin_positions,
        origin_velocities,
        destination_positions,
        destination_velocities,
    )


def _compute_delta_v_grids(
    origin_positions: np.ndarray,
    origin_velocities: np.ndarray,
    destination_positions: np.ndarray,
    destination_velocities: np.ndarray,
    times_of_flight: np.ndarray,
    n_depart: int,
    n_tof: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Solve a Lambert transfer for every grid point and record its delta-v cost.

    Returns (total_dv_grid, departure_dv_grid, arrival_dv_grid), each shape
    (n_depart, n_tof) with NaN where no transfer solution was found.
    """
    total_dv_grid = np.full((n_depart, n_tof), np.nan)
    departure_dv_grid = np.full((n_depart, n_tof), np.nan)
    arrival_dv_grid = np.full((n_depart, n_tof), np.nan)

    with np.errstate(invalid="ignore", divide="ignore"):
        for depart_index in range(n_depart):
            for tof_index in range(n_tof):
                try:
                    transfer_departure_velocity, transfer_arrival_velocity = izzo_v0(
                        MU_SUN,
                        origin_positions[depart_index],
                        destination_positions[depart_index, tof_index],
                        times_of_flight[tof_index],
                    )
                except (ValueError, RuntimeError):
                    continue

                departure_dv = np.linalg.norm(
                    transfer_departure_velocity - origin_velocities[depart_index]
                )
                arrival_dv = np.linalg.norm(
                    transfer_arrival_velocity
                    - destination_velocities[depart_index, tof_index]
                )
                departure_dv_grid[depart_index, tof_index] = departure_dv
                arrival_dv_grid[depart_index, tof_index] = arrival_dv
                total_dv_grid[depart_index, tof_index] = departure_dv + arrival_dv

    return total_dv_grid, departure_dv_grid, arrival_dv_grid


def _find_pareto_front(
    times_of_flight: np.ndarray, total_dv_grid: np.ndarray, n_tof: int
) -> list[_CandidateTransfer]:
    """Reduce the grid to its Pareto front of time-of-flight vs. total delta-v.

    First picks the best (lowest total dv) departure date for each time of
    flight, then keeps only the points where a longer time of flight buys a
    strictly lower delta-v.
    """
    best_per_tof = []
    for tof_index in range(n_tof):
        dv_column = total_dv_grid[:, tof_index]
        if np.all(np.isnan(dv_column)):
            continue
        best_depart_index = int(np.nanargmin(dv_column))
        best_per_tof.append(
            _CandidateTransfer(
                tof_days=times_of_flight[tof_index] / SECONDS_PER_DAY,
                total_dv=dv_column[best_depart_index],
                depart_index=best_depart_index,
                tof_index=tof_index,
            )
        )

    best_per_tof.sort(key=lambda candidate: candidate.tof_days)
    pareto_front = []
    best_dv_so_far = np.inf
    for candidate in best_per_tof:
        if candidate.total_dv < best_dv_so_far:
            pareto_front.append(candidate)
            best_dv_so_far = candidate.total_dv

    return pareto_front


def _select_top_n(
    pareto_front: list[_CandidateTransfer], top_n: int
) -> list[_CandidateTransfer]:
    """Evenly sample the Pareto front down to at most `top_n` candidates."""
    if len(pareto_front) <= top_n:
        return pareto_front

    sample_indices = sorted(
        set(np.linspace(0, len(pareto_front) - 1, top_n).round().astype(int))
    )
    return [pareto_front[index] for index in sample_indices]


def _build_transfer_results(
    candidates: list[_CandidateTransfer],
    depart_after: Time,
    departure_offsets: np.ndarray,
    departure_dv_grid: np.ndarray,
    arrival_dv_grid: np.ndarray,
) -> list[TransferResult]:
    """Convert selected grid points into fully described transfer results."""
    results = []
    for candidate in candidates:
        departure_time = (
            depart_after + departure_offsets[candidate.depart_index] * units.second
        )
        arrival_time = departure_time + candidate.tof_days * units.day
        results.append(
            TransferResult(
                departure_time=departure_time,
                arrival_time=arrival_time,
                tof_days=float(candidate.tof_days),
                departure_dv=float(
                    departure_dv_grid[candidate.depart_index, candidate.tof_index]
                ),
                arrival_dv=float(
                    arrival_dv_grid[candidate.depart_index, candidate.tof_index]
                ),
                total_dv=float(candidate.total_dv),
            )
        )
    return results
