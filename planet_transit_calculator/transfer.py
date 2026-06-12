"""Evaluation of a single Lambert transfer arc into a delta-v breakdown."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from astropy.time import Time

from .ephemeris import MU_SUN
from .lambert import izzo_v0

SECONDS_PER_DAY = 86400.0


@dataclass(frozen=True)
class TransferResult:
    departure_time: Time
    arrival_time: Time
    tof_days: float
    departure_dv: float  # km/s
    arrival_dv: float  # km/s
    total_dv: float  # km/s


def evaluate_transfer(
    r1: np.ndarray,
    v1_planet: np.ndarray,
    r2: np.ndarray,
    v2_planet: np.ndarray,
    departure_time: Time,
    arrival_time: Time,
    tof_days: float,
) -> TransferResult:
    """Solve Lambert's problem for r1 -> r2 and compute departure/arrival delta-v."""
    v1_transfer, v2_transfer = izzo_v0(MU_SUN, r1, r2, tof_days * SECONDS_PER_DAY)

    departure_dv = float(np.linalg.norm(v1_transfer - v1_planet))
    arrival_dv = float(np.linalg.norm(v2_transfer - v2_planet))

    return TransferResult(
        departure_time=departure_time,
        arrival_time=arrival_time,
        tof_days=tof_days,
        departure_dv=departure_dv,
        arrival_dv=arrival_dv,
        total_dv=departure_dv + arrival_dv,
    )
