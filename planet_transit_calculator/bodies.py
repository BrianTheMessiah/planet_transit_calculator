"""Registry of supported solar-system bodies and reference orbital data.

The semi-major axis and orbital period values below are standard
astronomical constants used only to size default search windows
(synodic period for departure dates, Hohmann time for flight duration).
They are not used in the delta-v calculations themselves, which rely on
real ephemeris positions/velocities (see ephemeris.py).
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Body:
    name: str
    ephemeris_key: str
    semi_major_axis_au: float
    orbital_period_days: float


BODIES: dict[str, Body] = {
    "mercury": Body("mercury", "mercury", 0.387, 87.97),
    "venus": Body("venus", "venus", 0.723, 224.70),
    "earth": Body("earth", "earth", 1.000, 365.25),
    "mars": Body("mars", "mars", 1.524, 686.98),
    "jupiter": Body("jupiter", "jupiter", 5.203, 4332.59),
    "saturn": Body("saturn", "saturn", 9.537, 10759.22),
    "uranus": Body("uranus", "uranus", 19.191, 30688.5),
    "neptune": Body("neptune", "neptune", 30.069, 60182.0),
}


def get_body(name: str) -> Body:
    key = name.strip().lower()
    if key not in BODIES:
        valid = ", ".join(sorted(BODIES))
        raise ValueError(f"Unknown body '{name}'. Valid options: {valid}")
    return BODIES[key]


def synodic_period_days(a: Body, b: Body) -> float:
    """Time between successive similar alignments of the two bodies."""
    inv_diff = abs(1 / a.orbital_period_days - 1 / b.orbital_period_days)
    if inv_diff == 0:
        raise ValueError("Bodies with identical orbital periods have no synodic period")
    return 1 / inv_diff


def hohmann_tof_days(a: Body, b: Body) -> float:
    """One-way transfer time (days) of an idealized Hohmann transfer ellipse."""
    transfer_sma_au = (a.semi_major_axis_au + b.semi_major_axis_au) / 2
    return 365.25 * transfer_sma_au**1.5 / 2
