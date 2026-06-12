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
class CelestialBodyOrbitalInfo:
    name: str
    ephemeris_key: str
    semi_major_axis_au: float
    orbital_period_days: float


MERCURY_DAYS_PER_YEAR = 87.97
VENUS_DAYS_PER_YEAR = 224.70
EARTH_DAYS_PER_YEAR = 365.25
MARS_DAYS_PER_YEAR = 686.98
JUPITER_DAYS_PER_YEAR = 4332.59
SATURN_DAYS_PER_YEAR = 10759.22
URANUS_DAYS_PER_YEAR = 30688.5
NEPTUNE_DAYS_PER_YEAR = 60182.0

CELESTIAL_BODIES: dict[str, CelestialBodyOrbitalInfo] = {
    "mercury": CelestialBodyOrbitalInfo("mercury", "mercury", 0.387, MERCURY_DAYS_PER_YEAR),
    "venus": CelestialBodyOrbitalInfo("venus", "venus", 0.723, VENUS_DAYS_PER_YEAR),
    "earth": CelestialBodyOrbitalInfo("earth", "earth", 1.000, EARTH_DAYS_PER_YEAR),
    "mars": CelestialBodyOrbitalInfo("mars", "mars", 1.524, MARS_DAYS_PER_YEAR),
    "jupiter": CelestialBodyOrbitalInfo("jupiter", "jupiter", 5.203, JUPITER_DAYS_PER_YEAR),
    "saturn": CelestialBodyOrbitalInfo("saturn", "saturn", 9.537, SATURN_DAYS_PER_YEAR),
    "uranus": CelestialBodyOrbitalInfo("uranus", "uranus", 19.191, URANUS_DAYS_PER_YEAR),
    "neptune": CelestialBodyOrbitalInfo("neptune", "neptune", 30.069, NEPTUNE_DAYS_PER_YEAR),
}


def get_body(name: str) -> CelestialBodyOrbitalInfo:
    key = name.strip().lower()
    if key not in CELESTIAL_BODIES:
        valid = ", ".join(sorted(CELESTIAL_BODIES))
        raise ValueError(f"Unknown body '{name}'. Valid options: {valid}")
    return CELESTIAL_BODIES[key]


def compute_synodic_period_days(
    origin_celestial_body: CelestialBodyOrbitalInfo,
    destination_celestial_body: CelestialBodyOrbitalInfo,
) -> float:
    """Time between successive similar alignments of the two bodies."""
    inverted_origin_body_period = 1 / origin_celestial_body.orbital_period_days
    inverted_destination_body_period = (
        1 / destination_celestial_body.orbital_period_days
    )

    orbital_frequency_difference = abs(
        inverted_origin_body_period - inverted_destination_body_period
    )

    if orbital_frequency_difference == 0:
        raise ValueError("Bodies with identical orbital periods have no synodic period")

    return 1 / orbital_frequency_difference


def compute_hohmann_time_of_flight_in_earth_days(
    origin_celestial_body: CelestialBodyOrbitalInfo,
    destination_celestial_body: CelestialBodyOrbitalInfo,
) -> float:
    """One-way transfer time (days) of an idealized Hohmann transfer ellipse."""
    transfer_semi_major_axis_au = (
        origin_celestial_body.semi_major_axis_au
        + destination_celestial_body.semi_major_axis_au
    ) / 2

    orbital_period_years = transfer_semi_major_axis_au**1.5

    return EARTH_DAYS_PER_YEAR * orbital_period_years / 2
