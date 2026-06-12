"""Registry of supported solar-system bodies and reference orbital data.

The semi-major axis and orbital period values below are standard
astronomical constants used only to size default search windows
(synodic period for departure dates, Hohmann time for flight duration).
They are not used in the delta-v calculations themselves, which rely on
real ephemeris positions/velocities (see ephemeris.py).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


AU_IN_KM = 149597870.7

MERCURY = "mercury"
VENUS = "venus"
EARTH = "earth"
MARS = "mars"
JUPITER = "jupiter"
SATURN = "saturn"
URANUS = "uranus"
NEPTUNE = "neptune"

MERCURY_ORBITAL_PERIOD_SECONDS = 87.97 * 86400
VENUS_ORBITAL_PERIOD_SECONDS = 224.70 * 86400
EARTH_ORBITAL_PERIOD_SECONDS = 365.25 * 86400
MARS_ORBITAL_PERIOD_SECONDS = 686.98 * 86400
JUPITER_ORBITAL_PERIOD_SECONDS = 4332.59 * 86400
SATURN_ORBITAL_PERIOD_SECONDS = 10759.22 * 86400
URANUS_ORBITAL_PERIOD_SECONDS = 30688.5 * 86400
NEPTUNE_ORBITAL_PERIOD_SECONDS = 60182.0 * 86400

MERCURY_SEMI_MAJOR_AXIS_AU = 0.387
VENUS_SEMI_MAJOR_AXIS_AU = 0.723
EARTH_SEMI_MAJOR_AXIS_AU = 1.000
MARS_SEMI_MAJOR_AXIS_AU = 1.524
JUPITER_SEMI_MAJOR_AXIS_AU = 5.203
SATURN_SEMI_MAJOR_AXIS_AU = 9.537
URANUS_SEMI_MAJOR_AXIS_AU = 19.191
NEPTUNE_SEMI_MAJOR_AXIS_AU = 30.069


@dataclass(frozen=True)
class CelestialBodyOrbitalInfo:
    name: str
    ephemeris_key: str
    semi_major_axis_au: float
    orbital_period_seconds: float
    mu: float = 0.0  # gravitational parameter of the central body, km^3/s^2, not used for planets but can be set for moons


CELESTIAL_BODIES: dict[str, CelestialBodyOrbitalInfo] = {
    MERCURY: CelestialBodyOrbitalInfo(
        MERCURY, MERCURY, MERCURY_SEMI_MAJOR_AXIS_AU, MERCURY_ORBITAL_PERIOD_SECONDS
    ),
    VENUS: CelestialBodyOrbitalInfo(
        VENUS, VENUS, VENUS_SEMI_MAJOR_AXIS_AU, VENUS_ORBITAL_PERIOD_SECONDS
    ),
    EARTH: CelestialBodyOrbitalInfo(
        EARTH, EARTH, EARTH_SEMI_MAJOR_AXIS_AU, EARTH_ORBITAL_PERIOD_SECONDS
    ),
    MARS: CelestialBodyOrbitalInfo(
        MARS, MARS, MARS_SEMI_MAJOR_AXIS_AU, MARS_ORBITAL_PERIOD_SECONDS
    ),
    JUPITER: CelestialBodyOrbitalInfo(
        JUPITER, JUPITER, JUPITER_SEMI_MAJOR_AXIS_AU, JUPITER_ORBITAL_PERIOD_SECONDS
    ),
    SATURN: CelestialBodyOrbitalInfo(
        SATURN, SATURN, SATURN_SEMI_MAJOR_AXIS_AU, SATURN_ORBITAL_PERIOD_SECONDS
    ),
    URANUS: CelestialBodyOrbitalInfo(
        URANUS, URANUS, URANUS_SEMI_MAJOR_AXIS_AU, URANUS_ORBITAL_PERIOD_SECONDS
    ),
    NEPTUNE: CelestialBodyOrbitalInfo(
        NEPTUNE, NEPTUNE, NEPTUNE_SEMI_MAJOR_AXIS_AU, NEPTUNE_ORBITAL_PERIOD_SECONDS
    ),
}


def get_body(name: str) -> CelestialBodyOrbitalInfo:
    key = name.strip().lower()
    if key not in CELESTIAL_BODIES:
        valid = ", ".join(sorted(CELESTIAL_BODIES))
        raise ValueError(f"Unknown body '{name}'. Valid options: {valid}")
    return CELESTIAL_BODIES[key]


def compute_synodic_period(
    origin: CelestialBodyOrbitalInfo,
    destination: CelestialBodyOrbitalInfo,
) -> float:
    """Time between successive similar alignments of the two bodies."""
    origin_mean_motion = 2 * np.pi / origin.orbital_period_seconds
    destination_mean_motion = 2 * np.pi / destination.orbital_period_seconds

    mean_motion_difference = abs(origin_mean_motion - destination_mean_motion)

    if mean_motion_difference == 0:
        raise ValueError("Identical orbital periods, synodic period is infinite")

    synodic_period_seconds = 2 * np.pi / mean_motion_difference
    return synodic_period_seconds


def compute_hohmann_time_of_flight(
    origin: CelestialBodyOrbitalInfo,
    destination: CelestialBodyOrbitalInfo,
    mu: float = 1.32712440018e11,  # gravitational parameter of the Sun in km^3/s^2
) -> float:
    """One-way transfer time (days) of an idealized Hohmann transfer ellipse."""
    origin_axis_km = origin.semi_major_axis_au * AU_IN_KM
    destination_axis_km = destination.semi_major_axis_au * AU_IN_KM

    combine_orbital_axis_km = origin_axis_km + destination_axis_km

    transfer_semi_major_axis_km = combine_orbital_axis_km / 2

    time_of_flight_seconds = np.pi * np.sqrt(transfer_semi_major_axis_km**3 / mu)

    return time_of_flight_seconds
