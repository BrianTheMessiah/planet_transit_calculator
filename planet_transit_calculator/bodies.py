"""Registry of supported solar-system bodies and reference orbital data.

The semi-major axis and orbital period values below are standard
astronomical constants used only to size default search windows
(synodic period for departure dates, Hohmann time for flight duration).
They are not used in the delta-v calculations themselves, which rely on
real ephemeris positions/velocities (see ephemeris.py).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from attrs import define

AU_IN_KM = 149_597_870.7
SUN_MU_KM3_S2 = 1.32712440018e11

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

MERCURY_PERIHELION_DISTANCE_KM = 46.00e6
MERCURY_APHELION_DISTANCE_KM = 69.82e6

VENUS_PERIHELION_DISTANCE_KM = 107.48e6
VENUS_APHELION_DISTANCE_KM = 108.94e6

EARTH_PERIHELION_DISTANCE_KM = 147.09e6
EARTH_APHELION_DISTANCE_KM = 152.10e6

MARS_PERIHELION_DISTANCE_KM = 206.62e6
MARS_APHELION_DISTANCE_KM = 249.23e6

JUPITER_PERIHELION_DISTANCE_KM = 740.52e6
JUPITER_APHELION_DISTANCE_KM = 816.62e6

SATURN_PERIHELION_DISTANCE_KM = 1352.55e6
SATURN_APHELION_DISTANCE_KM = 1514.50e6

URANUS_PERIHELION_DISTANCE_KM = 2741.30e6
URANUS_APHELION_DISTANCE_KM = 3001.39e6

NEPTUNE_PERIHELION_DISTANCE_KM = 4444.45e6
NEPTUNE_APHELION_DISTANCE_KM = 4545.67e6


@dataclass
class CelestialBodyOrbit:
    """
    Orbital data for a celestial body, used for sizing search windows.

    Parameters:
        - `name`: Common name of the body (e.g. "Earth")
        - `ephemeris_key`: Key to use when looking up ephemeris data (e.g. "earth")
        - `perihelion_distance_km`: Closest distance to the Sun, in kilometers
        - `aphelion_distance_km`: Farthest distance to the Sun, in kilometers
        - `orbital_period_seconds`: Orbital period of the body around the Sun, in seconds
        - `mu`: Gravitational parameter of the central body (e.g. the Sun) in km^3/s^2. Not
          used for planets but can be set for moons.
        - `semi_major_axis_km`/`semi_major_axis_au`: Semi-major axis of the body's orbit
          around the Sun, in kilometers/astronomical units. Derived from
          `perihelion_distance_km` and `aphelion_distance_km`.
    """

    name: str
    ephemeris_key: str
    perihelion_distance_km: float  # closest distance to the Sun
    aphelion_distance_km: float  # farthest distance to the Sun
    orbital_period_seconds: float
    mu: float = 0.0  # gravitational parameter of the central body, km^3/s^2, not used for planets but can be set for moons

    orbit_radius: float = field(init=False)
    semi_major_axis_au: float = field(init=False)

    def __post_init__(self) -> None:
        semi_major_axis_km = self._compute_semi_major_axis_km()
        self.orbit_radius = semi_major_axis_km
        self.major_axis_au = self._convert_km_to_au(semi_major_axis_km)

    def _compute_semi_major_axis_km(self) -> float:
        """Semi-major axis is the average of the perihelion and aphelion distances."""
        return (self.perihelion_distance_km + self.aphelion_distance_km) / 2

    def _convert_km_to_au(self, km: float) -> float:
        """Convert kilometers to astronomical units."""
        return km / AU_IN_KM


CELESTIAL_BODIES: dict[str, CelestialBodyOrbit] = {
    MERCURY: CelestialBodyOrbit(
        MERCURY,
        MERCURY,
        MERCURY_PERIHELION_DISTANCE_KM,
        MERCURY_APHELION_DISTANCE_KM,
        MERCURY_ORBITAL_PERIOD_SECONDS,
    ),
    VENUS: CelestialBodyOrbit(
        VENUS,
        VENUS,
        VENUS_PERIHELION_DISTANCE_KM,
        VENUS_APHELION_DISTANCE_KM,
        VENUS_ORBITAL_PERIOD_SECONDS,
    ),
    EARTH: CelestialBodyOrbit(
        EARTH,
        EARTH,
        EARTH_PERIHELION_DISTANCE_KM,
        EARTH_APHELION_DISTANCE_KM,
        EARTH_ORBITAL_PERIOD_SECONDS,
    ),
    MARS: CelestialBodyOrbit(
        MARS,
        MARS,
        MARS_PERIHELION_DISTANCE_KM,
        MARS_APHELION_DISTANCE_KM,
        MARS_ORBITAL_PERIOD_SECONDS,
    ),
    JUPITER: CelestialBodyOrbit(
        JUPITER,
        JUPITER,
        JUPITER_PERIHELION_DISTANCE_KM,
        JUPITER_APHELION_DISTANCE_KM,
        JUPITER_ORBITAL_PERIOD_SECONDS,
    ),
    SATURN: CelestialBodyOrbit(
        SATURN,
        SATURN,
        SATURN_PERIHELION_DISTANCE_KM,
        SATURN_APHELION_DISTANCE_KM,
        SATURN_ORBITAL_PERIOD_SECONDS,
    ),
    URANUS: CelestialBodyOrbit(
        URANUS,
        URANUS,
        URANUS_PERIHELION_DISTANCE_KM,
        URANUS_APHELION_DISTANCE_KM,
        URANUS_ORBITAL_PERIOD_SECONDS,
    ),
    NEPTUNE: CelestialBodyOrbit(
        NEPTUNE,
        NEPTUNE,
        NEPTUNE_PERIHELION_DISTANCE_KM,
        NEPTUNE_APHELION_DISTANCE_KM,
        NEPTUNE_ORBITAL_PERIOD_SECONDS,
    ),
}

def get_celestial_body(name: str) -> CelestialBodyOrbit:
    key = name.strip().lower()
    if key not in CELESTIAL_BODIES:
        valid = ", ".join(sorted(CELESTIAL_BODIES))
        raise ValueError(f"Unknown body '{name}'. Valid options: {valid}")
    return CELESTIAL_BODIES[key]
