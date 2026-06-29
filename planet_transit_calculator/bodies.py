"""Registry of supported solar-system bodies and reference orbital data.

The semi-major axis and orbital period values below are standard
astronomical constants used only to size default search windows
(synodic period for departure dates, Hohmann time for flight duration).
They are not used in the delta-v calculations themselves, which rely on
real ephemeris positions/velocities (see ephemeris.py).

Orbital radii (semi-major axes) derived from perihelion/aphelion:

    Body      Perihelion (km)   Aphelion (km)    Orbit Radius (km)   Orbit Radius (AU)
    --------  ---------------   -------------    -----------------   -----------------
    Mercury        46.00e6          69.82e6            57.91e6              0.387
    Venus         107.48e6         108.94e6            108.21e6             0.723
    Earth         147.09e6         152.10e6            149.60e6             1.000
    Mars          206.62e6         249.23e6            227.93e6             1.524
    Jupiter       740.52e6         816.62e6            778.57e6             5.203
    Saturn       1352.55e6        1514.50e6           1433.53e6             9.579
    Uranus       2741.30e6        3001.39e6           2871.35e6            19.191
    Neptune      4444.45e6        4545.67e6           4495.06e6            30.069
"""

from __future__ import annotations

from dataclasses import dataclass, field


SUN_MU_KM3_S2 = 1.32712440018e11

MERCURY = "mercury"
VENUS = "venus"
EARTH = "earth"
MARS = "mars"
JUPITER = "jupiter"
SATURN = "saturn"
URANUS = "uranus"
NEPTUNE = "neptune"

class OrbitalPeriods:
    """Sidereal orbital periods in seconds."""
    MERCURY = 87.97 * 86400
    VENUS = 224.70 * 86400
    EARTH = 365.25 * 86400
    MARS = 686.98 * 86400
    JUPITER = 4332.59 * 86400
    SATURN = 10759.22 * 86400
    URANUS = 30688.5 * 86400
    NEPTUNE = 60182.0 * 86400


class PerihelionDistances:
    """Closest approach to the Sun, in kilometers."""
    MERCURY = 46.00e6
    VENUS = 107.48e6
    EARTH = 147.09e6
    MARS = 206.62e6
    JUPITER = 740.52e6
    SATURN = 1352.55e6
    URANUS = 2741.30e6
    NEPTUNE = 4444.45e6


class AphelionDistances:
    """Farthest point from the Sun, in kilometers."""
    MERCURY = 69.82e6
    VENUS = 108.94e6
    EARTH = 152.10e6
    MARS = 249.23e6
    JUPITER = 816.62e6
    SATURN = 1514.50e6
    URANUS = 3001.39e6
    NEPTUNE = 4545.67e6


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
        - `orbit_radius`/`semi_major_axis_au`: Semi-major axis of the body's orbit
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

    def __post_init__(self) -> None:
        semi_major_axis_km = self.compute_orbit_radius_km(self.perihelion_distance_km, self.aphelion_distance_km)
        self.orbit_radius = semi_major_axis_km

    @staticmethod
    def compute_orbit_radius_km(perihelion_km: float, aphelion_km: float) -> float:
        """Semi-major axis in km: the average of perihelion and aphelion distances."""
        return (perihelion_km + aphelion_km) / 2


CELESTIAL_BODIES: dict[str, CelestialBodyOrbit] = {
    MERCURY: CelestialBodyOrbit(
        MERCURY,
        MERCURY,
        PerihelionDistances.MERCURY,
        AphelionDistances.MERCURY,
        OrbitalPeriods.MERCURY,
    ),
    VENUS: CelestialBodyOrbit(
        VENUS,
        VENUS,
        PerihelionDistances.VENUS,
        AphelionDistances.VENUS,
        OrbitalPeriods.VENUS,
    ),
    EARTH: CelestialBodyOrbit(
        EARTH,
        EARTH,
        PerihelionDistances.EARTH,
        AphelionDistances.EARTH,
        OrbitalPeriods.EARTH,
    ),
    MARS: CelestialBodyOrbit(
        MARS,
        MARS,
        PerihelionDistances.MARS,
        AphelionDistances.MARS,
        OrbitalPeriods.MARS,
    ),
    JUPITER: CelestialBodyOrbit(
        JUPITER,
        JUPITER,
        PerihelionDistances.JUPITER,
        AphelionDistances.JUPITER,
        OrbitalPeriods.JUPITER,
    ),
    SATURN: CelestialBodyOrbit(
        SATURN,
        SATURN,
        PerihelionDistances.SATURN,
        AphelionDistances.SATURN,
        OrbitalPeriods.SATURN,
    ),
    URANUS: CelestialBodyOrbit(
        URANUS,
        URANUS,
        PerihelionDistances.URANUS,
        AphelionDistances.URANUS,
        OrbitalPeriods.URANUS,
    ),
    NEPTUNE: CelestialBodyOrbit(
        NEPTUNE,
        NEPTUNE,
        PerihelionDistances.NEPTUNE,
        AphelionDistances.NEPTUNE,
        OrbitalPeriods.NEPTUNE,
    ),
}

def get_celestial_body(name: str) -> CelestialBodyOrbit:
    key = name.strip().lower()
    if key not in CELESTIAL_BODIES:
        valid = ", ".join(sorted(CELESTIAL_BODIES))
        raise ValueError(f"Unknown body '{name}'. Valid options: {valid}")
    return CELESTIAL_BODIES[key]
