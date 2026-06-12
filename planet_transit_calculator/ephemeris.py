"""Heliocentric ephemeris lookups backed by the JPL DE440s kernel via astropy.

Provides real planet position/velocity state vectors for any date (past,
present, or future within 1849-2150), which is the source of "real-time"
accuracy for the transfer calculations.
"""

from __future__ import annotations

import astropy.units as u
import numpy as np
from astropy.coordinates import get_body_barycentric_posvel, solar_system_ephemeris
from astropy.time import Time

MU_SUN = 1.32712440018e11  # km^3 / s^2

_configured = False


def configure() -> None:
    """Point astropy at the DE440s JPL ephemeris kernel (idempotent)."""
    global _configured
    if not _configured:
        solar_system_ephemeris.set("de440s")
        _configured = True


def get_state(body_key: str, time: Time) -> tuple[np.ndarray, np.ndarray]:
    """Heliocentric position (km) and velocity (km/s) of a body.

    `time` may be a scalar `Time` (returns shape (3,) arrays) or an array
    `Time` of length N (returns shape (N, 3) arrays).
    """
    configure()
    r_body, v_body = get_body_barycentric_posvel(body_key, time)
    r_sun, v_sun = get_body_barycentric_posvel("sun", time)

    r = (r_body - r_sun).get_xyz(xyz_axis=-1).to(u.km).value
    v = (v_body - v_sun).get_xyz(xyz_axis=-1).to(u.km / u.s).value
    return r, v
