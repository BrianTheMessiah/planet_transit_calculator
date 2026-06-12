"""Lambert's problem solver (Izzo 2015, zero-revolution case).

Given two position vectors and a time of flight, returns the velocity
vectors of the connecting Keplerian transfer orbit at each end.
"""

from __future__ import annotations

import numpy as np


def _hyp2f1b(x: float) -> float:
    """Hypergeometric function 2F1(3, 1, 5/2, x) via its series expansion."""
    if x >= 1.0:
        return np.inf
    res = term = 1.0
    ii = 0
    while True:
        term = term * (3 + ii) * (1 + ii) / (2.5 + ii) * x / (ii + 1)
        res_old = res
        res += term
        if res_old == res:
            return res
        ii += 1


def _compute_y(x: float, ll: float) -> float:
    return np.sqrt(1 - ll**2 * (1 - x**2))


def _compute_psi(x: float, y: float, ll: float) -> float:
    if -1 <= x < 1:
        return np.arccos(x * y + ll * (1 - x**2))
    if x > 1:
        return np.arcsinh((y - x * ll) * np.sqrt(x**2 - 1))
    return 0.0


def _tof_equation_y(x: float, y: float, t0: float, ll: float) -> float:
    if np.sqrt(0.6) < x < np.sqrt(1.4):
        eta = y - ll * x
        s1 = (1 - ll - x * eta) * 0.5
        q = 4 / 3 * _hyp2f1b(s1)
        t_ = (eta**3 * q + 4 * ll * eta) * 0.5
    else:
        psi = _compute_psi(x, y, ll)
        t_ = (psi / np.sqrt(np.abs(1 - x**2)) - x + ll * y) / (1 - x**2)
    return t_ - t0


def _tof_equation_p(x: float, y: float, t: float, ll: float) -> float:
    return (3 * t * x - 2 + 2 * ll**3 * x / y) / (1 - x**2)


def _tof_equation_p2(x: float, y: float, t: float, dt: float, ll: float) -> float:
    return (3 * t + 5 * x * dt + 2 * (1 - ll**2) * ll**3 / y**3) / (1 - x**2)


def _tof_equation_p3(x: float, y: float, dt: float, ddt: float, ll: float) -> float:
    return (7 * x * ddt + 8 * dt - 6 * (1 - ll**2) * ll**5 * x / y**5) / (1 - x**2)


def _initial_guess(t: float, ll: float) -> float:
    t0 = np.arccos(ll) + ll * np.sqrt(1 - ll**2)
    t1 = 2 * (1 - ll**3) / 3
    if t >= t0:
        return (t0 / t) ** (2 / 3) - 1
    if t < t1:
        return 5 / 2 * t1 / t * (t1 - t) / (1 - ll**5) + 1
    return np.exp(np.log(2) * np.log(t / t0) / np.log(t1 / t0)) - 1


def _householder(x0: float, t0: float, ll: float, atol: float, rtol: float, maxiter: int) -> float:
    x = x0
    for _ in range(maxiter):
        y = _compute_y(x, ll)
        fval = _tof_equation_y(x, y, t0, ll)
        t = fval + t0
        fder = _tof_equation_p(x, y, t, ll)
        fder2 = _tof_equation_p2(x, y, t, fder, ll)
        fder3 = _tof_equation_p3(x, y, fder, fder2, ll)
        x_new = x - fval * (
            (fder**2 - fval * fder2 / 2)
            / (fder * (fder**2 - fval * fder2) + fder3 * fval**2 / 6)
        )
        if abs(x_new - x) < rtol * abs(x) + atol:
            return x_new
        x = x_new
    raise RuntimeError("Lambert solver failed to converge")


def izzo_v0(
    mu: float,
    r1: np.ndarray,
    r2: np.ndarray,
    tof: float,
    prograde: bool = True,
    maxiter: int = 35,
    atol: float = 1e-8,
    rtol: float = 1e-10,
) -> tuple[np.ndarray, np.ndarray]:
    """Solve Lambert's problem for a single (zero-revolution) transfer arc.

    Parameters
    ----------
    mu : gravitational parameter of the central body (km^3/s^2)
    r1, r2 : departure/arrival position vectors (km)
    tof : time of flight (s)
    prograde : whether the transfer follows the system's prograde direction

    Returns
    -------
    v1, v2 : departure/arrival velocity vectors (km/s) on the transfer orbit
    """
    r1 = np.asarray(r1, dtype=float)
    r2 = np.asarray(r2, dtype=float)

    c = r2 - r1
    c_norm = np.linalg.norm(c)
    r1_norm = np.linalg.norm(r1)
    r2_norm = np.linalg.norm(r2)

    s = (r1_norm + r2_norm + c_norm) * 0.5

    i_r1 = r1 / r1_norm
    i_r2 = r2 / r2_norm
    i_h = np.cross(i_r1, i_r2)
    i_h = i_h / np.linalg.norm(i_h)

    ll = np.sqrt(max(0.0, 1 - min(1.0, c_norm / s)))

    if i_h[2] < 0:
        ll = -ll
        i_t1 = np.cross(i_r1, i_h)
        i_t2 = np.cross(i_r2, i_h)
    else:
        i_t1 = np.cross(i_h, i_r1)
        i_t2 = np.cross(i_h, i_r2)

    if not prograde:
        ll = -ll
        i_t1 = -i_t1
        i_t2 = -i_t2

    t = np.sqrt(2 * mu / s**3) * tof

    if abs(ll) >= 1:
        raise ValueError("Degenerate transfer geometry (collinear position vectors)")

    x0 = _initial_guess(t, ll)
    x = _householder(x0, t, ll, atol, rtol, maxiter)
    y = _compute_y(x, ll)

    gamma = np.sqrt(mu * s / 2)
    rho = (r1_norm - r2_norm) / c_norm
    sigma = np.sqrt(1 - rho**2)

    v_r1 = gamma * ((ll * y - x) - rho * (ll * y + x)) / r1_norm
    v_r2 = -gamma * ((ll * y - x) + rho * (ll * y + x)) / r2_norm
    v_t1 = gamma * sigma * (y + ll * x) / r1_norm
    v_t2 = gamma * sigma * (y + ll * x) / r2_norm

    v1 = v_r1 * i_r1 + v_t1 * i_t1
    v2 = v_r2 * i_r2 + v_t2 * i_t2
    return v1, v2
