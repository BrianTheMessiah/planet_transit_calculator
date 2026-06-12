"""Porkchop contour plot generation."""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

from .search import SearchResult


def save_porkchop(result: SearchResult, path: str) -> None:
    """Render a departure-date vs time-of-flight contour of total delta-v."""
    grid = result.grid
    departure_dates = grid.departure_times.datetime
    tofs = grid.tofs_days

    x, y = np.meshgrid(departure_dates, tofs)
    z = grid.total_dv.T

    fig, ax = plt.subplots(figsize=(10, 6))
    finite = z[np.isfinite(z)]
    levels = np.linspace(finite.min(), np.percentile(finite, 95), 20)
    contour = ax.contourf(x, y, z, levels=levels, cmap="viridis", extend="max")
    fig.colorbar(contour, ax=ax, label="Total delta-v (km/s)")

    for opt in result.options:
        ax.plot(opt.departure_time.datetime, opt.tof_days, "r*", markersize=12)

    ax.set_xlabel("Departure date")
    ax.set_ylabel("Time of flight (days)")
    ax.set_title("Porkchop plot: total delta-v")
    fig.autofmt_xdate()
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
