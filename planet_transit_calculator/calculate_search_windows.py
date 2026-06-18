from attrs import define
import numpy as np
from planet_transit_calculator.bodies import CelestialBodyOrbitalData


@define
class CalcluateSearchWindows:
    origin: CelestialBodyOrbitalData
    destination: CelestialBodyOrbitalData
    mu: float

    min_time_of_flight_seconds: float | None = None
    max_time_of_flight_seconds: float | None = None
    n_depart: int = 80
    n_tof: int = 50

    def calculate_depature_and_times_of_flight(self, window_seconds):
        if window_seconds is None:
            window_seconds = self._compute_synodic_period()
        if (
            self.min_time_of_flight_seconds is None
            or self.max_time_of_flight_seconds is None
        ):
            hohmann = self._compute_hohmann_time_of_flight()
            if self.min_time_of_flight_seconds is None:
                self.min_time_of_flight_seconds = max(10.0, 0.4 * hohmann)
            if self.max_time_of_flight_seconds is None:
                self.max_time_of_flight_seconds = 1.8 * hohmann

        departure_offsets = np.linspace(0, window_seconds, self.n_depart)
        tofs = np.linspace(
            self.min_time_of_flight_seconds, self.max_time_of_flight_seconds, self.n_tof
        )

        return departure_offsets, tofs

    def _compute_synodic_period(
        self,
    ) -> float:
        """Time between successive similar alignments of the two bodies."""
        origin_mean_motion = 2 * np.pi / self.origin.orbital_period_seconds
        destination_mean_motion = 2 * np.pi / self.destination.orbital_period_seconds

        mean_motion_difference = abs(origin_mean_motion - destination_mean_motion)

        if mean_motion_difference == 0:
            raise ValueError("Identical orbital periods, synodic period is infinite")

        synodic_period_seconds = 2 * np.pi / mean_motion_difference
        return synodic_period_seconds

    def _compute_hohmann_time_of_flight(
        self,
    ) -> float:
        """One-way transfer time (days) of an idealized Hohmann transfer ellipse."""
        combine_orbital_axis_km = (
            self.origin.semi_major_axis_km + self.destination.semi_major_axis_km
        )

        transfer_semi_major_axis_km = combine_orbital_axis_km / 2

        time_of_flight_seconds = np.pi * np.sqrt(
            transfer_semi_major_axis_km**3 / self.mu
        )

        return time_of_flight_seconds
