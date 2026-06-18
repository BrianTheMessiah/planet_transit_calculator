from attrs import define
import numpy as np
from planet_transit_calculator.bodies import CelestialBodyOrbit


@define
class CalcluateSearchWindows:
    origin_celestial_body: CelestialBodyOrbit
    destination_celestial_body: CelestialBodyOrbit
    mu: float

    min_time_of_flight_seconds: float | None = None
    max_time_of_flight_seconds: float | None = None
    n_depart: int = 80
    n_tof: int = 50

    def calculate_depature_and_times_of_flight(self, window_seconds):
        """
        Compute departure-time offsets and times of flight for a Lambert transfer.

        ### Parameters
        - window_seconds (float | None): Total search window in seconds. If None, uses the synodic period.

        ### Returns
        - departure_offsets (np.ndarray): Offsets from the initial epoch (s).
        - tofs (np.ndarray): Times of flight to evaluate (s).

        ### Notes
        Defines the departure/TOF grid for Lambert transfer search.
        """
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
        """
        Computes how often the two planets are in the same relative position, which is the natural period of transfer opportunities.

        ### Returns
        - synodic_period_seconds (float): Synodic period in seconds.

        ### Raises
        - ValueError: If the planets have identical orbital periods, resulting in an infinite synodic period.

        ### Notes
        The synodic period is calculated using the mean motions of the two planets.
        It represents the time it takes for the planets to realign in their orbits,
        which is crucial for planning interplanetary transfers. If the planets have identical orbital periods,
        the synodic period is infinite and a ValueError is raised.
        """
        origin_mean_motion = (
            2 * np.pi / self.origin_celestial_body.orbital_period_seconds
        )
        destination_mean_motion = (
            2 * np.pi / self.destination_celestial_body.orbital_period_seconds
        )

        mean_motion_difference = abs(origin_mean_motion - destination_mean_motion)

        if mean_motion_difference == 0:
            raise ValueError("Identical orbital periods, synodic period is infinite")

        synodic_period_seconds = 2 * np.pi / mean_motion_difference
        return synodic_period_seconds

    def _compute_hohmann_time_of_flight(
        self,
    ) -> float:
        """
        One-way transfer time (days) of an idealized Hohmann transfer ellipse.

        ### Returns
        - time_of_flight_seconds (float): Time of flight in seconds.

        ### Notes
        The Hohmann transfer is an efficient method for transferring between two circular orbits.
        """
        # Semi-major axis of the Hohmann transfer ellipse
        transfer_semi_major_axis = (
            self.origin_celestial_body.orbit_radius
            + self.destination_celestial_body.orbit_radius
        ) / 2

        return (
            np.sqrt(transfer_semi_major_axis**3 / self.mu) * np.pi
        )  # Time of flight in seconds
