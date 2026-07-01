"""Regression test: search grid offsets are in seconds, not days."""

import unittest

from astropy.time import Time

from planet_transit_calculator.bodies import get_celestial_body
from planet_transit_calculator.search import find_transfer_options

SECONDS_PER_DAY = 86400.0


class TestSearchWindowUnits(unittest.TestCase):
    def test_departure_and_arrival_times_stay_within_requested_window(self):
        origin = get_celestial_body("earth")
        destination = get_celestial_body("mars")
        depart_after = Time("2026-07-01")

        window_days = 30
        min_tof_days = 100
        max_tof_days = 300

        result = find_transfer_options(
            origin=origin,
            destination=destination,
            depart_after=depart_after,
            window_seconds=window_days * SECONDS_PER_DAY,
            min_time_of_flight_seconds=min_tof_days * SECONDS_PER_DAY,
            max_time_of_flight_seconds=max_tof_days * SECONDS_PER_DAY,
            n_depart=3,
            n_tof=3,
            top_n=3,
        )

        latest_departure = result.grid.departure_times[-1]
        departure_days_elapsed = (latest_departure - depart_after).to_value("day")

        self.assertLessEqual(
            departure_days_elapsed,
            window_days + 1,
            "departure_times should stay within the requested day-scale window, "
            "not be inflated by treating seconds as days",
        )

        self.assertTrue(result.options, "expected at least one transfer option")
        for option in result.options:
            option_departure_days_elapsed = (
                option.departure_time - depart_after
            ).to_value("day")
            self.assertLessEqual(
                option_departure_days_elapsed,
                window_days + 1,
                "option departure_time (what the CLI prints) should stay within "
                "the requested day-scale window",
            )
            option_tof_days_elapsed = (
                option.arrival_time - option.departure_time
            ).to_value("day")
            self.assertLessEqual(
                option_tof_days_elapsed,
                max_tof_days + 1,
                "option arrival_time should be within the requested TOF window",
            )


if __name__ == "__main__":
    unittest.main()
