"""Regression test: CLI day-denominated args must be converted to seconds."""

import unittest
from unittest import mock

import numpy as np
from astropy.time import Time

from planet_transit_calculator import cli
from planet_transit_calculator.search import PorkchopGrid, SearchResult

SECONDS_PER_DAY = 86400.0


def _fake_search_result():
    return SearchResult(
        grid=PorkchopGrid(
            departure_times=Time(["2026-07-01"]),
            tofs_days=np.array([100.0]),
            total_dv=np.array([[1.0]]),
        ),
        options=[],
    )


class TestCliUnitConversion(unittest.TestCase):
    def test_day_args_are_converted_to_seconds_before_forwarding(self):
        captured = {}

        def fake_find_transfer_options(**kwargs):
            captured.update(kwargs)
            return _fake_search_result()

        with mock.patch(
            "planet_transit_calculator.cli.find_transfer_options",
            side_effect=fake_find_transfer_options,
        ):
            cli.main(
                [
                    "--from", "earth",
                    "--to", "mars",
                    "--window-days", "30",
                    "--min-tof", "100",
                    "--max-tof", "300",
                ]
            )

        self.assertEqual(captured["window_seconds"], 30 * SECONDS_PER_DAY)
        self.assertEqual(captured["min_time_of_flight_seconds"], 100 * SECONDS_PER_DAY)
        self.assertEqual(captured["max_time_of_flight_seconds"], 300 * SECONDS_PER_DAY)


if __name__ == "__main__":
    unittest.main()
