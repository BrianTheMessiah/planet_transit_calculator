"""Command-line interface for the planet transit calculator."""

from __future__ import annotations

import argparse
from typing import Sequence

from astropy.time import Time

from .bodies import CELESTIAL_BODIES, get_celestial_body
from .search import find_transfer_options
from .transfer import SECONDS_PER_DAY


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="planet-transit",
        description=(
            "Find interplanetary transfer options between two solar system bodies, "
            "using real-time planetary ephemeris data."
        ),
    )
    body_names = ", ".join(sorted(CELESTIAL_BODIES))
    parser.add_argument("--from", dest="origin", required=True, help=f"Departure body ({body_names})")
    parser.add_argument("--to", dest="destination", required=True, help=f"Arrival body ({body_names})")
    parser.add_argument(
        "--depart-after",
        default=None,
        help="Earliest departure date, YYYY-MM-DD (default: now)",
    )
    parser.add_argument(
        "--window-days",
        type=float,
        default=None,
        help="Length of the departure-date search window in days (default: synodic period)",
    )
    parser.add_argument("--min-tof", type=float, default=None, help="Minimum time of flight in days")
    parser.add_argument("--max-tof", type=float, default=None, help="Maximum time of flight in days")
    parser.add_argument("--top", type=int, default=5, help="Number of options to display (default: 5)")
    parser.add_argument("--plot", default=None, metavar="FILE", help="Save a porkchop plot (PNG) to FILE")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    try:
        origin = get_celestial_body(args.origin)
        destination = get_celestial_body(args.destination)
    except ValueError as exc:
        parser.error(str(exc))

    if origin.name == destination.name:
        parser.error("Departure and arrival bodies must be different")

    depart_after = Time(args.depart_after) if args.depart_after else Time.now()

    window_seconds = (
        args.window_days * SECONDS_PER_DAY if args.window_days is not None else None
    )
    min_tof_seconds = (
        args.min_tof * SECONDS_PER_DAY if args.min_tof is not None else None
    )
    max_tof_seconds = (
        args.max_tof * SECONDS_PER_DAY if args.max_tof is not None else None
    )

    result = find_transfer_options(
        origin=origin,
        destination=destination,
        depart_after=depart_after,
        window_seconds=window_seconds,
        min_time_of_flight_seconds=min_tof_seconds,
        max_time_of_flight_seconds=max_tof_seconds,
        top_n=args.top,
    )

    if not result.options:
        print("No transfer solutions found in the search range.")
        return 1

    print(f"Transfer options: {origin.name.title()} -> {destination.name.title()}")
    print(f"Departure search window starts {depart_after.iso[:10]}\n")
    header = (
        f"{'#':>2}  {'Depart':<10}  {'Arrive':<10}  {'TOF (days)':>10}  "
        f"{'Depart dV':>10}  {'Arrival dV':>10}  {'Total dV':>10}"
    )
    print(header)
    print("-" * len(header))
    for idx, opt in enumerate(result.options, start=1):
        print(
            f"{idx:>2}  {opt.departure_time.iso[:10]:<10}  {opt.arrival_time.iso[:10]:<10}  "
            f"{opt.tof_days:>10.0f}  {opt.departure_dv:>9.2f}   {opt.arrival_dv:>9.2f}   "
            f"{opt.total_dv:>9.2f}"
        )
    print("\n(dV values in km/s)")

    if args.plot:
        from .plotting import save_porkchop

        save_porkchop(result, args.plot)
        print(f"Porkchop plot saved to {args.plot}")

    return 0
