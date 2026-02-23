# ---------------------------------------------------------------------------
# Title:        Command-Line Interface
# Author:       Dillon Breen
# Institute:    Dublin City University
# ---------------------------------------------------------------------------
#
# Single entry point for all commands.
#
# Usage:
#   uv run dm query   --type G --interval 100 --range 8000 --angle 5
#   uv run dm analyse --type G --sample 1000
#   uv run dm analyse --type all --sample all

import argparse

from .query import query_type, STELLAR_TYPES
from .visualise import stellar_distribution, rotational_curvature, dark_matter_distribution

VALID_TYPES = [*STELLAR_TYPES, "all"]


def cmd_query(args):
    """Download Gaia DR2 observations for the selected stellar type(s)."""
    targets = list(STELLAR_TYPES) if args.type == "all" else [args.type]
    for t in targets:
        query_type(t, args.interval, args.range, args.angle)
    print("\nAll queries complete.")


def cmd_analyse(args):
    """Run the full kinematic analysis and produce all three plots."""
    sample = args.sample if args.sample == "all" else int(args.sample)
    stellar_distribution(args.type, sample, fit=4)
    rotational_curvature(args.type, sample, fit=7)
    dark_matter_distribution(args.type, sample, fit=3)


def main():
    parser = argparse.ArgumentParser(
        prog="dm",
        description="Milky Way dark matter analysis using Gaia DR2 stellar kinematics.",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    # ── query ────────────────────────────────────────────────────────────────
    p_query = sub.add_parser("query", help="Download Gaia DR2 data to CSV files.")
    p_query.add_argument(
        "--type", required=True, choices=VALID_TYPES,
        metavar="TYPE", help="Stellar spectral type (O B A F G K M or all).",
    )
    p_query.add_argument(
        "--interval", required=True, type=float,
        metavar="PC", help="Shell thickness in parsec.",
    )
    p_query.add_argument(
        "--range", required=True, type=float,
        metavar="PC", help="Outer observation radius in parsec.",
    )
    p_query.add_argument(
        "--angle", required=True, type=float,
        metavar="DEG", help="Cone half-angle in degrees.",
    )
    p_query.set_defaults(func=cmd_query)

    # ── analyse ───────────────────────────────────────────────────────────────
    p_analyse = sub.add_parser("analyse", help="Compute kinematics and produce plots.")
    p_analyse.add_argument(
        "--type", required=True, choices=VALID_TYPES,
        metavar="TYPE", help="Stellar spectral type (O B A F G K M or all).",
    )
    p_analyse.add_argument(
        "--sample", default="all",
        metavar="N", help="Number of stars to sample, or 'all' (default: all).",
    )
    p_analyse.set_defaults(func=cmd_analyse)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
