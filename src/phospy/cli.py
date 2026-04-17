"""CLI placeholder for the rewrite package boundary."""

from __future__ import annotations


def main() -> int:
    """Run the rewrite-phase CLI entrypoint."""
    print(
        "phospy CLI is in rewrite mode. "
        "Use the Python API under src/phospy/ while workflows are rebuilt."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
