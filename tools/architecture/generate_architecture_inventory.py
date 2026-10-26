"""Generate or check the maintained architecture inventory."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

if __package__ in {None, ""}:
    # Keep direct-script execution available without importing from the test tree.
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from tools.architecture.inventory import render_inventory

OUTPUT = Path("docs/architecture_inventory.md")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--repo-root", type=Path, default=Path(__file__).resolve().parents[2]
    )
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args(argv)
    output = args.repo_root.resolve() / OUTPUT
    rendered = render_inventory(args.repo_root)
    if args.check:
        if not output.is_file() or output.read_text(encoding="utf-8") != rendered:
            print(f"architecture inventory is stale: {output}", file=sys.stderr)
            return 1
        print(f"architecture inventory is current: {output}")
        return 0
    output.write_text(rendered, encoding="utf-8", newline="\n")
    print(f"wrote {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
