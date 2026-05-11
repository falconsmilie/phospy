"""Thin CLI entrypoint: parse arguments and delegate to application services."""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

from phospy.io.cli_errors import (
    CLI_EXIT_SUCCESS,
    exit_code_from_system_exit,
    present_cli_error,
)
from phospy.io.cli_parser import build_parser
from phospy.io.cli_request_factory import build_command
from phospy.io.cli_runner import CliCommandRunner


def main(argv: Sequence[str] | None = None) -> int:
    """Run the CLI."""

    parser = build_parser()
    try:
        args = parser.parse_args(list(argv) if argv is not None else None)
    except SystemExit as exc:
        return exit_code_from_system_exit(exc)
    try:
        command = build_command(args)
        result = CliCommandRunner().run(command)
    except Exception as exc:  # pragma: no cover - boundary safety
        return present_cli_error(exc)
    _print_written_summary(result.command_name, result.written)
    return CLI_EXIT_SUCCESS


def _print_written_summary(command: str, written: dict[str, Path]) -> None:
    print(f"phospy {command} completed.")
    for key in sorted(written):
        print(f"{key}: {written[key]}")
