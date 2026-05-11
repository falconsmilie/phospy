"""CLI error presentation and exit code mapping."""

from __future__ import annotations

import sys
from typing import TextIO

from phospy.errors import PhosPyError

CLI_EXIT_SUCCESS = 0
CLI_EXIT_INTERNAL_ERROR = 1
CLI_EXIT_USER_ERROR = 2


def exit_code_from_system_exit(exc: SystemExit) -> int:
    """Normalize argparse/system-exit termination into explicit CLI exit codes."""

    code = exc.code
    if code is None:
        return CLI_EXIT_SUCCESS
    if isinstance(code, int):
        return code
    return CLI_EXIT_USER_ERROR


def present_cli_error(exc: Exception, *, stream: TextIO | None = None) -> int:
    """Render an exception into user-facing CLI stderr output and exit code."""

    target = sys.stderr if stream is None else stream
    if isinstance(exc, PhosPyError):
        print(f"{exc.__class__.__name__}: {exc}", file=target)
        return CLI_EXIT_USER_ERROR
    print(f"UnhandledError: {exc}", file=target)
    return CLI_EXIT_INTERNAL_ERROR
