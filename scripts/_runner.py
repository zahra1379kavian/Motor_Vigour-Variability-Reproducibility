#!/usr/bin/env python3
"""Small helpers for running copied analysis modules from repo-level scripts."""

from __future__ import annotations

import importlib
import runpy
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
MODULE_DIR = REPO_ROOT / "analysis" / "modules"
GVS_DIR = REPO_ROOT / "analysis" / "gvs_connectivity_coactivation"


def _ensure_paths() -> None:
    for path in (MODULE_DIR, GVS_DIR, REPO_ROOT):
        if str(path) not in sys.path:
            sys.path.insert(0, str(path))


def _with_argv(argv: list[str] | None):
    class _ArgvContext:
        def __enter__(self):
            self.previous = sys.argv[:]
            if argv is not None:
                sys.argv = argv[:]

        def __exit__(self, exc_type, exc, tb):
            sys.argv = self.previous

    return _ArgvContext()


def run_module(module_name: str, argv: list[str] | None = None) -> int | None:
    _ensure_paths()
    with _with_argv(argv):
        module = importlib.import_module(module_name)
        return module.main()


def run_script(relative_path: str, argv: list[str] | None = None) -> None:
    _ensure_paths()
    with _with_argv(argv):
        runpy.run_path(str(REPO_ROOT / relative_path), run_name="__main__")
