#!/usr/bin/env python3
"""Thin wrapper around the standalone ``fortc`` package.

The transpiler logic now lives in ``fortc`` (a standalone Python package kept
in ``../fortc`` relative to this repo). This shim keeps the in-repo build
scripts (run_difftest.sh, run_e2e.sh, generate.sh) working whether or not
``fortc`` has been ``pip install``-ed, by falling back to the sibling checkout.

    python3 transpile.py foo.f            # prints foo.c
    python3 transpile.py foo.f --header   # prints foo.h
    python3 transpile.py foo.f --only a,b
"""
import os
import sys

try:
    from fort2c.cli import main
except ModuleNotFoundError:
    # fall back to the sibling source checkout: <repo>/../fortc
    _pkg = os.path.abspath(os.path.join(
        os.path.dirname(__file__), "..", "..", "..", "fortc"))
    sys.path.insert(0, _pkg)
    from fort2c.cli import main

if __name__ == "__main__":
    sys.exit(main())
