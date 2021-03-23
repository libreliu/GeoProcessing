#!/usr/bin/env python3

import sys, importlib

if len(sys.argv) < 2:
    print(f"Usage: {sys.argv[0]} package_name")
    sys.exit(1)

pkg_main = importlib.import_module(f'{sys.argv[1]}.main', )
pkg_main.main(sys.argv[2:])