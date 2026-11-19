#!/usr/bin/env python3
"""Report the build configuration of one OpenBLAS shared library."""
import argparse
import ctypes
import json
from pathlib import Path
import re


def inspect(path):
    path = Path(path).resolve(strict=True)
    lib = ctypes.CDLL(str(path))
    lib.openblas_get_config.argtypes = []
    lib.openblas_get_config.restype = ctypes.c_char_p
    lib.openblas_get_parallel.argtypes = []
    lib.openblas_get_parallel.restype = ctypes.c_int
    config = lib.openblas_get_config().decode()
    mode = lib.openblas_get_parallel()
    match = re.search(r"MAX_THREADS=(\d+)", config)
    return dict(library=str(path), configuration=config, threading_mode=mode,
        threading={0: "serial", 1: "pthreads", 2: "OpenMP"}.get(mode, "unknown"),
        uses_64bit_blas_integers="USE64BITINT" in config,
        build_max_threads=int(match.group(1)) if match else None)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("library", type=Path)
    args = parser.parse_args()
    try:
        print(json.dumps(inspect(args.library), indent=2))
    except (OSError, AttributeError) as exc:
        parser.exit(1, "OpenBLAS inspection failed: " + str(exc) + "\n")
