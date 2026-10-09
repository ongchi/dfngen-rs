#!/usr/bin/env python3
"""Export the dfnGen input files of a pydfnworks example without running dfnGen.

Usage: export_dfngen_input.py <example_dir> <work_dir> [driver_script]

Copies the example to <work_dir>, runs its driver script there, and stops the
driver when it starts the network generation. At that point pydfnworks has
written the dfnGen input file (`*_clean.dat`) and the user-defined fracture
files below <work_dir>. The script prints the path of the input file on a line
which starts with "dfnGen input file:".

`driver_script` is the file name of the driver in the example folder. It is
necessary only if the folder contains more than one Python script.

Requires pydfnworks, for example in a virtual environment:

    uv venv .venv
    uv pip install --python .venv/bin/python <dfnWorks>/pydfnworks
    .venv/bin/python scripts/export_dfngen_input.py <dfnWorks>/examples/TPL /tmp/TPL

The exported file contains absolute paths to the user-defined fracture files.

dfngen-rs reads TOML. Convert the exported file with:

    python3 scripts/convert_dat_to_toml.py <work_dir>/..._clean.dat <name>.toml
"""

from __future__ import annotations

import os
import runpy
import shutil
import sys
from pathlib import Path

from pydfnworks import DFNWORKS


class NetworkGenerationReached(Exception):
    pass


def stop_driver(self, *args, **kwargs):
    raise NetworkGenerationReached


def find_driver(driver_name: str | None) -> Path:
    if driver_name:
        return Path(driver_name)

    drivers = sorted(Path(".").glob("*.py"))
    if len(drivers) != 1:
        names = ", ".join(str(driver) for driver in drivers)
        sys.exit(f"expected one driver script, found: {names}. Give the driver script name.")

    return drivers[0]


def main() -> int:
    if len(sys.argv) not in (3, 4):
        print(__doc__, file=sys.stderr)
        return 2

    example_dir, work_dir = Path(sys.argv[1]), Path(sys.argv[2]).resolve()
    driver_name = sys.argv[3] if len(sys.argv) == 4 else None

    shutil.rmtree(work_dir, ignore_errors=True)
    shutil.copytree(example_dir, work_dir)
    os.chdir(work_dir)

    # The dfnWorks executables and the flow and transport input files are not
    # necessary to write the dfnGen input files
    DFNWORKS.define_paths = lambda self: None
    DFNWORKS.check_input_paths = lambda self: None
    DFNWORKS.create_network = stop_driver

    try:
        runpy.run_path(str(find_driver(driver_name)), run_name="__main__")
    except NetworkGenerationReached:
        input_files = sorted(work_dir.rglob("*_clean.dat"))
        if len(input_files) != 1:
            print(f"expected one dfnGen input file, found {len(input_files)}", file=sys.stderr)
            return 1

        print(f"dfnGen input file: {input_files[0]}")
        return 0

    print("the driver ended before network generation", file=sys.stderr)
    return 1


if __name__ == "__main__":
    sys.exit(main())
