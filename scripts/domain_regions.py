"""Cells of a domain mesh (.vtu) by their regions (cell arrays region__<name>).

Usage: uv run --no-project --with numpy --with defusedxml python -I domain_regions.py <domain.vtu>
"""
import sys
from collections import Counter

import numpy as np
from defusedxml import ElementTree as ET


def main(path):
    root = ET.parse(path).getroot()
    cell_data = next(root.iter("CellData"))
    regions = {
        array.get("Name")[len("region__"):]: np.array(array.text.split(), dtype=float) > 0
        for array in cell_data.iter("DataArray")
        if array.get("Name", "").startswith("region__")
    }
    names = list(regions)
    cell_count = len(next(iter(regions.values())))
    combinations = Counter(
        tuple(name for name in names if regions[name][cell]) for cell in range(cell_count)
    )
    print(f"cells: {cell_count}; regions: {names}")
    for combination, count in combinations.most_common():
        print(f"  {combination or '(no region)'}: {count}")


main(sys.argv[1])
