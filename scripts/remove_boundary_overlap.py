"""Test only: sets boundary__<first> to 0 at the points where boundary__<second> is 1.

Usage: uv run --no-project --with defusedxml python -I remove_boundary_overlap.py <in.vtu> <out.vtu> <first> <second>
"""
import sys

from defusedxml import ElementTree as ET

source, target, first, second = sys.argv[1:5]
tree = ET.parse(source)
arrays = {array.get("Name"): array for array in tree.getroot().iter("DataArray")}
first_array = arrays[f"boundary__{first}"]
second_values = arrays[f"boundary__{second}"].text.split()
first_values = first_array.text.split()
changed = 0
for i, value in enumerate(second_values):
    if value != "0" and first_values[i] != "0":
        first_values[i] = "0"
        changed += 1
first_array.text = "\n" + "\n".join(first_values) + "\n"
tree.write(target)
print(f"{changed} points: {first} set to 0")
