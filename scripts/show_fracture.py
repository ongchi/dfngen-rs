"""Shape of one fracture of fractures.vtp and its intersections.

Usage: uv run --no-project --with numpy --with defusedxml python -I show_fracture.py <output_folder> <fracture_id>
"""
import sys

import numpy as np
from defusedxml import ElementTree as ET


def arrays(path):
    root = ET.parse(path).getroot()
    return {
        data_array.get("Name"): np.array(data_array.text.split(), dtype=float)
        for data_array in root.iter("DataArray")
    }


def main(folder, fracture):
    data = arrays(f"{folder}/fractures.vtp")
    points = data["points"].reshape(-1, 3)
    offsets = data["offsets"].astype(int)
    start = 0 if fracture == 1 else offsets[fracture - 2]
    vertices = points[data["connectivity"].astype(int)[start:offsets[fracture - 1]]]

    print(f"fracture {fracture}: user_fracture_id {int(data['user_fracture_id'][fracture - 1])}")
    for i, vertex in enumerate(vertices):
        edge = np.linalg.norm(vertices[(i + 1) % len(vertices)] - vertex)
        print(f"  vertex {vertex}, edge to next {edge:.2f} m")

    lines = [line.split() for line in open(f"{folder}/intersection_list.dat")][1:]
    mine = [line for line in lines if int(line[0]) == fracture or int(line[1]) == fracture]
    lengths = sorted(float(line[5]) for line in mine)
    print(f"  intersections: {len(mine)}; lengths: shortest {lengths[:5]}, longest {lengths[-3:]}")


main(sys.argv[1], int(sys.argv[2]))
