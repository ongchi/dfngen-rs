"""Triangles of dfn_mesh.vtu which fehm-rs rejects as having no area.

fehm-rs (coefficients::cotangents) rejects a triangle with
2 × area <= 1e-12 × (longest edge)². Prints their number, the fractures
with most of them, and the first few with their edge lengths.

Usage: uv run --no-project --with numpy --with defusedxml python -I degenerate_triangles.py <output_folder> [ratio]
"""
import sys

import numpy as np
from defusedxml import ElementTree as ET


def main(folder, ratio):
    root = ET.parse(f"{folder}/dfn_mesh.vtu").getroot()
    arrays = {a.get("Name"): a.text for a in root.iter("DataArray")}
    points = np.array(arrays["points"].split(), dtype=float).reshape(-1, 3)
    triangles = np.array(arrays["connectivity"].split(), dtype=int).reshape(-1, 3)
    fracture_ids = np.array(arrays["fracture_id"].split(), dtype=int)

    corners = points[triangles]
    double_areas = np.linalg.norm(np.cross(corners[:, 1] - corners[:, 0], corners[:, 2] - corners[:, 0]), axis=1)
    lengths = np.stack([np.linalg.norm(corners[:, (i + 1) % 3] - corners[:, i], axis=1) for i in range(3)], axis=1)
    longest_squared = lengths.max(axis=1) ** 2
    is_degenerate = double_areas <= ratio * longest_squared

    print(f"triangles: {len(triangles)}; degenerate (2 A <= {ratio} L²): {is_degenerate.sum()}")
    if not is_degenerate.any():
        return
    unique, counts = np.unique(fracture_ids[is_degenerate], return_counts=True)
    order = np.argsort(counts)[::-1][:10]
    print("fractures with most of them: " + ", ".join(f"{unique[i]} ({counts[i]})" for i in order))
    for cell in np.flatnonzero(is_degenerate)[:5]:
        print(f"  cell {cell}: fracture {fracture_ids[cell]}, points {triangles[cell]}, edges {lengths[cell]}")


main(sys.argv[1], float(sys.argv[2]) if len(sys.argv) > 2 else 1e-12)
