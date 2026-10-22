"""Lengths of the boundary line segments of fracture_boundaries.vtp.

Usage: uv run --no-project --with numpy --with defusedxml python -I boundary_segments.py <output_folder> <min_feature_size>
"""
import sys

import numpy as np
from defusedxml import ElementTree as ET


def main(folder, min_feature_size):
    root = ET.parse(f"{folder}/fracture_boundaries.vtp").getroot()
    arrays = {a.get("Name"): a for a in root.iter("DataArray")}
    points = np.array(arrays["points"].text.split(), dtype=float).reshape(-1, 2, 3)
    lengths = np.linalg.norm(points[:, 1] - points[:, 0], axis=1)
    mesh_surface = np.array(arrays["mesh_surface"].text.split(), dtype=int) > 0

    print(f"segments: {len(lengths)}, total length {lengths.sum():.0f} m")
    for limit in (0.1, 0.5, 1.0):
        short = lengths < limit * min_feature_size
        print(f"  shorter than {limit} h: {short.sum()} ({np.mean(short) * 100:.1f} %); "
              f"on the mesh surface: {np.sum(short & mesh_surface)}")
    print(f"  median length {np.median(lengths):.2f} m")


main(sys.argv[1], float(sys.argv[2]))
