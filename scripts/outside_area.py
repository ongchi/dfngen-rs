"""Part of the triangle mesh outside the domain mesh, by the vertical extent.

Reports the bounds of the domain mesh (.vtu) and the part of the triangle
area of dfn_mesh.vtu above the top and below the bottom of the domain. The
domain can have a surface which is not flat: the check uses the bounds only.

Usage: uv run --no-project --with numpy --with defusedxml python -I outside_area.py <domain.vtu> <output_folder>
"""
import sys

import numpy as np
from defusedxml import ElementTree as ET


def points_of(path):
    root = ET.parse(path).getroot()
    points_element = next(root.iter("Points"))
    data_array = next(points_element.iter("DataArray"))
    return np.array(data_array.text.split(), dtype=float).reshape(-1, 3)


def main(domain_path, folder):
    domain = points_of(domain_path)
    low, high = domain.min(axis=0), domain.max(axis=0)
    print(f"domain bounds: {low} to {high}")

    root = ET.parse(f"{folder}/dfn_mesh.vtu").getroot()
    arrays = {a.get("Name"): a for a in root.iter("DataArray")}
    points = np.array(arrays["points"].text.split(), dtype=float).reshape(-1, 3)
    triangles = np.array(arrays["connectivity"].text.split(), dtype=int).reshape(-1, 3)
    corners = points[triangles]
    areas = np.linalg.norm(np.cross(corners[:, 1] - corners[:, 0], corners[:, 2] - corners[:, 0]), axis=1) / 2
    centroids = corners.mean(axis=1)
    print(f"mesh bounds: {points.min(axis=0)} to {points.max(axis=0)}")

    outside = np.any((centroids < low) | (centroids > high), axis=1)
    print(f"triangles outside the domain bounds: {np.mean(outside) * 100:.1f} % of the triangles, "
          f"{areas[outside].sum() / areas.sum() * 100:.1f} % of the area")
    above = centroids[:, 2] > high[2]
    print(f"above the top of the bounds: {areas[above].sum() / areas.sum() * 100:.1f} % of the area")


main(sys.argv[1], sys.argv[2])
