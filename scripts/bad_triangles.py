"""Where the triangles with a small angle of one fracture are.

Prints, for the triangles of the fracture with an angle below the limit:
their height (z of the centroid) in bins, their shortest edge, and how many
of their points are on a boundary array.

Usage: uv run --no-project --with numpy --with defusedxml python -I bad_triangles.py <output_folder> <fracture_id> [limit_degrees]
"""
import sys

import numpy as np
from defusedxml import ElementTree as ET


def main(folder, fracture, limit):
    root = ET.parse(f"{folder}/dfn_mesh.vtu").getroot()
    arrays = {a.get("Name"): np.array(a.text.split(), dtype=float) for a in root.iter("DataArray")}
    points = arrays["points"].reshape(-1, 3)
    triangles = arrays["connectivity"].astype(int).reshape(-1, 3)
    mine = triangles[arrays["fracture_id"].astype(int) == fracture]
    on_boundary = np.zeros(len(points), dtype=bool)
    for name, values in arrays.items():
        if name.startswith("boundary__"):
            on_boundary |= values > 0

    corners = points[mine]
    edges = [corners[:, (i + 1) % 3] - corners[:, i] for i in range(3)]
    lengths = np.stack([np.linalg.norm(edge, axis=1) for edge in edges], axis=1)
    angles = []
    for i in range(3):
        a, b = -edges[(i + 2) % 3], edges[i]
        cosine = np.sum(a * b, axis=1) / (np.linalg.norm(a, axis=1) * np.linalg.norm(b, axis=1))
        angles.append(np.degrees(np.arccos(np.clip(cosine, -1, 1))))
    smallest = np.stack(angles, axis=1).min(axis=1)
    bad = smallest < limit

    print(f"fracture {fracture}: {len(mine)} triangles, {bad.sum()} below {limit}°")
    if not bad.any():
        return
    z = corners.mean(axis=1)[:, 2]
    print(f"  z of all triangles: {z.min():.0f} to {z.max():.0f}")
    counts, bins = np.histogram(z[bad], bins=8, range=(z.min(), z.max()))
    print("  bad triangles by z: " + ", ".join(f"{bins[i]:.0f}: {c}" for i, c in enumerate(counts)))
    print(f"  shortest edge of bad triangles: quantiles 10/50/90 %: {np.quantile(lengths[bad].min(axis=1), [0.1, 0.5, 0.9])}")
    print(f"  longest edge of bad triangles: quantiles 10/50/90 %: {np.quantile(lengths[bad].max(axis=1), [0.1, 0.5, 0.9])}")
    boundary_counts = on_boundary[mine[bad]].sum(axis=1)
    print(f"  boundary points per bad triangle (0, 1, 2, 3): {np.bincount(boundary_counts, minlength=4)}")


main(sys.argv[1], int(sys.argv[2]), float(sys.argv[3]) if len(sys.argv) > 3 else 10.0)
