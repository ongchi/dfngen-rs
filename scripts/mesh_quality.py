"""Triangle quality of a dfn_mesh.vtu of dfngen, by kind of fracture.

Usage: uv run --no-project --with numpy --with defusedxml python -I mesh_quality.py <output_folder>
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


def triangle_angles(points, triangles):
    corners = points[triangles]
    angles = []
    for i in range(3):
        a = corners[:, (i + 1) % 3] - corners[:, i]
        b = corners[:, (i + 2) % 3] - corners[:, i]
        cosine = np.sum(a * b, axis=1) / (np.linalg.norm(a, axis=1) * np.linalg.norm(b, axis=1))
        angles.append(np.degrees(np.arccos(np.clip(cosine, -1, 1))))
    return np.stack(angles, axis=1)


def report(name, angles):
    smallest = angles.min(axis=1)
    largest = angles.max(axis=1)
    print(
        f"{name}: {len(angles)} triangles; smallest angle < 10°: {np.mean(smallest < 10) * 100:.1f} %, "
        f"< 20°: {np.mean(smallest < 20) * 100:.1f} %, < 30°: {np.mean(smallest < 30) * 100:.1f} %; "
        f"largest angle > 120°: {np.mean(largest > 120) * 100:.1f} %, > 90°: {np.mean(largest > 90) * 100:.1f} %"
    )


def main(folder):
    mesh = arrays(f"{folder}/dfn_mesh.vtu")
    points = mesh["points"].reshape(-1, 3)
    triangles = mesh["connectivity"].astype(int).reshape(-1, 3)
    fracture_ids = mesh["fracture_id"].astype(int)
    families = arrays(f"{folder}/fractures.vtp")["family"].astype(int)
    angles = triangle_angles(points, triangles)

    report("all", angles)
    is_user = families[fracture_ids - 1] == 0
    if is_user.any():
        report("user-defined", angles[is_user])
    if (~is_user).any():
        report("stochastic", angles[~is_user])

    corners = points[triangles]
    lengths = np.concatenate([np.linalg.norm(corners[:, (i + 1) % 3] - corners[:, i], axis=1) for i in range(3)])
    quantiles = np.quantile(lengths, [0.01, 0.1, 0.5, 0.9, 0.99])
    print("edge length quantiles 1 %, 10 %, 50 %, 90 %, 99 %: " + ", ".join(f"{q:.3g}" for q in quantiles))

    # Fractures with the most bad triangles
    is_bad = angles.min(axis=1) < 10
    bad_counts = np.bincount(fracture_ids[is_bad], minlength=len(families) + 1)
    totals = np.bincount(fracture_ids, minlength=len(families) + 1)
    worst = np.argsort(bad_counts)[::-1][:8]
    print("fractures with most triangles below 10°: (fracture, family, bad, all)")
    for fracture in worst:
        print(f"  {fracture} {families[fracture - 1]} {bad_counts[fracture]} {totals[fracture]}")


main(sys.argv[1])
