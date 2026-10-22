"""Fractures of fractures.vtp which are thin, and their part of the bad triangles.

The width of a polygon is 2 × area / longest edge. Prints the number of
fractures with a width below each limit (as a part of min_feature_size) and
the part of the triangles of dfn_mesh.vtu with an angle below 10° which they
have.

Usage: uv run --no-project --with numpy --with defusedxml python -I thin_pieces.py <output_folder> <min_feature_size>
"""
import sys

import numpy as np
from defusedxml import ElementTree as ET


def arrays(path):
    root = ET.parse(path).getroot()
    return {a.get("Name"): np.array(a.text.split(), dtype=float) for a in root.iter("DataArray")}


def main(folder, min_feature_size):
    fractures = arrays(f"{folder}/fractures.vtp")
    points = fractures["points"].reshape(-1, 3)
    connectivity = fractures["connectivity"].astype(int)
    widths = []
    corner_angles = []
    start = 0
    for end in fractures["offsets"].astype(int):
        polygon = points[connectivity[start:end]]
        count = len(polygon)
        cross = sum(np.cross(polygon[i], polygon[(i + 1) % count]) for i in range(count))
        area = np.linalg.norm(cross) / 2
        longest = max(np.linalg.norm(polygon[(i + 1) % count] - polygon[i]) for i in range(count))
        widths.append(2 * area / longest)
        angles = []
        for i in range(count):
            a = polygon[(i + 1) % count] - polygon[i]
            b = polygon[i - 1] - polygon[i]
            cosine = np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b))
            angles.append(np.degrees(np.arccos(np.clip(cosine, -1, 1))))
        corner_angles.append(min(angles))
        start = end
    widths = np.array(widths)
    corner_angles = np.array(corner_angles)

    mesh = arrays(f"{folder}/dfn_mesh.vtu")
    mesh_points = mesh["points"].reshape(-1, 3)
    triangles = mesh["connectivity"].astype(int).reshape(-1, 3)
    fracture_ids = mesh["fracture_id"].astype(int)
    corners = mesh_points[triangles]
    smallest = np.full(len(triangles), 180.0)
    for i in range(3):
        a = corners[:, (i + 1) % 3] - corners[:, i]
        b = corners[:, (i + 2) % 3] - corners[:, i]
        cosine = np.sum(a * b, axis=1) / (np.linalg.norm(a, axis=1) * np.linalg.norm(b, axis=1))
        smallest = np.minimum(smallest, np.degrees(np.arccos(np.clip(cosine, -1, 1))))
    is_bad = smallest < 10

    for name, selection in (
        ("thin, smallest corner below 10°", (widths < min_feature_size) & (corner_angles < 10)),
        ("thin, smallest corner 10° or more", (widths < min_feature_size) & (corner_angles >= 10)),
        ("not thin, smallest corner below 10°", (widths >= min_feature_size) & (corner_angles < 10)),
    ):
        in_selection = selection[fracture_ids - 1]
        print(f"{name}: {selection.sum()} fractures; bad triangles in them: "
              f"{np.sum(is_bad & in_selection)} of {is_bad.sum()} of their {in_selection.sum()} triangles")

    for limit in (0.1, 0.5, 1.0):
        is_thin = widths < limit * min_feature_size
        in_thin = is_thin[fracture_ids - 1]
        print(f"width below {limit} h: {is_thin.sum()} fractures; their triangles: {in_thin.sum()} "
              f"({np.mean(in_thin) * 100:.2f} % of all); bad triangles in them: {np.sum(is_bad & in_thin)} "
              f"of {is_bad.sum()} ({np.sum(is_bad & in_thin) / max(is_bad.sum(), 1) * 100:.0f} %)")


main(sys.argv[1], float(sys.argv[2]))
