"""Edges of dfn_mesh.vtu with a negative connection coefficient, by kind.

The coefficient of an edge is the sum over its triangles of
1/2 × cot(angle opposite the edge) × aperture × permeability, as in fehm-rs
(coefficients.rs). Kinds: an edge of one triangle (on the outline of the
mesh), an edge of triangles of two or more fractures (shared line), or an
edge of two triangles of one fracture. Also the nodes with a negative
volume (sum of their Voronoi areas × aperture).

Usage: uv run --no-project --with numpy --with defusedxml python -I negative_coefficients.py <output_folder> [min_feature_size]
"""
import sys
from collections import Counter, defaultdict

import numpy as np
from defusedxml import ElementTree as ET


def main(folder, min_feature_size):
    root = ET.parse(f"{folder}/dfn_mesh.vtu").getroot()
    arrays = {a.get("Name"): a.text for a in root.iter("DataArray")}
    points = np.array(arrays["points"].split(), dtype=float).reshape(-1, 3)
    triangles = np.array(arrays["connectivity"].split(), dtype=int).reshape(-1, 3)
    fracture_ids = np.array(arrays["fracture_id"].split(), dtype=int)
    apertures = np.array(arrays["aperture"].split(), dtype=float)
    permeabilities = np.array(arrays["permeability"].split(), dtype=float)
    weights = apertures * permeabilities

    corners = points[triangles]
    double_areas = np.linalg.norm(np.cross(corners[:, 1] - corners[:, 0], corners[:, 2] - corners[:, 0]), axis=1)
    coefficients = defaultdict(float)
    triangle_count = Counter()
    fractures_of_edge = defaultdict(set)
    volumes = np.zeros(len(points))
    obtuse_at = Counter()

    for corner in range(3):
        i, j, k = corner, (corner + 1) % 3, (corner + 2) % 3
        to_j = corners[:, j] - corners[:, i]
        to_k = corners[:, k] - corners[:, i]
        cotangents = np.sum(to_j * to_k, axis=1) / double_areas
        for t in range(len(triangles)):
            edge = tuple(sorted((triangles[t, j], triangles[t, k])))
            coefficients[edge] += 0.5 * cotangents[t] * weights[t]
            triangle_count[edge] += 1
            fractures_of_edge[edge].add(fracture_ids[t])
            length_squared = np.sum((corners[t, j] - corners[t, k]) ** 2)
            # Voronoi area part of the two nodes of the edge
            part = length_squared * cotangents[t] / 8 * apertures[t]
            volumes[triangles[t, j]] += part
            volumes[triangles[t, k]] += part
            if cotangents[t] < 0:
                obtuse_at[fracture_ids[t]] += 1

    negative = [edge for edge, value in coefficients.items() if value < 0]
    kinds = Counter()
    for edge in negative:
        if len(fractures_of_edge[edge]) > 1:
            kinds["shared (two or more fractures)"] += 1
        elif triangle_count[edge] == 1:
            kinds["outline (one triangle)"] += 1
        else:
            kinds["inside one fracture"] += 1
    print(f"edges: {len(coefficients)}; negative coefficients: {len(negative)}")
    for kind, count in kinds.most_common():
        print(f"  {kind}: {count}")
    print(f"nodes with a negative volume: {int(np.sum(volumes < 0))}")
    fractures_of_negative = Counter(f for edge in negative for f in fractures_of_edge[edge])
    print("fractures with most negative edges: " + ", ".join(f"{f} ({c})" for f, c in fractures_of_negative.most_common(8)))

    if min_feature_size is not None:
        slivers = sliver_fractures(folder, min_feature_size)
        on_sliver = sum(1 for edge in negative if fractures_of_edge[edge] & slivers)
        print(f"slivers (narrower than h, a corner below 10°): {len(slivers)}; "
              f"negative edges on a sliver: {on_sliver} of {len(negative)}")


def sliver_fractures(folder, min_feature_size):
    """Numbers of the fractures of fractures.vtp which are slivers"""
    root = ET.parse(f"{folder}/fractures.vtp").getroot()
    arrays = {a.get("Name"): a.text for a in root.iter("DataArray")}
    points = np.array(arrays["points"].split(), dtype=float).reshape(-1, 3)
    connectivity = np.array(arrays["connectivity"].split(), dtype=int)
    slivers = set()
    start = 0
    for number, end in enumerate(np.array(arrays["offsets"].split(), dtype=int), start=1):
        polygon = points[connectivity[start:end]]
        count = len(polygon)
        area = np.linalg.norm(sum(np.cross(polygon[i], polygon[(i + 1) % count]) for i in range(count))) / 2
        longest = max(np.linalg.norm(polygon[(i + 1) % count] - polygon[i]) for i in range(count))
        smallest_corner = min(
            np.degrees(np.arccos(np.clip(
                np.dot(polygon[(i + 1) % count] - polygon[i], polygon[i - 1] - polygon[i])
                / (np.linalg.norm(polygon[(i + 1) % count] - polygon[i]) * np.linalg.norm(polygon[i - 1] - polygon[i])),
                -1, 1)))
            for i in range(count)
        )
        if 2 * area / longest < min_feature_size and smallest_corner < 10:
            slivers.add(number)
        start = end
    return slivers


main(sys.argv[1], float(sys.argv[2]) if len(sys.argv) > 2 else None)
