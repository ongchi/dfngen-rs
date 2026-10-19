"""Checks a dfn_mesh.vtu of dfngen: angles, areas, duplicate points, connectivity.

Usage: uv run --no-project --with numpy --with defusedxml python -I check_mesh.py <output_folder>
"""
import sys
from collections import defaultdict

import numpy as np
from defusedxml import ElementTree as ET


def arrays(path):
    root = ET.parse(path).getroot()
    result = {}
    for data_array in root.iter("DataArray"):
        name = data_array.get("Name")
        result[name] = np.array(data_array.text.split(), dtype=float)
    return result


def polygon_areas(path):
    data = arrays(path)
    points = data["points"].reshape(-1, 3)
    offsets = data["offsets"].astype(int)
    connectivity = data["connectivity"].astype(int)
    areas = []
    start = 0
    for end in offsets:
        polygon = points[connectivity[start:end]]
        cross = np.zeros(3)
        for i in range(len(polygon)):
            cross += np.cross(polygon[i], polygon[(i + 1) % len(polygon)])
        areas.append(np.linalg.norm(cross) / 2)
        start = end
    return np.array(areas)


def main(folder):
    data = arrays(f"{folder}/dfn_mesh.vtu")
    points = data["points"].reshape(-1, 3)
    triangles = data["connectivity"].astype(int).reshape(-1, 3)
    fractures = data["fracture_id"].astype(int)
    print(f"points {len(points)}, triangles {len(triangles)}")

    corners = points[triangles]
    edges = [corners[:, (i + 1) % 3] - corners[:, i] for i in range(3)]
    lengths = np.stack([np.linalg.norm(edge, axis=1) for edge in edges], axis=1)
    angles = []
    for i in range(3):
        a = -edges[(i + 2) % 3]
        b = edges[i]
        cosine = np.sum(a * b, axis=1) / (np.linalg.norm(a, axis=1) * np.linalg.norm(b, axis=1))
        angles.append(np.degrees(np.arccos(np.clip(cosine, -1, 1))))
    angles = np.stack(angles, axis=1)
    print(f"min angle {angles.min():.2f}°, 1% quantile of min angles {np.quantile(angles.min(axis=1), 0.01):.2f}°")
    print(f"obtuse triangles {np.sum(angles.max(axis=1) > 90)} ({100 * np.mean(angles.max(axis=1) > 90):.1f} %)")
    print(f"edge length min {lengths.min():.4g}, max {lengths.max():.4g}")

    triangle_areas = np.linalg.norm(np.cross(edges[0], -edges[2]), axis=1) / 2
    mesh_areas = np.bincount(fractures - 1, weights=triangle_areas)
    fracture_areas = polygon_areas(f"{folder}/fractures.vtp")
    relative = np.abs(mesh_areas - fracture_areas) / fracture_areas
    print(f"area: largest relative difference {relative.max():.2e} (fracture {relative.argmax() + 1})")

    rounded = np.round(points / (1e-9 * np.ptp(points)))
    _, counts = np.unique(rounded, axis=0, return_counts=True)
    print(f"points at the same position: {np.sum(counts > 1)}")

    parent = list(range(len(points)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for a, b, c in triangles:
        for x, y in ((a, b), (b, c)):
            parent[find(x)] = find(y)
    components = defaultdict(set)
    for triangle, fracture in zip(triangles, fractures):
        components[find(triangle[0])].add(fracture)
    print(f"connected parts of the mesh: {len(components)}, fractures {len(set(fractures))}")
    sizes = sorted((len(fracture_set) for fracture_set in components.values()), reverse=True)
    print(f"fractures in the parts: {sizes[:10]}")

    neighbors = [list(map(int, line.split())) for line in open(f"{folder}/connectivity.dat")]
    cluster = list(range(len(neighbors) + 1))

    def find_cluster(i):
        while cluster[i] != i:
            i = cluster[i]
        return i

    for fracture, others in enumerate(neighbors, start=1):
        for other in others:
            cluster[find_cluster(fracture)] = find_cluster(other)
    # The pieces of a user-defined fracture are one fracture
    user_ids = arrays(f"{folder}/fractures.vtp")["user_fracture_id"].astype(int)
    first_piece = {}
    for fracture, user_id in enumerate(user_ids, start=1):
        if user_id >= 0:
            if user_id in first_piece:
                cluster[find_cluster(fracture)] = find_cluster(first_piece[user_id])
            else:
                first_piece[user_id] = fracture
    expected = {frozenset(f for f in range(1, len(neighbors) + 1) if find_cluster(f) == root)
                for root in {find_cluster(f) for f in range(1, len(neighbors) + 1)}}
    actual = {frozenset(fracture_set) for fracture_set in components.values()}
    print(f"mesh parts equal the clusters of connectivity.dat: {expected == actual}"
          f" ({len(expected)} clusters)")
    if expected != actual:
        split = [sorted(group) for group in expected if group not in actual]
        print(f"clusters which are not one mesh part: {len(split)}; sizes {sorted(map(len, split), reverse=True)[:10]}")
        for group in split[:3]:
            parts = sorted((sorted(part & set(group)) for part in actual if part & set(group)), key=len)
            print(f"  cluster of {len(group)} in {len(parts)} parts: {[part[:12] for part in parts[:5]]}")

    for name in data:
        if name.startswith("boundary__"):
            print(f"{name}: {int(np.sum(data[name] > 0))} points")


main(sys.argv[1])
