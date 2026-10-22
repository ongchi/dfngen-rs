"""Independent check: triangles of dfn_mesh.vtu with their centroid outside the domain mesh.

The cells of the domain mesh (.vtu) are split into tetrahedra (VTK types 10
tetrahedron, 12 hexahedron, 13 wedge, 14 pyramid), and the centroid of each
sampled triangle is tested against the nearest tetrahedra (k-d tree of the
tetrahedron centroids) with barycentric coordinates. The output coordinates
of dfn_mesh.vtu are those of the domain mesh.

Usage: uv run --no-project --with numpy --with scipy --with defusedxml python -I triangles_outside.py <domain.vtu> <output_folder> [sample_size]
"""
import sys

import numpy as np
from defusedxml import ElementTree as ET
from scipy.spatial import cKDTree

# Corner indices of the tetrahedra of each cell type
SPLITS = {
    10: [(0, 1, 2, 3)],
    12: [(0, 1, 3, 4), (1, 2, 3, 6), (1, 4, 5, 6), (3, 4, 6, 7), (1, 3, 4, 6)],
    13: [(0, 1, 2, 3), (1, 2, 3, 4), (2, 3, 4, 5)],
    14: [(0, 1, 2, 4), (0, 2, 3, 4)],
}


def arrays(path):
    root = ET.parse(path).getroot()
    return {a.get("Name"): a.text for a in root.iter("DataArray")}, root


def domain_tetrahedra(path):
    data, root = arrays(path)
    points_element = next(root.iter("Points"))
    points = np.array(next(points_element.iter("DataArray")).text.split(), dtype=float).reshape(-1, 3)
    connectivity = np.array(data["connectivity"].split(), dtype=int)
    offsets = np.array(data["offsets"].split(), dtype=int)
    types = np.array(data["types"].split(), dtype=int)
    tetrahedra = []
    start = 0
    for end, cell_type in zip(offsets, types):
        corners = connectivity[start:end]
        for split in SPLITS[cell_type]:
            tetrahedra.append(corners[list(split)])
        start = end
    return points[np.array(tetrahedra)]


def contains(tetrahedron, point, tolerance=1e-9):
    a, b, c, d = tetrahedron
    matrix = np.column_stack([b - a, c - a, d - a])
    try:
        weights = np.linalg.solve(matrix, point - a)
    except np.linalg.LinAlgError:
        return False
    return weights.min() >= -tolerance and weights.sum() <= 1 + tolerance


def main(domain_path, folder, sample_size):
    tetrahedra = domain_tetrahedra(domain_path)
    tree = cKDTree(tetrahedra.mean(axis=1))

    data, _ = arrays(f"{folder}/dfn_mesh.vtu")
    points = np.array(data["points"].split(), dtype=float).reshape(-1, 3)
    triangles = np.array(data["connectivity"].split(), dtype=int).reshape(-1, 3)
    fracture_ids = np.array(data["fracture_id"].split(), dtype=int)

    generator = np.random.default_rng(1)
    sample = generator.choice(len(triangles), size=min(sample_size, len(triangles)), replace=False)
    centroids = points[triangles[sample]].mean(axis=1)
    _, nearest = tree.query(centroids, k=64)

    outside = [
        index
        for index, centroid, candidates in zip(sample, centroids, nearest)
        if not any(contains(tetrahedra[t], centroid) for t in candidates)
    ]
    print(f"tetrahedra: {len(tetrahedra)}; sampled triangles: {len(sample)}; centroid outside: {len(outside)} "
          f"({len(outside) / len(sample) * 100:.2f} %)")
    if outside:
        unique, counts = np.unique(fracture_ids[outside], return_counts=True)
        order = np.argsort(counts)[::-1][:10]
        print("fractures with most of them: " + ", ".join(f"{unique[i]} ({counts[i]})" for i in order))
        example = points[triangles[outside[0]]].mean(axis=0)
        print(f"example centroid: {example}")


main(sys.argv[1], sys.argv[2], int(sys.argv[3]) if len(sys.argv) > 3 else 20000)
