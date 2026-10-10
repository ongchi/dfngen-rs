#!/usr/bin/env python3
"""Convert a dfnGen input file of pydfnworks (`.dat`) to the TOML input of dfngen-rs.

Usage: convert_dat_to_toml.py <input.dat> <output.toml>

The user-defined fracture files to which the `.dat` file refers are read too.
Their content goes into the TOML file. A relative path in the `.dat` file is
relative to the folder of the `.dat` file.

The domain of dfngen-rs is a mesh. The script writes the box of `domainSize`
as a mesh of one hexahedron to `<output>.vtu`, next to the TOML file. Layers,
regions and a polygon boundary have no equivalent, and the script stops with
an error if the `.dat` file uses them. `domainSizeIncrease` is dropped: the
fracture centers of dfngen-rs are always in the domain.

Keys which dfngen-rs does not use are dropped (`domainCenter`, `framOn`,
`ecpmOutput`, `nFracFam`, `minimum_fracture_size`).

A rectangle becomes a polygon with 4 points. Its radius in dfnGen is half the
side; in dfngen-rs it is the distance from the center to a corner of the
square with that side. Thus the radius values of a rectangle are multiplied
by sqrt(2).

All angles in the TOML file are in degrees. Angles in radians (`angleOption: 0`,
or `AngleOption: radian` in a user-defined fracture file) are changed to degrees.

Needs only the Python standard library.
"""

from __future__ import annotations

import math
import re
import sys
from pathlib import Path

BRACKET_GROUP = re.compile(r"[{\[(]([^}\])]*)[}\])]")
SQRT_2 = math.sqrt(2.0)
RECTANGLE_POINT_COUNT = 4


class DatFile:
    """A file of `label: value` entries. A value is a word or a bracketed list."""

    def __init__(self, path: Path):
        self.path = path
        self.text = path.read_text()

    def _after(self, label: str) -> str:
        match = re.search(rf"(?:^|\s){re.escape(label)}:(?=\s|$)", self.text)
        if match is None:
            sys.exit(f"{self.path}: variable not found: {label}")
        return self.text[match.end() :]

    def word(self, label: str) -> str:
        words = self._after(label).split(maxsplit=1)
        if not words:
            sys.exit(f"{self.path}: no value for {label}")
        return words[0]

    def words(self, label: str, count: int) -> list[str]:
        words = self._after(label).split()[:count]
        if len(words) != count:
            sys.exit(f"{self.path}: expected {count} values for {label}")
        return words

    def flag(self, label: str) -> bool:
        return int(self.word(label)) != 0

    def list(self, label: str) -> list[str]:
        """The values of the bracketed list which follows the label."""
        return self.groups(label, 1)[0]

    def groups(self, label: str, count: int) -> list[list[str]]:
        """The values of the `count` bracketed lists which follow the label."""
        return read_groups(self.path, self._after(label), count, label)

    def rest(self, label: str) -> str:
        return self._after(label)


def read_groups(path: Path, text: str, count: int, label: str) -> list[list[str]]:
    groups = []
    position = 0
    for _ in range(count):
        match = BRACKET_GROUP.match(text, skip_whitespace(text, position))
        if match is None:
            sys.exit(f"{path}: expected {count} bracketed lists for {label}")
        groups.append([value.strip() for value in match.group(1).split(",") if value.strip()])
        position = match.end()
    return groups


def skip_whitespace(text: str, position: int) -> int:
    while position < len(text) and text[position].isspace():
        position += 1
    return position


def number(value: str) -> str:
    """A TOML float with the same binary value as the text of the `.dat` file."""
    return repr(float(value))


def rectangle_radius(value: str) -> str:
    """The radius of the 4-point polygon for a rectangle with this half side."""
    return repr(float(value) * SQRT_2)


def numbers(values: list[str]) -> str:
    return "[" + ", ".join(number(value) for value in values) + "]"


def boolean(is_set: bool) -> str:
    return "true" if is_set else "false"


def angle(value: str, is_degree: bool) -> str:
    """A TOML float with the angle in degrees."""
    return number(value) if is_degree else repr(math.degrees(float(value)))


def angles(values: list[str], is_degree: bool) -> str:
    return "[" + ", ".join(angle(value, is_degree) for value in values) + "]"


def stop_if_not_convertible(dat: DatFile) -> None:
    """Stops if the file uses a domain option which a box mesh cannot give."""
    for label, description in (
        ("numOfLayers", "layers"),
        ("numOfRegions", "regions"),
        ("polygonBoundaryFlag", "a polygon boundary"),
    ):
        if int(dat.word(label)) != 0:
            sys.exit(f"{dat.path}: {description} cannot be converted ({label} is not 0)")


def box_mesh(size: list[float]) -> str:
    """An ASCII `.vtu` file with one hexahedron: a box with its center at the origin."""
    half_x, half_y, half_z = (length / 2.0 for length in size)
    points = [
        f"{x!r} {y!r} {z!r}"
        for z in (-half_z, half_z)
        for x, y in ((-half_x, -half_y), (half_x, -half_y), (half_x, half_y), (-half_x, half_y))
    ]
    return "\n".join(
        [
            '<?xml version="1.0"?>',
            '<VTKFile type="UnstructuredGrid" version="1.0" byte_order="LittleEndian">',
            "<UnstructuredGrid>",
            '<Piece NumberOfPoints="8" NumberOfCells="1">',
            "<Points>",
            '<DataArray type="Float64" Name="Points" NumberOfComponents="3" format="ascii">',
            *points,
            "</DataArray>",
            "</Points>",
            "<Cells>",
            '<DataArray type="Int64" Name="connectivity" format="ascii">',
            "0 1 2 3 4 5 6 7",
            "</DataArray>",
            '<DataArray type="Int64" Name="offsets" format="ascii">',
            "8",
            "</DataArray>",
            '<DataArray type="UInt8" Name="types" format="ascii">',
            "12",
            "</DataArray>",
            "</Cells>",
            "</Piece>",
            "</UnstructuredGrid>",
            "</VTKFile>",
            "",
        ]
    )


def convert_general_options(dat: DatFile, mesh_file_name: str) -> list[str]:
    stop_condition = int(dat.word("stopCondition"))
    if stop_condition not in (0, 1):
        sys.exit(f"{dat.path}: stopCondition must be 0 or 1")

    lines = [
        f"seed = {int(dat.word('seed'))}",
        f"h = {number(dat.word('h'))}",
        "",
        "[domain]",
        f'mesh = "{mesh_file_name}"',
        "truncate_fractures = true",
    ]

    lines += ["", "[generation]"]
    if stop_condition == 0:
        lines += ['stop_condition = "n_poly"', f"n_poly = {int(dat.word('nPoly'))}"]
    else:
        lines.append('stop_condition = "p32"')
    lines += [
        f"radii_list_increase = {number(dat.word('radiiListIncrease'))}",
        f"rejects_per_fracture = {int(dat.word('rejectsPerFracture'))}",
        f"force_large_fractures = {boolean(dat.flag('forceLargeFractures'))}",
        f"remove_fractures_less_than = {number(dat.word('removeFracturesLessThan'))}",
        "",
        "[fram]",
        f"disable = {boolean(dat.flag('disableFram'))}",
        f"relaxed = {boolean(dat.flag('rFram'))}",
        f"triple_intersections = {boolean(dat.flag('tripleIntersections'))}",        "",
        "[clusters]",
        "boundary_faces = ["
        + ", ".join(boolean(int(face) != 0) for face in dat.list("boundaryFaces"))
        + "]",
        f"ignore_boundary_faces = {boolean(dat.flag('ignoreBoundaryFaces'))}",
        f"keep_only_largest = {boolean(dat.flag('keepOnlyLargestCluster'))}",
        f"keep_isolated_fractures = {boolean(dat.flag('keepIsolatedFractures'))}",
        "",
        "[output]",
        f"print_reject_reasons = {boolean(dat.flag('printRejectReasons'))}",
        f"rejected_fractures = {boolean(dat.flag('outputAllRadii'))}",
    ]
    return lines


def resolve(dat: DatFile, file_name: str) -> Path:
    path = Path(file_name)
    return path if path.is_absolute() else dat.path.parent / path


def convert_families(dat: DatFile, prefix: str, probabilities: list[str]) -> list[str]:
    """The `[[family]]` tables of the ellipse (`e`) or rectangle (`r`) families.

    Takes the probabilities of these families from the front of `probabilities`.
    """
    family_count = int(dat.word("nFamEll" if prefix == "e" else "nFamRect"))
    if family_count == 0:
        return []

    def values(label: str) -> list[str]:
        return dat.list(prefix + label)

    is_p32 = int(dat.word("stopCondition")) == 1
    is_degree = dat.flag("angleOption")
    angle_names = [("theta", "phi"), ("trend", "plunge"), ("dip", "strike")][
        int(dat.word("orientationOption"))
    ]

    has_beta = [int(value) != 0 for value in values("betaDistribution")]
    # These lists have a value only for the families which use them
    betas = iter(values("beta"))
    distribution_values = {
        label: iter(values(label))
        for label in (
            "LogMean", "sd", "LogMin", "LogMax", "alpha", "min", "max",
            "ExpMean", "ExpMin", "ExpMax", "const",
        )
    }  # fmt: skip

    def next_value(label: str) -> str:
        try:
            return number(next(distribution_values[label]))
        except StopIteration:
            sys.exit(f"{dat.path}: not sufficient values for {prefix}{label}")

    def next_length(label: str) -> str:
        """The next value of a parameter which has the unit of the radius."""
        value = next_value(label)
        return rectangle_radius(value) if prefix == "r" else value

    def next_log_mean() -> str:
        value = next_value("LogMean")
        return repr(float(value) + math.log(SQRT_2)) if prefix == "r" else value

    lines = []
    for i in range(family_count):
        lines += ["", "[[family]]"]
        point_count = int(values("numPoints")[i]) if prefix == "e" else RECTANGLE_POINT_COUNT
        lines += [
            f"number_of_points = {point_count}",
            f"probability = {number(probabilities.pop(0))}",
            f"aspect_ratio = {number(values('aspect')[i])}",
            "region = 0",
        ]
        if has_beta[i]:
            lines.append(f"beta = {angle(next(betas), is_degree)}")
        if is_p32:
            lines.append(f"p32_target = {number(values('_p32Targets')[i])}")

        distribution = int(values("distr")[i])
        if distribution == 1:
            radius = (
                f'distribution = "log_normal", log_mean = {next_log_mean()}, '
                f'sd = {next_value("sd")}, min = {next_length("LogMin")}, '
                f'max = {next_length("LogMax")}'
            )
        elif distribution == 2:
            radius = (
                f'distribution = "truncated_power_law", alpha = {next_value("alpha")}, '
                f'min = {next_length("min")}, max = {next_length("max")}'
            )
        elif distribution == 3:
            radius = (
                f'distribution = "exponential", mean = {next_length("ExpMean")}, '
                f'min = {next_length("ExpMin")}, max = {next_length("ExpMax")}'
            )
        elif distribution == 4:
            radius = f'distribution = "constant", value = {next_length("const")}'
        else:
            sys.exit(f"{dat.path}: unknown distribution type: {distribution}")
        lines.append(f"radius = {{ {radius} }}")

        lines.append(
            f"orientation = {{ {angle_names[0]} = {angle(values(angle_names[0])[i], is_degree)}, "
            f"{angle_names[1]} = {angle(values(angle_names[1])[i], is_degree)}, "
            f"kappa = {number(values('kappa')[i])} }}"
        )
    return lines


def convert_fractures_by_parameters(path: Path, shape: str) -> list[str]:
    dat = DatFile(path)
    count = int(dat.word("nUserEll" if shape == "ellipse" else "nUserRect"))

    angle_unit = dat.word("AngleOption")
    if angle_unit not in ("degree", "radian"):
        sys.exit(f"{path}: AngleOption must be degree or radian")
    is_degree = angle_unit == "degree"

    radii = dat.words("Radii", count)
    betas = dat.words("Beta", count)
    aspect_ratios = dat.words("Aspect_Ratio", count)
    translations = dat.groups("Translation", count)
    orientation_key, orientation_label = [
        ("normal", "Normal"),
        ("trend_plunge", "Trend_Plunge"),
        ("dip_strike", "Dip_Strike"),
    ][int(dat.word("userOrientationOption"))]
    orientations = dat.groups(orientation_label, count)
    if shape == "ellipse":
        point_counts = dat.words("Number_of_Vertices", count)

    lines = []
    for i in range(count):
        lines += ["", "[[user_fracture]]"]
        if shape == "ellipse":
            point_count, radius = int(point_counts[i]), number(radii[i])
        else:
            point_count, radius = RECTANGLE_POINT_COUNT, rectangle_radius(radii[i])
        lines += [
            f"number_of_points = {point_count}",
            f"radius = {radius}",
            f"aspect_ratio = {number(aspect_ratios[i])}",
            f"beta = {angle(betas[i], is_degree)}",
            f"translation = {numbers(translations[i])}",
        ]
        if orientation_key == "normal":
            lines.append(f"normal = {numbers(orientations[i])}")
        else:
            lines.append(f"{orientation_key} = {angles(orientations[i], is_degree)}")
    return lines


def fracture_by_vertices(vertices: list[list[str]]) -> list[str]:
    lines = ["", "[[user_fracture]]", "vertices = ["]
    lines += [f"    {numbers(vertex)}," for vertex in vertices]
    lines.append("]")
    return lines


def convert_fractures_by_vertices(path: Path, shape: str) -> list[str]:
    dat = DatFile(path)
    if shape == "ellipse":
        count = int(dat.word("nEllipses"))
        point_count = int(dat.word("nNodes"))
    else:
        count = int(dat.word("nRectangles"))
        point_count = RECTANGLE_POINT_COUNT

    vertices = dat.groups("Coordinates", count * point_count)
    lines = []
    for i in range(count):
        lines += fracture_by_vertices(vertices[i * point_count : (i + 1) * point_count])
    return lines


def convert_polygons(path: Path) -> list[str]:
    dat = DatFile(path)
    count = int(dat.word("nPolygons"))
    text = dat.rest("nPolygons").split(maxsplit=1)[1]

    # Each polygon is its number of nodes, then one {x, y, z} group for each node
    lines = []
    for _ in range(count):
        point_count, text = text.split(maxsplit=1)
        vertices = []
        for _ in range(int(point_count)):
            match = BRACKET_GROUP.match(text, skip_whitespace(text, 0))
            if match is None:
                sys.exit(f"{path}: expected {point_count} nodes for a polygon")
            vertices.append([value.strip() for value in match.group(1).split(",")])
            text = text[match.end() :]
        lines += fracture_by_vertices(vertices)
    return lines


def convert_user_fractures(dat: DatFile) -> list[str]:
    def path_if_enabled(flag_label: str, path_label: str) -> Path | None:
        return resolve(dat, dat.word(path_label)) if dat.flag(flag_label) else None

    # dfngen-rs inserts the fractures in the order of the file. dfnGen inserts
    # the polygons first, then the rectangles and the ellipses in the order
    # which insertUserRectanglesFirst gives, and in each of these groups the
    # fractures by parameters before those by coordinates.
    polygons = []
    if path := path_if_enabled("userPolygonByCoord", "PolygonByCoord_Input_File_Path"):
        polygons += convert_polygons(path)

    rectangles = []
    if path := path_if_enabled("userRectanglesOnOff", "UserRect_Input_File_Path"):
        rectangles += convert_fractures_by_parameters(path, "rectangle")
    if path := path_if_enabled("userRecByCoord", "RectByCoord_Input_File_Path"):
        rectangles += convert_fractures_by_vertices(path, "rectangle")

    ellipses = []
    if path := path_if_enabled("userEllipsesOnOff", "UserEll_Input_File_Path"):
        ellipses += convert_fractures_by_parameters(path, "ellipse")
    if path := path_if_enabled("userEllByCoord", "EllByCoord_Input_File_Path"):
        ellipses += convert_fractures_by_vertices(path, "ellipse")

    if dat.flag("insertUserRectanglesFirst"):
        return polygons + rectangles + ellipses
    return polygons + ellipses + rectangles


def convert(dat: DatFile, mesh_file_name: str) -> str:
    stop_if_not_convertible(dat)

    lines = convert_general_options(dat, mesh_file_name)

    family_count = int(dat.word("nFamEll")) + int(dat.word("nFamRect"))
    probabilities = dat.list("famProb") if family_count > 0 else []
    if len(probabilities) != family_count:
        sys.exit(f"{dat.path}: famProb must have {family_count} values")
    lines += convert_families(dat, "e", probabilities)
    lines += convert_families(dat, "r", probabilities)

    lines += convert_user_fractures(dat)

    return "\n".join(lines) + "\n"


def main() -> int:
    if len(sys.argv) != 3:
        print(__doc__, file=sys.stderr)
        return 2

    dat = DatFile(Path(sys.argv[1]))
    toml_path = Path(sys.argv[2])
    mesh_path = toml_path.with_suffix(".vtu")

    toml_path.write_text(convert(dat, mesh_path.name))
    mesh_path.write_text(box_mesh([float(length) for length in dat.list("domainSize")]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
