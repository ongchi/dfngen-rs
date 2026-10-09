# dfngen-rs

Rust implementation of discrete fracture networks (DFN) generator based on DFNGen from [dfnWorks](https://github.com/lanl/dfnWorks).

## Usage

```bash
cargo run --release -- <input_file.toml> <output_folder>
```

## Input file

The input is one TOML file. `tests/data` contains examples. All keys are
necessary, except those marked as optional below. An unknown key is an error.

```toml
seed = 12345                       # 0 = use the system time
h = 0.1                            # minimum feature size (FRAM)
angle_unit = "radian"              # or "degree"; applies to all angles in the file

[domain]
size = [10.0, 10.0, 10.0]          # x, y, z; the domain center is the origin
size_increase = [1.0, 1.0, 1.0]    # larger volume for the fracture centers
layers = [[-5.0, 0.0], [0.0, 5.0]]             # optional; -z, +z of each layer
regions = [[-3.0, 3.0, -3.0, 3.0, -3.0, 3.0]]  # optional; -x, +x, -y, +y, -z, +z
polygon_boundary = [[-4.0, -4.0], [4.0, -4.0], [0.0, 4.0]]  # optional; x, y of each vertex

[generation]
stop_condition = "n_poly"          # or "p32"
n_poly = 60                        # necessary only for "n_poly"
radii_list_increase = 0.1
rejects_per_fracture = 10
force_large_fractures = false
remove_fractures_less_than = 0.0
insert_user_rectangles_first = true

[fram]
disable = false
relaxed = false
triple_intersections = true

[clusters]
boundary_faces = [true, true, false, false, false, false]  # +x, -x, +y, -y, +z, -z
ignore_boundary_faces = true
keep_only_largest = false
keep_isolated_fractures = false

[output]
print_reject_reasons = false
rejected_fractures = false

# One table for each stochastic fracture family (optional)
[[family]]
shape = "ellipse"                  # or "rectangle"
number_of_points = 8               # ellipse only
probability = 0.4
aspect_ratio = 1.0
layer = 0                          # 0 = whole domain, 1 = first layer
region = 0                         # 0 = whole domain, 1 = first region
beta = 0.5                         # optional; without it the rotation is uniform
p32_target = 0.5                   # necessary only for stop_condition "p32"
radius = { distribution = "truncated_power_law", alpha = 2.5, min = 1.0, max = 5.0 }
orientation = { theta = 0.0, phi = 0.0, kappa = 8.0 }

# One table for each user-defined fracture (optional)
[[user_fracture]]
shape = "rectangle"                # or "ellipse", with number_of_points
angle_unit = "degree"              # optional; replaces the angle unit of the file
radius = 4.0
aspect_ratio = 1.0
beta = 0.0
translation = [0.0, 0.0, 0.0]
normal = [0.0, 0.0, 1.0]

[[user_fracture]]
shape = "polygon"                  # or "ellipse", or "rectangle" (4 vertices)
vertices = [[2.5, 2.0, 1.0], [2.5, 0.0, 3.0], [2.5, -2.0, 1.0]]
```

Radius distributions:

| `distribution` | Parameters |
|----------------|------------|
| `log_normal` | `log_mean`, `sd`, `min`, `max` |
| `truncated_power_law` | `alpha`, `min`, `max` |
| `exponential` | `mean`, `min`, `max` |
| `constant` | `value` |

Orientation of a family: `theta` and `phi`, or `trend` and `plunge`, or `dip`
and `strike`, each with `kappa` (Fisher distribution).

Orientation of a user-defined fracture: `normal = [x, y, z]`, or
`trend_plunge = [trend, plunge]`, or `dip_strike = [dip, strike]`.

Order of the fractures: all ellipse families come before all rectangle
families. User-defined polygons are inserted first, then rectangles and
ellipses in the order which `insert_user_rectangles_first` gives. Each group
keeps the order of the file.

### dfnWorks input files

`scripts/convert_dat_to_toml.py <input.dat> <output.toml>` converts a dfnGen
input file of pydfnworks, with its user-defined fracture files and its polygon
boundary file, to this format.

## Output files

All files are in `<output_folder>`. A "final" fracture is an accepted
fracture which was not removed after the generation (isolated fracture, or
not in a kept cluster). Final fractures have the numbers 1, 2, ... in the
order of acceptance. All files are text files.

The log of the run goes to the console (stdout), not to a file. The level is
`info`. Set the environment variable `RUST_LOG` to change it, for example
`RUST_LOG=warn`.

| File | Content |
|---|---|
| `fractures.vtp` | All accepted fractures as polygons with their properties (VTK PolyData, for ParaView or PyVista). See below. |
| `connectivity.dat` | Line N: numbers of the final fractures which intersect final fracture N. |
| `intersection_list.dat` | One line for each intersection: the two fracture numbers, mid point and length. A negative number is a domain face: -1 top, -2 bottom, -3 left, -4 front, -5 right, -6 back. |
| `triple_points.dat` | Coordinates of the triple intersection points. Written only if `triple_intersections = true`. |
| `rejections.dat` | Number of rejected fractures for each rejection reason. |
| `rejectsPerAttempt.dat` | Line N: number of rejected attempts before accepted fracture N. |
| `userFractureRejections.dat` | Number and type of each rejected user-defined fracture. Written only if a user-defined fracture was rejected. |
| `rejected_fractures.vtp` | The rejected fractures as polygons with their properties (VTK PolyData). Written only if `rejected_fractures = true`. See below. |

### Cell data of `fractures.vtp`

The file has one polygon cell for each accepted fracture, with its vertices
in 3D.

| Array | Content |
|---|---|
| `fracture_id` | Number of the final fracture, 0 if the fracture was removed. |
| `accepted_id` | Number in the order of acceptance, starts at 1. |
| `removed` | 1 if the fracture was removed, else 0. |
| `family` | -2 user polygon, -1 user rectangle, 0 user ellipse, > 0 stochastic family. |
| `cluster` | Number of the connected group of fractures. |
| `truncated` | 1 if the domain cut the fracture, else 0. |
| `intersection_count` | Number of intersections with other accepted fractures. |
| `boundary_face_count` | Number of domain faces which the fracture meets (at least two vertices on the face). For a final fracture, this is the number of its lines with a domain face in `intersection_list.dat`. |
| `x_radius`, `y_radius`, `aspect_ratio` | Size before truncation. |
| `area` | Area after truncation. |
| `normal` | Unit normal vector. |
| `translation` | Translation of the fracture center from the origin. |

### Cell data of `rejected_fractures.vtp`

The file has one polygon cell for each rejected fracture, at the position of
its last attempt. A stochastic fracture is rejected when FRAM rejects it and
it has no re-translation left (`rejects_per_fracture`). A user-defined
fracture has only one attempt. The file can be very large.

A fracture has no cell if its last attempt was outside the domain, because it
then has no valid polygon.

| Array | Content |
|---|---|
| `family` | Same numbers as in `fractures.vtp`. |
| `reject_code` | Reason for the rejection of the last attempt. See the table below. |
| `truncated` | 1 if the domain cut the fracture, else 0. |
| `x_radius`, `y_radius`, `aspect_ratio` | Size before truncation. |
| `normal` | Unit normal vector. |
| `translation` | Translation of the fracture center from the origin. |

| `reject_code` | Reason |
|---|---|
| -2 | Intersection shorter than `h`. |
| -5 | Intersection too near to another intersection. |
| -6 | Fracture too near to an edge of another fracture. |
| -7 | Intersecting fractures are in the same plane. |
| -10 | Triple intersection, and `triple_intersections = false`. |
| -11 | Intersection too near to a previous intersection. |
| -12 | Triple intersection with an angle that is too small. |
| -13 | Triple intersection point too near to an end point of an intersection. |
| -14 | Triple intersection point too near to another triple intersection point. |

`fractures.vtp` and `rejected_fractures.vtp` together do not have all
generated fractures. An accepted fracture which is removed later because it
is smaller than `remove_fractures_less_than`, or outside the polygon boundary,
is in neither file.
