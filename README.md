# dfngen-rs

Rust implementation of discrete fracture networks (DFN) generator based on DFNGen from [dfnWorks](https://github.com/lanl/dfnWorks).

## Usage

```bash
cargo run --release -- <input_file.toml> <output_folder>
```

## Input file

The input is one TOML file, one mesh file for the domain, and an optional file
with the user-defined fractures. Each subfolder of
`tests/data` is an example. All keys are necessary, except those marked as
optional below. An unknown key is an error.

```toml
seed = 12345                       # 0 = use the system time
h = 0.1                            # minimum feature size (FRAM)
user_fractures = "fractures.vtp"   # optional; the file with the user-defined fractures (see below)

[domain]
mesh = "domain.vtu"                # the mesh file which gives the domain (see below)
truncate_fractures = true          # false = keep each fracture whole (see below)

[generation]
stop_condition = "n_poly"          # or "p32"
n_poly = 60                        # necessary only for "n_poly"
radii_list_increase = 0.1
rejects_per_fracture = 10
force_large_fractures = false
remove_fractures_less_than = 0.0

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
number_of_points = 8               # 3 or more; 4 gives a rectangle
probability = 0.4
aspect_ratio = 1.0
region = 0                         # 0 = whole domain, N = the cells with the region N
beta = 30.0                        # optional; without it the rotation is uniform
p32_target = 0.5                   # necessary only for stop_condition "p32"
radius = { distribution = "truncated_power_law", alpha = 2.5, min = 1.0, max = 5.0 }
orientation = { theta = 0.0, phi = 0.0, kappa = 8.0 }
```

Radius distributions:

| `distribution` | Parameters |
|----------------|------------|
| `log_normal` | `log_mean`, `sd`, `min`, `max`, `sampling` (optional) |
| `truncated_power_law` | `alpha`, `min`, `max` |
| `exponential` | `mean`, `min`, `max` |
| `constant` | `value` |

`sampling` of `log_normal` is `"rejection"` (default) or `"inverse_cdf"`.
`"rejection"` samples the full distribution until a value is in `[min, max]`. It
is fast when `[min, max]` contains most of the distribution, and it stops with an
error after 1000 rejected values. `"inverse_cdf"` uses one random number for each
radius, also when `[min, max]` is in a tail of the distribution. The two methods
give different networks for the same seed.

`domain.truncate_fractures = false` keeps each fracture whole. Then:

- A fracture can have vertices outside the domain. A fracture with no part in
  the domain is rejected, as with truncation.
- A fracture which crosses a domain face is connected to that face for
  `clusters.boundary_faces`.
- FRAM also checks the intersections outside the domain, and they are written.
- The area and the P32 values include the part outside the domain. Thus a `p32`
  run stops with a lower intensity in the domain than `p32_target`.
- `truncated` and `boundary_face_count` are 0, and `intersection_list.dat` has no
  line between a fracture and a domain face.

### Domain mesh

`domain.mesh = "domain.vtu"` gives the domain as a VTK UnstructuredGrid file.
A relative path starts at the folder of the input file. The mesh is the only
definition of the domain. `tests/data/mesh_domain` is an example, and
`tests/data/stochastic_npoly` has a box domain (a mesh of one hexahedron).

- The file must be an ASCII `.vtu` file. In ParaView, select the data mode
  "Ascii" in the save dialog. In PyVista, use `mesh.save("domain.vtu", binary=False)`.
- The cells can be tetrahedra, voxels, hexahedra, wedges and pyramids (VTK
  types 10 to 14).
- The bounding box of the mesh is the box domain: `truncate_fractures` and
  `clusters.boundary_faces` use this box. The mesh can be at any position.
- The center of each stochastic fracture is in a cell of the mesh. The centers
  are sampled in the cells, uniform by volume. A fracture with its center in
  the mesh can have a part outside the mesh; the fracture is cut only at the
  bounding box.
- A user-defined fracture is rejected only if it intersects no cell. Its center
  can be outside the mesh.
- The volume for P32 is the sum of the cell volumes.
- All coordinates in the input file and in the output files are the coordinates
  of the mesh file.

Regions of the mesh: the cell data array `region` (optional) gives one integer
of 0 or more for each cell. `region = N` of a family (N > 0) puts the fracture
centers of the family only in the cells with the value N, and the P32 volume of
the family is the volume of these cells. `region = 0` of a family is the whole
mesh. Without the `region` array, all cells have the region 0. A family with a
region value which no cell has is an error.

A cell is split into tetrahedra for the sampling of the centers, the check of
the user-defined fractures and the volume. This is approximate for a
hexahedron with a face which is not flat.

Orientation of a family: `theta` and `phi`, or `trend` and `plunge`, or `dip`
and `strike`, each with `kappa` (Fisher distribution).

All angles (`theta`, `phi`, `trend`, `plunge`, `dip`, `strike`, `beta`) are in
degrees.

Shape of a stochastic fracture: a family makes
a regular polygon with `number_of_points` vertices on a circle of the given
radius. The polygon is then scaled in one direction by `aspect_ratio`.
Many points give an ellipse. 4 points give a rectangle: its half sides are
`radius / sqrt(2)` and `radius * aspect_ratio / sqrt(2)`. Because of rounding,
a corner of such a rectangle can be a little outside a domain face which it
must only touch. Then the domain cuts the rectangle.

Order of the fractures: the user-defined fractures are inserted first, in the
order of their file. The number of a family is its position in the input file.

### User-defined fractures

`user_fractures = "fractures.vtp"` gives the user-defined fractures as a VTK
PolyData file. A relative path starts at the folder of the input file. The key
must be before the first table (`[domain]`) of the input file.
`tests/data/user_rectangles` is an example.

- The file must be an ASCII `.vtp` file. In ParaView, select the data mode
  "Ascii" in the save dialog. In PyVista, use `mesh.save("fractures.vtp", binary=False)`.
- Each polygon cell is one fracture. Its points are the vertices of the
  fracture, in the order of the cell. Give points in one plane and a convex
  polygon; the code does not check this.
- A polygon must have 3 points or more. The fractures are inserted in the
  order of the cells.
- The coordinates are those of the mesh file.
- The point data and the cell data of the file are ignored. Thus
  `fractures.vtp` of a run can be the user-defined fractures of a new run.
- Cells which are not polygons (vertices, lines, strips) are ignored with a
  warning.
- A user-defined fracture has no radius in its file. `x_radius`, `y_radius` and
  `aspect_ratio` in the output are estimates from the vertices: `x_radius` is
  half the length of the first edge, and `y_radius` is half the distance
  between the vertices at 1/4 and 3/4 of the polygon.
  `generation.remove_fractures_less_than` uses this `x_radius`.

### dfnWorks input files

`scripts/convert_dat_to_toml.py <input.dat> <output.toml>` converts a dfnGen
input file of pydfnworks, with its user-defined fracture files, to this format.
It writes the box domain of dfnGen as a mesh of one hexahedron to
`<output>.vtu`, and the user-defined fractures as polygons to
`<output>_user_fractures.vtp`. For an ellipse or a rectangle which dfnGen
defines by parameters, it calculates the vertices as for a stochastic fracture. It stops with an error if the dfnGen file has layers, regions
or a polygon boundary. A rectangle of dfnGen becomes a polygon with 4
points, and its radius values are multiplied by sqrt(2) (in dfnGen, the radius
of a rectangle is half the side). Angles in radians are changed to degrees. The ellipses do not have the same vertices
as in dfnGen, and the fracture centers are always in the domain (dfnGen can
put them in a larger box, `domainSizeIncrease`). Thus the network is not the
same as that of dfnGen for the same seed.

## Output files

All files are in `<output_folder>`. A "final" fracture is an accepted
fracture which was not removed after the generation (isolated fracture, or
not in a kept cluster). Final fractures have the numbers 1, 2, ... in the
order of acceptance. All files are text files.

The log of the run goes to the console (stdout), not to a file. The level is
`info`. Set the environment variable `RUST_LOG` to change it, for example
`RUST_LOG=warn`.

During the insertion of the stochastic fractures, a progress bar is shown on
stderr, if stderr is a terminal:

```
[00:00:12] [################>-----------------------]  42 %  25/60 fractures  ETA 00:00:17
```

With `stop_condition = "p32"`, the percent is the part of the total target P32
which the families have, and the bar shows the number of accepted fractures
without a total. The large fractures are inserted first, thus the remaining time
of such a run is approximate.

| File | Content |
|---|---|
| `fractures.vtp` | The final fractures as polygons with their properties (VTK PolyData, for ParaView or PyVista). See below. |
| `fractures_removed.vtp` | The accepted fractures which are not final fractures, as polygons with their properties (VTK PolyData). Written only if a fracture was removed. See below. |
| `connectivity.dat` | Line N: numbers of the final fractures which intersect final fracture N. |
| `intersection_list.dat` | One line for each intersection: the two fracture numbers, mid point and length. A negative number is a domain face: -1 top, -2 bottom, -3 left, -4 front, -5 right, -6 back. |
| `triple_points.dat` | Coordinates of the triple intersection points. Written only if `triple_intersections = true`. |
| `rejections.dat` | Number of rejected fractures for each rejection reason. |
| `rejectsPerAttempt.dat` | Line N: number of rejected attempts before accepted fracture N. |
| `userFractureRejections.dat` | Number of each rejected user-defined fracture (the position of its polygon in the `.vtp` file, from 1). Written only if a user-defined fracture was rejected. |
| `fractures_rejected.vtp` | The rejected fractures as polygons with their properties (VTK PolyData). Written only if `rejected_fractures = true`. See below. |

### Cell data of `fractures.vtp`

The file has one polygon cell for each final fracture, with its vertices
in 3D. Cell N is final fracture N.

| Array | Content |
|---|---|
| `fracture_id` | Number of the final fracture, starts at 1. |
| `accepted_id` | Number in the order of acceptance, starts at 1. |
| `family` | 0 user-defined fracture, > 0 stochastic family. |
| `cluster` | Number of the connected group of fractures. |
| `truncated` | 1 if the domain cut the fracture, else 0. |
| `intersection_count` | Number of intersections with other accepted fractures. |
| `boundary_face_count` | Number of domain faces which the fracture meets (at least two vertices on the face). For a final fracture, this is the number of its lines with a domain face in `intersection_list.dat`. |
| `x_radius`, `y_radius`, `aspect_ratio` | Size before truncation. |
| `area` | Area after truncation. |
| `normal` | Unit normal vector. |
| `translation` | Translation of the fracture center from the origin. |

### Cell data of `fractures_removed.vtp`

The file has one polygon cell for each accepted fracture which is not a final
fracture (isolated fracture, or not in a kept cluster), in the order of
acceptance. The arrays are those of `fractures.vtp` without `fracture_id`.

### Cell data of `fractures_rejected.vtp`

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

`fractures.vtp`, `fractures_removed.vtp` and `fractures_rejected.vtp` together
do not have all generated fractures. An accepted fracture which is removed later because it
is smaller than `remove_fractures_less_than` is in none of these files.
