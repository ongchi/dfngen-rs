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
min_feature_size = 0.1             # FRAM length limit; minimum radius = 3 times it; eps = 1e-8 times it
domain = "domain.vtu"              # the mesh file which gives the domain (see below)
user_fractures = "fractures.vtp"   # optional; the file with the user-defined fractures (see below)

[generation]
stop_condition = "n_poly"          # or "p32"
n_poly = 60                        # necessary only for "n_poly"
radii_list_increase = 0.1
rejects_per_fracture = 10
force_large_fractures = false
remove_fractures_less_than = 0.0

[fram]
enable = true
relaxed = false
triple_intersections = true
check_user_fractures = true        # optional, default true; false = no FRAM for two user-defined fractures (see below)

[clusters]
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
region = "marble"                  # optional; a region of the mesh (see below). Without it: the whole mesh
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

Family probability: for each new fracture, the generator selects a family at
random, in proportion to its `probability`. The effect is different for each
stop condition:

| | `stop_condition = "n_poly"` | `stop_condition = "p32"` |
|---|---|---|
| Family selection | By `probability` only. | By `probability`, but only among the families which are below their `p32_target`. |
| Final quantity of each family | `probability` controls it. | `p32_target` controls it. `probability` does not. |
| What `probability` changes | The number of fractures of each family. | The order in which the families get their fractures, and thus the network. |

- With `n_poly`, `probability` is the fraction of the insertion attempts, not
  of the accepted fractures. FRAM rejects more fractures of some families, thus
  the final fractions can be different.
- With `p32`, a family which reaches its `p32_target` is removed from the
  selection. Its probability is added in equal parts to the other families.
- The sum of the probabilities must be 1. A sum outside 0.999 to 1.001 gives a
  warning, not an error, and the selection is then not in proportion to the
  values.

Each fracture is whole: the domain does not cut it. Thus:

- A fracture can have vertices outside the domain. A fracture with no part in
  the domain is rejected.
- FRAM also checks the intersections outside the domain, and they are written.
- The area and the P32 values include the part outside the domain. Thus a `p32`
  run stops with a lower intensity in the domain than `p32_target`.

### Domain mesh

`domain = "domain.vtu"` gives the domain as a VTK UnstructuredGrid file.
A relative path starts at the folder of the input file. The key must be
before the first table (`[generation]`) of the input file. The mesh is the only
definition of the domain. `tests/data/mesh_domain` is an example, and
`tests/data/stochastic_npoly` has a box domain (a mesh of one hexahedron).

- The file must be an ASCII `.vtu` file. In ParaView, select the data mode
  "Ascii" in the save dialog. In PyVista, use `mesh.save("domain.vtu", binary=False)`.
- The cells can be tetrahedra, voxels, hexahedra, wedges and pyramids (VTK
  types 10 to 14).
- The bounding box of the mesh is the box domain: a fracture with no part
  in this box is rejected. The mesh can be at any position.
- The center of each stochastic fracture is in a cell of the mesh. The centers
  are sampled in the cells, uniform by volume. A fracture with its center in
  the mesh can have a part outside the mesh; the fracture is not cut.
- A user-defined fracture is rejected only if it intersects no cell (or no
  cell of its region, see "User-defined fractures"). Its center can be outside
  the mesh.
- The volume for P32 is the sum of the cell volumes.
- All coordinates in the input file and in the output files are the coordinates
  of the mesh file.

Regions of the mesh: each cell data array with a name `region__<name>`
(optional) is one region with the name `<name>`. The array has one value for
each cell (one component). A cell with a value above 0 is in the region. Thus
each region is one scalar array, which ParaView can show and use in a
Threshold filter.

`region = "<name>"` of a family puts the fracture centers of the family only in
the cells of that region, and the P32 volume of the family is the volume of
these cells. A family without the `region` key uses the whole mesh. A name
which the mesh does not have, and a region with no cell, are errors.

A cell can be in more than one region, and in none. The volume of a region is
the sum of the volumes of all its cells, thus the sum of the region volumes can
be larger than the volume of the mesh. In this example, the first cell is in
the regions `marble` and `north`, the second cell is in the region `north`, and
the third cell has no region (`tests/data/all_input_options` is an example):

```xml
<CellData>
<DataArray type="Int32" Name="region__marble" format="ascii">1 0 0</DataArray>
<DataArray type="Int32" Name="region__north" format="ascii">1 1 0</DataArray>
</CellData>
```

A cell data array with the name `region` (the format of the versions before
0.7.0) is ignored, as all other cell data.

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
`radius / sqrt(2)` and `radius * aspect_ratio / sqrt(2)`.

Order of the fractures: the user-defined fractures are inserted first, in the
order of their file. The number of a family is its position in the input file.

### User-defined fractures

`user_fractures = "fractures.vtp"` gives the user-defined fractures as a VTK
PolyData file. A relative path starts at the folder of the input file. The key
must be before the first table (`[generation]`) of the input file.
`tests/data/user_rectangles` is an example.

- The file must be an ASCII `.vtp` file. In ParaView, select the data mode
  "Ascii" in the save dialog. In PyVista, use `mesh.save("fractures.vtp", binary=False)`.
- Each polygon cell is one fracture. Its points are the vertices of the
  fracture, in the order of the cell. Give points in one plane and a convex
  polygon; the code does not check this.
- A polygon must have 3 points or more. The fractures are inserted in the
  order of the cells.
- The coordinates are those of the mesh file.
- The cell data arrays `region__<name>` (optional) give the regions of each
  fracture, as in the mesh file: a polygon with a value above 0 in the array
  `region__<name>` has the region `<name>`. A fracture can have many regions,
  and none. A fracture with regions is accepted only if it intersects a cell
  of the mesh in one or more of its regions. If not, it is rejected as a
  fracture outside the domain. The fracture is not cut at a region: a part of
  it can be in other regions. A fracture with no region is checked against all
  cells. A region of a fracture which the mesh does not have, or which has no
  cell, is an error. `tests/data/mesh_domain` is an example.
- The cell data array `fracture_id` (optional) gives the fracture of each
  polygon: one integer of 0 or more for each cell of the file. Polygons with
  the same value are the pieces of one fracture, for example a surface from
  field data. Without the array, each polygon is one fracture. For two pieces
  of one fracture:
  - FRAM does not check them against each other, and their contact is not an
    intersection in the output.
  - They are always in the same cluster, also if they do not touch. Thus a
    fracture of 2 or more accepted pieces is not removed as an isolated
    fracture.
  - A fracture which intersects the two pieces can cross their shared edge.
    FRAM does not compare its intersections with the two pieces to each other.

  - The domain check is for the fracture. If one piece has a part in the
    domain (it intersects a cell of the mesh, or of its regions), no piece of
    the fracture is rejected as outside the domain. Thus an accepted piece
    can be fully outside the mesh or in cells of another region, and its
    area is in the total area and in P32.

  The pieces must not cross each other; the code does not check this. FRAM
  accepts or rejects each piece alone: with `check_user_fractures = true`, it
  can reject one piece of a fracture and accept another. Each piece has its
  own line in the output files.
- `fram.check_user_fractures = false` turns off FRAM for an intersection of
  two user-defined fractures. Then each user-defined fracture in the domain
  (and in its region) is accepted. The intersection is recorded as FRAM finds
  it: it is not shortened, and such a pair gives no triple intersection point.
  An intersection of a stochastic fracture with a user-defined fracture keeps
  all FRAM checks. The network then does not obey the minimum feature size `min_feature_size`
  at the user-defined fractures.
- The point data and the other cell data of the file are ignored. Thus
  `fractures.vtp` of a run can be the user-defined fractures of a new run. It
  has no `region__<name>` array, so its fractures have no region. Its `fracture_id`
  array has a different number for each polygon, so each polygon is one
  fracture.
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
| `intersection_list.dat` | One line for each intersection: the two fracture numbers, mid point and length. |
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
| `user_fracture_id` | Value of the `fracture_id` array of the user-defined fracture file. -1 for a stochastic fracture and for a file without the array. |
| `cluster` | Number of the connected group of fractures. |
| `intersection_count` | Number of intersections with other accepted fractures. |
| `x_radius`, `y_radius`, `aspect_ratio` | Size. |
| `area` | Area. |
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
| `user_fracture_id` | Same values as in `fractures.vtp`. |
| `reject_code` | Reason for the rejection of the last attempt. See the table below. |
| `x_radius`, `y_radius`, `aspect_ratio` | Size. |
| `normal` | Unit normal vector. |
| `translation` | Translation of the fracture center from the origin. |

| `reject_code` | Reason |
|---|---|
| -2 | Intersection shorter than `min_feature_size`. |
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
