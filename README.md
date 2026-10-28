# dfngen-rs

Rust implementation of discrete fracture networks (DFN) generator based on DFNGen from [dfnWorks](https://github.com/lanl/dfnWorks).

## Usage

```bash
cargo run --release -- <input_file.toml> <output_folder>
```

## Input file

The input is one TOML file, one mesh file for the domain, and an optional file
with the user-defined fractures. Each subfolder of
`tests/data` is an example of one feature; the comment at the top of its input
file tells which. All keys are necessary, except those marked as
optional below. An unknown key is an error.

```toml
seed = 12345                       # 0 = use the system time
min_feature_size = 0.1             # > 0; FRAM length limit; minimum radius = 3 times it; eps = 1e-8 times it
domain = "domain.vtu"              # the mesh file which gives the domain (see below)
user_fractures = "fractures.vtp"   # optional; the file with the user-defined fractures (see below)

[generation]
n_poly = 60                        # optional; stop at 60 stochastic fractures (user-defined fractures are not included). Without it: each family stops at its p32_target
force_large_fractures = false      # true = each family first gets one fracture with its maximum radius
remove_fractures_less_than = 0.0   # after the generation, remove the accepted fractures with x_radius below it; 0 = none

# Optional section; each key is optional. Without it: FRAM is off.
[fram]
enable = true                      # default false; false = no FRAM check, fractures can intersect with no limit, and the other keys have no effect
relaxed = false                    # default false; true = FRAM rejects only short intersections and fractures in the same plane
triple_intersections = true        # default false; false = reject a fracture which makes a triple intersection; no effect with relaxed = true
check_user_fractures = true        # default false; false = no FRAM for two user-defined fractures (see below)
write_rejected_fractures = false   # default false; true = also write the fractures which FRAM rejected to fractures_rejected.vtp (see "Output files")
rejects_per_fracture = 10          # default 10; FRAM rejects of a stochastic fracture before it is rejected; it moves after each reject

# Necessary section and key. fehm-rs needs the lines: they give the boundaries of the mesh.
[region_boundaries]
enable = true                      # true = write the lines where the boundary of the region of each fracture crosses it (see "Boundary lines of the regions"); false = no lines
# Necessary section; each key is necessary. With enable = false, the other keys have no effect.
[mesh]
enable = true                      # true = write the triangle mesh of the final fractures for fehm-rs (see "Triangle mesh for fehm-rs"); needs region_boundaries.enable = true; false = no mesh
max_size_factor = 10.0             # >= 1; largest triangle size = max_size_factor × min_feature_size
slope = 0.5                        # > 0; growth of the triangle size with the distance from an intersection or a boundary line

[clusters]
keep_only_largest = false          # true = keep only the largest cluster; no effect with keep_isolated_fractures = true
keep_isolated_fractures = false    # true = keep all accepted fractures; false = remove the fractures with no intersection
keep_only_connected_to_surfaces = false # true = keep only the clusters which touch each boundary surface of the domain mesh (see below); no effect with keep_isolated_fractures = true

# One table for each stochastic fracture family (optional)
[[family]]
number_of_points = 8               # 3 or more; 4 gives a rectangle
probability = 0.4                  # chance to select this family; the sum of all families must be 1 (see below)
aspect_ratio = 1.0                 # > 0; y radius / x radius
region = "marble"                  # optional; a region of the mesh (see below). Without it: the whole mesh
beta = 30.0                        # optional; without it the rotation is uniform
p32_target = 0.5                   # necessary without n_poly; not allowed with n_poly
aperture = { distribution = "constant", value = 1e-4 } # necessary; only the mesh uses it; see "Apertures"
radius = { distribution = "truncated_power_law", alpha = 2.5, min = 1.0, max = 5.0 } # see "Radius distributions" below
orientation = { theta = 0.0, phi = 0.0, kappa = 8.0 } # see "Orientation of a family" below
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

Radii lists: before the insertion, each family gets a list of radii, sampled
from its distribution. The list is sorted from large to small, and the
fractures of the family use the radii in this order. Thus the large fractures
go in first. When the list is empty, each new fracture samples its radius from
the distribution.

- With `n_poly`, the list has `ceil(probability × n_poly)` radii. With
  `force_large_fractures = true`, it has `max` and
  `ceil(probability × (n_poly − number of families))` radii.
- Without `n_poly`, the family samples radii until their area gives its
  `p32_target`: `sum(area) × 2 / volume ≥ p32_target`, with the volume of the
  family's region. The area of a fracture with `n` points and radius `r` is
  `(n / 2) × sin(2π / n) × aspect_ratio × r²`. With
  `force_large_fractures = true`, the list starts with `max`. (dfnWorks makes
  this list with a dry run of the generation, with no FRAM check.)
- A fracture which FRAM rejects also uses a radius of the list. Thus with FRAM
  rejections, the list ends before the generation stops. (dfnWorks adds extra
  radii with `radiiListIncrease`; dfngen-rs has no such key.)

Family probability: for each new fracture, the generator selects a family at
random, in proportion to its `probability`. The effect is different for each
stop condition:

| | With `n_poly` | Without `n_poly` (P32) |
|---|---|---|
| Family selection | By `probability` only. | By `probability`, but only among the families which are below their `p32_target`. |
| Final quantity of each family | `probability` controls it. | `p32_target` controls it. `probability` does not. |
| What `probability` changes | The number of fractures of each family. | The order in which the families get their fractures, and thus the network. |

- With `n_poly`, `probability` is the fraction of the insertion attempts, not
  of the accepted fractures. FRAM rejects more fractures of some families, thus
  the final fractions can be different.
- Without `n_poly`, a family which reaches its `p32_target` is removed from the
  selection. Its probability is added in equal parts to the other families.
- The sum of the probabilities must be 1. A sum outside 0.999 to 1.001 gives a
  warning, not an error, and the selection is then not in proportion to the
  values.

Each fracture is whole: the domain does not cut it. Thus:

- A fracture can have vertices outside the domain. A fracture with no part in
  the domain is rejected.
- FRAM also checks the intersections outside the domain, and they are written.
- The area and the P32 values include the part outside the domain. Thus a run
  without `n_poly` stops with a lower intensity in the domain than `p32_target`.

### Boundary lines of the regions

`region_boundaries.enable = true` finds, for each fracture, the lines where
the boundary of its region crosses it, and the lines where the mesh surface
crosses it. The fracture is not cut: it stays whole. `fracture_boundaries.vtp` has the lines of the final fractures, with
the regions on the other side of each line (see "Output files"). A mesher
can mesh the fractures with these lines, and make boundary nodes from them,
for example the point arrays `boundary__<name>` of fehm-rs.
`tests/data/region_boundaries_mesh`, `tests/data/region_boundaries_family_region`
and `tests/data/region_boundaries_user_fractures` are examples.

- The region of a fracture:
  - A stochastic fracture: all cells of the mesh, thus the lines are at the
    mesh surface. The `region` of its family gives only its center.
  - A user-defined fracture: the cells of all its regions together (the arrays
    `region__<name>` of the fracture file). The faces between two of these
    cells give no line. A fracture with no region: all cells of the mesh.
- The plane of the fracture cuts each tetrahedron of these cells (see "Domain
  mesh"). The boundary of the part of the fracture in the cross-sections,
  without its edges on the fracture edge, gives the lines. Where the fracture
  crosses the mesh surface in the cells of other regions, these lines are
  also lines (with `mesh_surface` = 1). Thus each fracture which crosses the
  mesh surface has a line there. A line ends on the fracture edge or on
  another line (where the boundary of the region meets the mesh surface), or
  it is a closed loop around a hole of the region. A
  line has one or more straight segments; a segment ends where the line
  bends, where the cell on its other side changes, or where its boundary
  surfaces change. The vertices are on the cell faces.
- A segment on a surface cell of the mesh with a boundary surface (see
  "Domain mesh") has that boundary surface: its output array
  `boundary__<name>` is 1. The surface cells only give names to the lines:
  a surface cell inside the cells of the region of a fracture gives no line.
  `tests/data/boundary_surfaces` is an example.
- With FRAM, FRAM checks the lines of a new fracture (h = `min_feature_size`):
  - `ShortBoundaryLine`: a line is shorter than h.
  - `BoundaryCloseToEdge`: a vertex of a line which is not on the fracture
    edge is nearer than h to the fracture edge; or, where a line ends on the
    edge, the point at the distance h along the line is nearer than 0.4 h to
    the edge (the line meets the edge at an angle of less than about 24°).
  - `IntersectionCloseToBoundary`: an intersection is nearer than h to a line
    of one of its two fractures, and does not cross it; or it crosses the
    line nearer than h to an end of the intersection, or at an angle outside
    of 47° to 133° (the angle rule of a triple intersection).

  With `fram.relaxed = true`, only `ShortBoundaryLine`. The checks of the
  lines of a user-defined fracture (the first two) apply only with
  `fram.check_user_fractures = true`; thus by default each user-defined
  fracture in the domain is accepted, as without the option. A rejected
  stochastic fracture moves, as for the other FRAM rejects.
- Without FRAM, the network is the same as without the option: the option
  adds only the lines.
- The steps of a region boundary in a fine mesh give short segments. A short
  segment alone is not a reject; only a line shorter than h is.
- `rejections.dat` and the log have the number of rejects for the lines:
  `Boundary Lines of the Regions`.
- The lines are calculated at each attempt of a fracture. On a mesh with many
  cells, this takes the most time of the run.

### Triangle mesh for fehm-rs

With `mesh.enable = true`, the run also writes `dfn_mesh.vtu`: one triangle
mesh of the final fractures, the mesh file of [fehm-rs](../fehm-rs). The
network is the same as without the option.

- Only the part of each fracture in the domain mesh, and in the regions of
  the fracture, is in the mesh: a user-defined fracture with regions (its
  `region__<name>` arrays) has no triangle outside them; a stochastic
  fracture can cross the boundary of the `region` of its family. The fractures stay whole in the network; the
  mesh is cut along their boundary lines. The log gives the part of the
  area which is cut off.
- Fractures which intersect share the points of their intersection. Two
  fractures which only touch at a point share that point. The pieces of a
  user-defined fracture share the points of their common edges and corners.
  Thus each cluster is one connected mesh (the cut at the mesh surface can
  split a cluster).
- The mesh of each fracture follows its outline, its intersections and its
  boundary lines (constrained Delaunay triangulation of the crate `spade`).
- Triangle size: `min_feature_size` (h) along the intersections and the
  boundary lines, and h + `slope` × distance from them, at most
  `max_size_factor` × h, and at most the width of the fracture at the point
  (not below h / 2), thus a narrow fracture gets triangles of a good shape
  across it.
- Angles: a triangle with an angle below 20° gets a new point (Ruppert
  refinement), except where the small angle is an angle of the input (two
  lines of the mesh meet at it). A line of one fracture (outline, boundary
  line) is divided where needed; a shared line is divided in all its
  fractures, which are then meshed again (at most 10 rounds; a warning gives
  the edges which still need a division). No new point is nearer than 0.7 h
  to another point.
- Tiny features of user-defined fractures which FRAM did not check (FRAM off,
  or `fram.check_user_fractures = false`): their intersections shorter than
  h / 2 are point contacts; the points of their shared lines within h / 2 of
  each other are one point; their shared lines which are almost parallel and
  within h / 2 of each other get the same points; an outline edge with such
  a line or point within h / 2 goes through it. The corners of the pieces do
  not move. A line can then cross another line; the log gives the number of
  edges which the mesh does not follow for that reason.
- A sliver (a fracture narrower than h with a corner below 10°) has small
  angles near that corner: no refinement can change an angle of the input.
  The log gives the number of slivers.
- Smoothing (as LaGriT `smooth` and `recon` in dfnWorks, and the smoothing
  of the MeshLab filter "Isotropic Explicit Remeshing"): three passes; each
  point which is on no line, with a triangle around it with an angle below
  30° or above 90°, moves to the area-weighted centroid of its triangles if
  no triangle around it then gets a smaller smallest angle. The mesh stays
  Delaunay.
- A line edge with an angle above 90° opposite it (a negative connection
  coefficient in fehm-rs) is divided, a shared line in all its fractures,
  then the refinement and the smoothing again. Remaining negative
  coefficients are mostly in slivers.
- A point on a line is a point of that line (also of a shared line, in all
  its fractures). A triangle with no area (fehm-rs rejects it) is thus not
  made; if one is found, it is not in the mesh, with a warning.
- fehm-rs needs a boundary point in each cluster. A cluster with no point on
  a boundary array (see `dfn_mesh.vtu` below; for example, fully in one
  region and away from the mesh surface) is not in the mesh; the log gives
  the number of these fractures. The other output files keep them.
- A line with a cell of no region on the other side, and on no boundary
  surface, gives no array: its points are not boundary points.
- A point where a line of a region meets a line on the mesh surface (or of
  another region) is in both arrays. fehm-rs gives such a point the pressure
  of the first of the two boundaries in its input file, with a warning.

Apertures: each family has an `aperture` distribution, and each user-defined
fracture has the value of the cell data array `aperture` of its file (see
"User-defined fractures"). The distributions (as the transmissivity
options of dfnWorks, with the aperture in place of the transmissivity; r is
the x radius):

| `distribution` | Parameters | Aperture (m) |
|---|---|---|
| `constant` | `value` (> 0) | `value` |
| `log_normal` | `log_mean`, `sd` (>= 0) | exp(`log_mean` + `sd` × N(0, 1)) |
| `correlated` | `alpha` (> 0), `beta` | `alpha` × r^`beta` |
| `semi_correlated` | `alpha` (> 0), `beta`, `sd` (>= 0) | `alpha` × r^`beta` × exp(`sd` × N(0, 1)) |

The permeability is
aperture² / 12 (cubic law). The random numbers of the apertures come from a
separate generator with the seed, thus they do not change the network.

### Domain mesh

`domain = "domain.vtu"` gives the domain as a VTK UnstructuredGrid file.
A relative path starts at the folder of the input file. The key must be
before the first table (`[generation]`) of the input file. The mesh is the only
definition of the domain. `tests/data/mesh_domain` is an example, and
`tests/data/box.vtu` is a box domain (a mesh of one hexahedron).

- The file must be an ASCII `.vtu` file. In ParaView, select the data mode
  "Ascii" in the save dialog. In PyVista, use `mesh.save("domain.vtu", binary=False)`.
- The volume cells can be tetrahedra, voxels, hexahedra, wedges and pyramids
  (VTK types 10 to 14). The file can also have surface cells, triangles and
  quads (VTK types 5 and 9), for the boundary surfaces (see below). The
  surface cells have no effect on the volume, the bounding box, the regions
  or the fracture centers.
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
the third cell has no region (`tests/data/regions.vtu` is an example):

```xml
<CellData>
<DataArray type="Int32" Name="region__marble" format="ascii">1 0 0</DataArray>
<DataArray type="Int32" Name="region__north" format="ascii">1 1 0</DataArray>
</CellData>
```

Boundary surfaces of the mesh: each cell data array with a name
`boundary__<name>` (optional) is one boundary surface with the name `<name>`.
A surface cell with a value above 0 is on the boundary surface; a surface
cell can be on more than one boundary surface, and on none. The surface
cells must be on faces of the volume cells: on the mesh surface, or between
cells of two regions. With `[region_boundaries]`, a boundary line segment on
such a cell gets the boundary surface, and the output files have the array
`boundary__<name>` (see "Output files"); for fehm-rs, a boundary surface is a
boundary as a region is. In this example, the first cell is a hexahedron, the
second cell (a quad) is on the boundary surface `inflow`, and the third cell
(a quad) is on no boundary surface (`tests/data/boundary_surfaces` is an
example):

```xml
<CellData>
<DataArray type="Int32" Name="boundary__inflow" format="ascii">0 1 0</DataArray>
</CellData>
```

These are errors: a volume cell with a value above 0 in a `boundary__<name>`
array, a surface cell with a value above 0 in a `region__<name>` array, a
boundary surface with the name of a region (both give the output array
`boundary__<name>`), and the name `mesh_surface`. The value 0 is allowed
in both kinds of array (ParaView writes 0 for a missing value when it merges
two meshes).

`clusters.keep_only_connected_to_surfaces = true` keeps only the clusters
which touch each boundary surface of the mesh (as `boundaryFaces` of
dfnGen, with all surfaces selected). A cluster touches a surface if a
boundary line segment of one of its fractures is on it. Thus the option
needs `[region_boundaries]` with `enable = true` and a mesh with boundary
surfaces; else it is an error. `keep_only_largest` then keeps the largest
of these clusters. If there are clusters but none touches each surface, all
clusters stay (as in dfnGen), with a warning in the log, and the library
result `Network` has `is_without_cluster_on_all_surfaces = true`.

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
`tests/data/user_fractures` is an example.

- The file must be an ASCII `.vtp` file. In ParaView, select the data mode
  "Ascii" in the save dialog. In PyVista, use `mesh.save("fractures.vtp", binary=False)`.
- Each polygon cell is one fracture. Its points are the vertices of the
  fracture, in the order of the cell. Give points in one plane and a convex
  polygon; the code does not check this.
- A polygon must have 3 points or more. The fractures are inserted in the
  order of the cells.
- The coordinates are those of the mesh file.
- The cell data array `aperture` (necessary) gives the aperture (m) of each
  polygon: one number above 0 for each cell of the file. The pieces of one
  fracture (same `fracture_id`, see below) must have the same aperture.
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
- `fram.check_user_fractures = false` (the default) turns off FRAM for an intersection of
  two user-defined fractures. Then each user-defined fracture in the domain
  (and in its region) is accepted. The intersection is recorded as FRAM finds
  it: it is not shortened, and such a pair gives no triple intersection point.
  An intersection of a stochastic fracture with a user-defined fracture keeps
  all FRAM checks. The network then does not obey the minimum feature size `min_feature_size`
  at the user-defined fractures.
- The point data and the other cell data of the file are ignored. Thus
  `fractures.vtp` of a run can be the user-defined fractures of a new run,
  after an `aperture` array is added to it. It
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
of a rectangle is half the side). Angles in radians are changed to degrees.
The box mesh has no boundary surfaces, thus `boundaryFaces` is not converted
(`keep_only_connected_to_surfaces = false`). The ellipses do not have the same vertices
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

The numbers of fractures in the statistics of the log count the pieces of a
user-defined fracture (polygons with the same `fracture_id`) as one fracture.
Each number also gives the number of polygons, for example
`Final Number of Fractures: 1357 (5003 polygons)`. P30 uses the number of
fractures. The section of the user-defined fractures gives the accepted, the
rejected (no polygon accepted) and the final user-defined fractures. The
warning for a rejected user-defined fracture and `userFractureRejections.dat`
give the number of the polygon in its file.

During the insertion of the stochastic fractures, a progress bar is shown on
stderr, if stderr is a terminal:

```
[00:00:12] [################>-----------------------]  42 %  25/60 fractures  ETA 00:00:17
```

Without `n_poly`, the percent is the part of the total target P32
which the families have, and the bar shows the number of accepted fractures
without a total. The large fractures are inserted first, thus the remaining time
of such a run is approximate.

| File | Content |
|---|---|
| `fractures.vtp` | The final fractures as polygons with their properties (VTK PolyData, for ParaView or PyVista). See below. |
| `fractures_removed.vtp` | The accepted fractures which are not final fractures, as polygons with their properties (VTK PolyData). Written only if a fracture was removed. See below. |
| `connectivity.dat` | Line N: numbers of the final fractures which intersect final fracture N. |
| `intersection_list.dat` | One line for each intersection: the two fracture numbers, mid point and length. |
| `triple_points.dat` | Coordinates of the triple intersection points. Written only if FRAM is on and `fram.triple_intersections = true`. |
| `rejections.dat` | Number of rejected fractures for each rejection reason. |
| `rejectsPerAttempt.dat` | Line N: number of rejected attempts before accepted fracture N. |
| `userFractureRejections.dat` | Number of each rejected user-defined fracture (the position of its polygon in the `.vtp` file, from 1). Written only if a user-defined fracture was rejected. |
| `fractures_rejected.vtp` | The rejected fractures as polygons with their properties (VTK PolyData). Written only if `fram.write_rejected_fractures = true`. See below. |
| `dfn_mesh.vtu` | The triangle mesh of the final fractures for fehm-rs (VTK UnstructuredGrid). Written only if `mesh.enable = true`. See below. |
| `fracture_boundaries.vtp` | The boundary lines of the regions on the final fractures, as line segments (VTK PolyData). Written only if `region_boundaries.enable = true`. See below. |

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
it has no re-translation left (`fram.rejects_per_fracture`). A user-defined
fracture has only one attempt. The file can be very large.

A fracture has no cell if its last attempt was outside the domain, because it
then has no valid polygon.

| Array | Content |
|---|---|
| `family` | Same numbers as in `fractures.vtp`. |
| `user_fracture_id` | Same values as in `fractures.vtp`. |
| `reject_reason` | Reason for the rejection of the last attempt, as text (VTK `String` array; ParaView can color by it). See the table below. |
| `x_radius`, `y_radius`, `aspect_ratio` | Size. |
| `normal` | Unit normal vector. |
| `translation` | Translation of the fracture center from the origin. |

| `reject_reason` | Reason |
|---|---|
| `ShortIntersection` | Intersection shorter than `min_feature_size`. |
| `CloseToIntersection` | Intersection too near to another intersection. |
| `CloseToEdge` | Fracture too near to an edge of another fracture. |
| `SamePlane` | Intersecting fractures are in the same plane. |
| `TripleIntersectionNotAllowed` | Triple intersection, and `triple_intersections = false`. |
| `CloseToPreviousIntersection` | Intersection too near to a previous intersection. |
| `TripleAngleTooSmall` | Triple intersection with an angle that is too small. |
| `TriplePointCloseToEndpoint` | Triple intersection point too near to an end point of an intersection. |
| `TriplePointCloseToTriplePoint` | Triple intersection point too near to another triple intersection point. |
| `ShortBoundaryLine` | A boundary line of the region on the fracture is shorter than `min_feature_size`. |
| `BoundaryCloseToEdge` | A boundary line of the region is too near to the fracture edge. |
| `IntersectionCloseToBoundary` | An intersection is too near to a boundary line of the region, or crosses it near its end or at a small angle. |

The log gives the same name for each FRAM reject. The reject of a
user-defined fracture is a warning with the number of the fracture. The
reject of each attempt of a stochastic fracture is at the level `debug`; to
see it, run with `RUST_LOG=debug` (or `RUST_LOG=dfngen::fram=debug` for only
these lines).

`fractures.vtp`, `fractures_removed.vtp` and `fractures_rejected.vtp` together
do not have all generated fractures. An accepted fracture which is removed later because it
is smaller than `remove_fractures_less_than` is in none of these files.

### Cell data of `fracture_boundaries.vtp`

The file has one line cell for each segment of a boundary line of a final
fracture, in the order of the fractures and of their lines. The points of a
segment are in the plane of the fracture; the ends of a line are on the
fracture edge in `fractures.vtp`.

| Array | Content |
|---|---|
| `fracture_id` | Number of the final fracture, as in `fractures.vtp`. |
| `line_id` | Number of the boundary line on its fracture, starts at 1. The segments of one line are in the order of the line. |
| `mesh_surface` | 1 if no cell of the mesh is on the other side of the segment (the segment is on the surface of the mesh), else 0. |
| `boundary__<name>` | One array for each region of the mesh, in the order of the mesh file: 1 if the cell on the other side of the segment is in the region `<name>`, else 0. All arrays are 0 for a segment with a cell of no region on the other side. |
| `boundary__<name>` | Then one array for each boundary surface of the mesh, in the order of the mesh file: 1 if the segment is on a surface cell of the boundary surface `<name>`, else 0. |

For example, a segment of a fracture of the region `south` with 1 in
`boundary__north` is on the boundary between the regions `south` and `north`.

### Data of `dfn_mesh.vtu`

The points are in the coordinates of the domain mesh. Each cell is a
triangle (VTK type 5), counter-clockwise around the normal of its fracture.

| Array | Content |
|---|---|
| `fracture_id` (cell) | Number of the final fracture, as in `fractures.vtp`. |
| `aperture` (cell) | Aperture of the fracture (m). |
| `permeability` (cell) | aperture² / 12 (m²). |
| `boundary__mesh_surface` (point) | 1 for a point on a boundary line on the mesh surface, else 0. |
| `boundary__<name>` (point) | One array for each region of the mesh: 1 for a point on a boundary line with a cell of the region `<name>` on the other side, else 0. |
| `boundary__<name>` (point) | Then one array for each boundary surface of the mesh: 1 for a point on a boundary line on a surface cell of the boundary surface `<name>`, else 0. |

In the input file of fehm-rs, each `[[boundary]]` gives the pressure of one of
these arrays by its name (`mesh_surface` or `<name>`).

## Library

The package is also the library `dfngen`. The `dfngen` binary is a small
program which uses it. The library reads no file and writes no file, unless
you call `read_input` or `write_output_to_folder`:

```rust
use std::io::{self, Write};

use dfngen::{build_input_from_text, generate, read_input_file_names, write_output, DfngenError};

fn generate_from_text(toml_text: &str, load: impl Fn(&str) -> String) -> Result<(), DfngenError> {
    // The TOML input names the mesh file and the user-defined fracture file
    let names = read_input_file_names(toml_text)?;
    let mesh_text = load(&names.domain);
    let fractures_text = names.user_fractures.as_deref().map(&load);

    let (input, families) =
        build_input_from_text(toml_text, &mesh_text, fractures_text.as_deref())?;

    let network = generate(input, families, &mut |progress| {
        eprintln!("{} fractures, {:.0} %", progress.accepted_count, 100. * progress.fraction);
    })?;

    // `create_file` gives the writer of each output file (see "Output files")
    write_output(&network, &mut |file_name| {
        println!("--- {file_name}");
        Ok(Box::new(io::stdout()) as Box<dyn Write>)
    })?;

    Ok(())
}
```

- `read_input(path)` does the same as `read_input_file_names` and
  `build_input_from_text`, with the files of the file system.
- `generate` uses `seed` as it is. Only the binary replaces `seed = 0` with
  the system time.
- `generate` calls the progress function after each accepted stochastic
  fracture. `fraction` goes from 0 to 1.
- `Network` has the input, the families, the accepted and rejected fractures,
  the intersections, the statistics, and the indices of the final fractures.
- The library writes its log with the `tracing` crate. Your program installs
  the subscriber.
- Only the items at the crate root are a stable interface. The modules are
  public, but they can change.

Without the default feature `cli`, the library has no dependency on `clap`,
`indicatif` and `tracing-subscriber`, and it compiles for WebAssembly:

```toml
dfngen-rs = { version = "0.11", default-features = false }
```

```bash
cargo check --lib --no-default-features --target wasm32-unknown-unknown
```
