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
visualization_mode = false

[clusters]
boundary_faces = [true, true, false, false, false, false]  # +x, -x, +y, -y, +z, -z
ignore_boundary_faces = true
keep_only_largest = false
keep_isolated_fractures = false

[output]
print_reject_reasons = false
all_radii = false
final_radii_per_family = false
accepted_radii_per_family = false

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
