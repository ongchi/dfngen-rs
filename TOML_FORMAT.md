# TOML Configuration Format for dfngen-rs

## Overview

dfngen-rs now supports TOML (Tom's Obvious, Minimal Language) as a modern, structured alternative to the original custom text-based configuration format. TOML provides:

- **Better Readability**: Clear hierarchical structure with section headings
- **Type Safety**: Automatic type validation during parsing
- **Error Messages**: Specific line numbers and error descriptions on parse failures
- **Industry Standard**: Wide tool support and documentation
- **Less Error-Prone**: No ambiguity about value formats or ordering

## File Format

Configuration files use the `.toml` extension:
```bash
dfngen input.toml output/
```

## Configuration Sections

### [domain]
Defines the 3D domain for fracture generation.

**Parameters:**
- `size` (required): Domain dimensions `{x, y, z}` in meters
- `h` (required): Minimum feature size (FRAM parameter)
- `size_increase` (optional): Size extension for fractures beyond boundaries `{x, y, z}`

**Example:**
```toml
[domain]
size = { x = 100, y = 100, z = 100 }
h = 1.0
size_increase = { x = 0, y = 0, z = 0 }
```

### [stopping_condition]
Defines when fracture generation stops.

**Parameters:**
- `mode` (required): `"nPoly"` (generate N fractures) or `"P32"` (generate until intensity target)
- `n_poly` (required for nPoly mode): Number of fractures to generate
- `p32_target` (required for P32 mode): Fracture intensity target per family
- `radii_list_increase` (optional): Percentage increase to radii list (default: 0.2 = 20%)

**Example:**
```toml
[stopping_condition]
mode = "nPoly"
n_poly = 1000
radii_list_increase = 0.2
```

### [fram]
Configures the FRAM (Feature Rejection Algorithm for Meshing) algorithm.

**Parameters:**
- `enable` (optional, default: true): Enable/disable FRAM
- `relaxed` (optional, default: false): Use relaxed FRAM (less strict)
- `triple_intersections` (optional, default: true): Accept/reject triple intersections
- `visualization_mode` (optional, default: false): Only first triangulation round

**Example:**
```toml
[fram]
enable = true
relaxed = false
triple_intersections = true
visualization_mode = false
```

### [output]
Configures output file generation.

**Parameters:**
- `print_reject_reasons` (optional): Print rejection reasons to console
- `output_all_radii` (optional): Output all generated radii (accepted + rejected)
- `output_accepted_radii_per_family` (optional): Output radii before removal of isolated fractures
- `output_final_radii_per_family` (optional): Output radii after removal of isolated fractures

**Example:**
```toml
[output]
print_reject_reasons = false
output_all_radii = false
output_accepted_radii_per_family = false
output_final_radii_per_family = false
```

### [fractures]
Configures general fracture generation parameters.

**Parameters:**
- `rejects_per_fracture` (optional, default: 1): Number of re-translation attempts
- `force_large_fractures` (optional): Force insertion of largest fractures per family
- `remove_smaller_than` (optional): Remove fractures below this radius after generation
- `orientation_option` (optional, default: "spherical"): Angle definition system
  - `"spherical"`: Theta/Phi angles
  - `"trend_plunge"`: Trend/Plunge angles
  - `"dip_strike"`: Dip/Strike angles (RHR)

**Example:**
```toml
[fractures]
rejects_per_fracture = 1
force_large_fractures = false
remove_smaller_than = 0.0
orientation_option = "spherical"
```

### [boundaries]
Configures boundary face constraints.

**Parameters:**
- `faces` (required): Array of 6 booleans `[+x, -x, +y, -y, +z, -z]`
- `keep_only_largest_cluster` (optional): Keep only largest cluster meeting boundary constraints
- `keep_isolated_fractures` (optional, default: true): Keep fractures not connected to boundaries
- `ignore_boundary_faces` (optional): Ignore boundary constraints entirely

**Example:**
```toml
[boundaries]
faces = { faces = [true, true, true, true, true, true] }
keep_only_largest_cluster = false
keep_isolated_fractures = true
ignore_boundary_faces = false
```

### [[ellipse_families]] and [[rectangle_families]]
Defines stochastic fracture families (array sections, can have multiple).

**Common Parameters:**
- `name` (optional): Family name for reference
- `probability` (required): Probability weight for family selection
- `aspect_ratio` (required): y-radius / x-radius ratio
- `beta` (required): Rotation angle about normal (degrees)
- `beta_distribution` (required): Beta is distributed (true) or constant (false)
- `layer` (optional, default: 0): Layer ID (0=whole domain)
- `region` (optional, default: 0): Region ID (0=whole domain)

**Ellipse-Specific:**
- `num_points` (optional, default: 8): Vertices for ellipse approximation

**Sub-sections:**
- `[..families.orientation]`: Angle distribution parameters
  - `angle1`: Primary angle (theta/trend/dip)
  - `angle2`: Secondary angle (phi/plunge/strike)
  - `kappa`: Concentration parameter (Fisher distribution)

- `[..families.radius]`: Radius distribution parameters
  - `distribution_type`: `"lognormal"`, `"power_law"`, `"exponential"`, or `"constant"`
  - Distribution-specific parameters (see below)

**Radius Distribution Types:**

*Lognormal:*
```toml
[ellipse_families.radius]
distribution_type = "lognormal"
log_mean = 1.0
log_std = 0.5
log_min = 0.5
log_max = 5.0
```

*Power-Law:*
```toml
[ellipse_families.radius]
distribution_type = "power_law"
alpha = 2.0
min = 0.5
max = 10.0
```

*Exponential:*
```toml
[ellipse_families.radius]
distribution_type = "exponential"
exp_mean = 2.0
exp_min = 0.5
exp_max = 10.0
```

*Constant:*
```toml
[ellipse_families.radius]
distribution_type = "constant"
constant = 2.5
```

**Example Family:**
```toml
[[ellipse_families]]
name = "Ellipse Family 1"
probability = 0.6
num_points = 8
aspect_ratio = 0.5
beta = 0.0
beta_distribution = false
layer = 0
region = 0
p32_target = 5.0

[ellipse_families.orientation]
angle1 = 0.5
angle2 = 0.5
kappa = 0.5

[ellipse_families.radius]
distribution_type = "lognormal"
log_mean = 1.0
log_std = 0.5
log_min = 0.5
log_max = 5.0
```

### [optional]
Optional configuration section containing seed and advanced features.

**Parameters:**
- `seed` (optional): Random seed for reproducible results (0=system time)

**Subsections:**
- `[optional.layers]`: Layer boundary definition
  - `boundaries`: Array of z-values defining layer boundaries
- `[optional.regions]`: Region boundary definition
  - `boundaries`: Array of coordinates [xmin, xmax, ymin, ymax, zmin, zmax, ...]
- `[optional.polygon_boundary]`: 2D polygon domain boundary
  - `vertices_file`: Path to file with boundary vertices

**Example:**
```toml
[optional]
seed = 12345

[optional.layers]
boundaries = [10, 20, 30, 50, 70, 80]

[optional.regions]
boundaries = [0, 100, 0, 100, 0, 100]

[optional.polygon_boundary]
vertices_file = "domain_boundary.txt"
```

### [user_fractures]
Configures user-defined fractures.

**Parameters:**
- `ellipses_file` (optional): Path to user ellipses file
- `rectangles_file` (optional): Path to user rectangles file
- `polygons_by_coord_file` (optional): Path to polygon coordinates file
- `ellipses_by_coord_file` (optional): Path to ellipse coordinates file
- `rectangles_by_coord_file` (optional): Path to rectangle coordinates file
- `insert_rectangles_first` (optional): Insert rectangles before ellipses

**Example:**
```toml
[user_fractures]
ellipses_file = "user_ellipses.txt"
rectangles_file = "user_rectangles.txt"
insert_rectangles_first = true
```

## Complete Example

See `examples/example_config.toml` for a fully commented example configuration file.

## Migration from Text Format

The legacy text format is still supported. To convert:

1. **Domain section**: `domainSize:` → `[domain]` section with `size`
2. **Stopping conditions**: `stopCondition:` and `nPoly:` → `[stopping_condition]`
3. **FRAM settings**: `disableFram:`, `rFram:` → `[fram]` with `enable` (inverted)
4. **Families**: Flat list of parameters → Array sections `[[ellipse_families]]` or `[[rectangle_families]]`

## Error Handling

If the TOML file is invalid, you'll see:
```
Failed to parse TOML file 'input.toml': error at line X: [detailed error message]
```

Common errors:
- **Duplicate sections**: `[[family]]` declared twice without array
- **Missing required fields**: Check section documentation
- **Type mismatches**: `size = 100` should be `size = { x = 100, y = 100, z = 100 }`
- **Invalid keys**: Typos in configuration keys

## Notes

- All coordinates and sizes are in meters
- All angles are in degrees (internally converted to radians)
- Comments starting with `#` are supported
- Array values use square brackets `[...]` or table syntax `{...}`
- Probabilities should sum to ~1.0 across families (normalized automatically)
