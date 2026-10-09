# Upstream examples that dfngen-rs cannot run yet

These inputs were exported from the dfnWorks examples with
`scripts/export_dfngen_input.py`. `scripts/compare_output.sh` does not run them,
because it reads only the `.dat` files directly in `tests/data`.

Run one from its folder: `cd <case> && dfngen input.dat <output_folder>`.

| Case | Upstream driver | Result (binary at f02d8ab) | Cause |
|------|-----------------|----------------------------|-------|
| `calibrate_p32_from_p10` | `make_network.py` | Panic in `zip_eq` | The user-fracture file has `[x, y, z]` rows and one `AngleOption` line for each fracture. The reader expects `(x,y,z)` rows and one value. |
| `dfm_4_frac` | `driver.py` | Panic in `zip_eq` | Same as above. |
| `pruning` | `driver.py` | Does not end | P32 mode with 2 families. The log repeats "P32 For Family 1 Completed". |
| `finite_matrix_diffusion` | `forsmark_percolation_comparison.py` | Does not end in 60 s | Not known. P32 mode with 3 families in a 1000 m domain. This can be the same loop as `pruning`, or a long run. The family probabilities in the input sum to 2.43. |

Move a case to `tests/data/upstream_<name>.dat` when it runs.
