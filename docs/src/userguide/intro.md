# General Usage

Running a simulation requires two inputs:

1. A **velocity model** — a Julia array (or file path to `.jld2` / `.npy` / `.npz`).
2. A **configuration** — a Julia `Dict` (or file path to a `.yaml` file).

The simulation dimension (2D or 3D) is **automatically detected** from the shape of the velocity model array:
a 3-dimensional array (`7 × nx × nz`) triggers a 2D simulation,
and a 4-dimensional array (`13 × nx × ny × nz`) triggers a 3D simulation.

## Running a simulation

```julia
using ElasticFDSG

# config and velmod can each be a file path (String) or a Julia object (Dict / Array)
result = runsim(config, velmod)
```

If the configuration has `output_file = nothing`, `runsim` returns the populated `FDSG` struct directly, which contains all field data.
If an `output_file` path is provided (must end in `.h5`), the results are written to HDF5 and `runsim` returns `nothing`.

### Example — passing objects directly

```julia
using ElasticFDSG

velmod = zeros(7, 200, 200)   # fill with your data
# ... (see Velocity Models for details)

config = config_template_2d(
    device = "cpu", precision = "Float32",
    # ... (see Configurations for details)
)

fdsg = runsim(config, velmod)
```

### Example — passing file paths

```julia
using ElasticFDSG

fdsg = runsim("path/to/config.yaml", "path/to/velmod.jld2")
```

## Loading results from HDF5

When results are saved to disk, they can be loaded back into a nested Julia dictionary:

```julia
using ElasticFDSG

data = load_results("path/to/output.h5")
```

The dictionary mirrors the HDF5 group structure:

| Key | Contents |
|-----|----------|
| `"grid"` | `x_coordinates`, (`y_coordinates`,) `z_coordinates` — inner-domain coordinate vectors |
| `"time"` | `t0`, `tend`, `dt`, `time` vector |
| `"source"` | wavelet, moment tensor components, source location |
| `"geophones"` | `geophone_i/data` — `(ncomp, nt)` array; `geophone_i/location` |
| `"das"` | `x_aligned/fiber_i/data` — `(nch, nt)`; analogous for `z_aligned` (and `y_aligned` in 3D) |
| `"snapshots"` | `XZ` — `(ntime, nfields, nx, nz)` in 2D; per-plane datasets in 3D |

## Examples

A fully worked 2D example is available in
[`examples/demo2d.ipynb`](https://github.com/wtegtow/ElasticFDSG.jl/tree/main/examples).
A 3D example (`demo3d`) is planned and will be added to the same folder.
