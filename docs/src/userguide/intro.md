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
runsim(config, velmod)
```

The configuration must set `output_file` to a path ending in `.h5`. `runsim` writes all results to
that file and returns `nothing` — results can be retrieved afterwards with `load_results`.

### Example — passing objects directly

```julia
using ElasticFDSG

velmod = zeros(7, 200, 200)   # fill with your data
# ... (see Velocity Models for details)

config = Dict(
    "settings" => Dict(
        "device" => "cpu", "precision" => "Float32",
        "spatial_derivative_order" => 4, "verbose" => true,
        "output_file" => "path/to/output.h5",
    ),
    # ... (see Configurations for the full schema)
)

runsim(config, velmod)
results = load_results(config["settings"]["output_file"])
```

### Example — passing file paths

```julia
using ElasticFDSG

runsim("path/to/config.yaml", "path/to/velmod.jld2")
results = load_results("path/to/output.h5")
```

The results dictionary mirrors the HDF5 group structure.