# Configurations

Configurations are plain Julia dictionaries or paths to `.yaml` files that follow a fixed schema.

Top-level keys `"settings"`, `"time"`, `"source"`, `"boundaries"`, and `"receivers"` are all required.
Within `"receivers"`, the `"geophones"`, `"das"` (and its `"x_aligned"` / `"y_aligned"` / `"z_aligned"`
sub-lists), and `"snapshots"` entries are all optional — simply omit any that are not needed.

---

## 2D configuration

```julia
using ElasticFDSG

config = Dict(
    "settings" => Dict(
        "device" => "cpu",                     # "cpu" | "cuda" | "metal" | "amd" | "oneapi"
        "precision" => "Float32",              # "Float32" | "Float64"
        "spatial_derivative_order" => 4,       # stencil half-width 1–10
        "verbose" => true,                     # print simulation summary and progress
        "output_file" => "path/to/output.h5",  # must end in .h5
    ),
    "time" => Dict(
        "start" => 0.0,
        "end" => 1.0,
        "timestep" => 0.001,                   # reduced automatically if CFL is violated
    ),
    "source" => Dict(
        "dominant_frequency" => 30.0,          # [Hz]
        "wavelet_type" => "ricker",            # "ricker" | "gauss1d"
        "wavelet_center" => 0.05,              # peak time of the wavelet [s] (≥ 1.25/fdom)
        "seismic_moment" => 1e9,               # scalar seismic moment M₀ [N·m]
        "location" => Dict(
            "x" => 500.0,
            "z" => 500.0,
        ),
        "moment_tensor" => Dict(
            "Mxx" => 0.0,
            "Mxz" => 1.0,
            "Mzz" => 0.0,
        ),
    ),
    "boundaries" => Dict(
        # "absorbing" | "none"
        "xstart" => "absorbing", "xend" => "absorbing",
        "zstart" => "none",      "zend" => "absorbing",     # reflecting top
        "pml_layer" => 10,                                  
    ),
    "receivers" => Dict(
        # Geophones — list of dicts with x,z locations
        "geophones" => [
            Dict("x" => 800.0, "z" => 300.0),
            Dict("x" => 900.0, "z" => 300.0),
            # add more here ...
        ],

        # DAS — axis-aligned strain-rate profiles
        "das" => Dict(
            # x_aligned: list of dicts for fibers running along x at fixed z
            "x_aligned" => [
                Dict("x" => Dict("start"=>100.0, "step"=>5.0, "end"=>900.0), "z"=>400.0),
                # add more here ...
            ],
            # z_aligned: list of dicts for fibers running along z at fixed x 
            "z_aligned" =>[
                Dict("x" => 100.0, "z"=>Dict("start"=>100.0, "step"=>5.0, "end"=>900.0)),
                # add more here ...
            ],
        ),
        "snapshots" => Dict(
            "times" => [0.25, 0.5, 0.75, 1.0],
            "fields" => ["vx", "vz", "sxx", "sxz", "szz"],
        ),
    ),
)
```

### Available snapshot field names (2D)

| Name | Description |
|------|-------------|
| `"vx"` | Particle velocity, x-component |
| `"vz"` | Particle velocity, z-component |
| `"sxx"` | Normal stress $\sigma_{xx}$ |
| `"szz"` | Normal stress $\sigma_{zz}$ |
| `"sxz"` | Shear stress $\sigma_{xz}$ |

---

## 3D configuration

```julia
using ElasticFDSG

config = Dict(
    "settings" => Dict(
        "device" => "cuda",
        "precision" => "Float32",
        "spatial_derivative_order" => 4,
        "verbose" => true,
        "output_file" => "path/to/output.h5",
    ),
    "time" => Dict(
        "start" => 0.0,
        "end" => 0.8,
        "timestep" => 0.0005,
    ),
    "source" => Dict(
        "dominant_frequency" => 40.0,
        "wavelet_type" => "ricker",
        "wavelet_center" => 0.04,
        "seismic_moment" => 1e10,
        "location" => Dict(
            "x" => 500.0,
            "y" => 125.0,
            "z" => 250.0,
        ),
        "moment_tensor" => Dict(
            "Mxx" => -1.0, "Mxy" => 0.0, "Mxz" => 0.0,
            "Myy" =>  0.0, "Myz" => 0.0, "Mzz" => 1.0,
        ),
    ),
    "boundaries" => Dict(
        "xstart" => "absorbing", "xend" => "absorbing",
        "ystart" => "absorbing", "yend" => "absorbing",
        "zstart" => "none",      "zend" => "absorbing",
        "pml_layer" => 10,
    ),
    "receivers" => Dict(
        "geophones" => [
            Dict("x"=>950.0, "y"=>20.0, "z"=>250.0),
            Dict("x"=>750.0, "y"=>20.0, "z"=>250.0),
            Dict("x"=>550.0, "y"=>20.0, "z"=>250.0),
            # add more here ...
        ],
        "das" => Dict(
            "x_aligned" => [
                Dict("x"=>Dict("start"=>0.0, "step"=>5.0, "end"=>500.0), "y"=>50.0, "z"=>500.0),
                # add more here ...
            ],
            "y_aligned" => [
                Dict("x"=>950.0, "y"=>Dict("start"=>0.0, "step"=>5.0, "end"=>500.0), "z"=>500.0)
                # add more here ...
            ],
            "z_aligned" => [
                Dict("x"=>950.0, "y"=>50.0, "z"=>Dict("start"=>0.0, "step"=>5.0, "end"=>500.0)),
                Dict("x"=>250.0, "y"=>50.0, "z"=>Dict("start"=>0.0, "step"=>5.0, "end"=>500.0)),
                # add more here ...
            ],
        ),
        "snapshots" => Dict(
            # one entry per centre point; each centre produces an XY-, XZ-, and YZ-plane snapshot
            "plane_positions" => [
                Dict("x"=>500.0, "y"=>125.0, "z"=>250.0),
                Dict("x"=>250.0, "y"=>250.0, "z"=>250.0),
                # ...
            ],
            "times" => [0.4, 0.8],
            "fields" => ["vx", "vy", "vz"],
        ),
    ),
)

```

### Available snapshot field names (3D)

| Name | Description |
|------|-------------|
| `"vx"` / `"vy"` / `"vz"` | Particle velocity components |
| `"sxx"` / `"syy"` / `"szz"` | Normal stress components |
| `"sxy"` / `"sxz"` / `"syz"` | Shear stress components |

---

## Using a YAML file

Instead of constructing the dict in Julia you can write a `.yaml` file and pass its path to `runsim`:

```yaml
# config2d.yaml
settings:
    device: cpu
    precision: Float32
    spatial_derivative_order: 4
    verbose: true
    output_file: /path/to/output.h5

time:
    start: 0.0
    end:   1.0
    timestep: 0.001

source:
    dominant_frequency: 30
    wavelet_type: ricker
    wavelet_center: 0.05
    seismic_moment: 1.0e9
    location:
        x: 500.0
        z: 500.0
    moment_tensor:
        Mxx: 0.0
        Mxz: 1.0
        Mzz: 0.0

boundaries:
    xstart: absorbing
    xend:   absorbing
    zstart: none
    zend:   absorbing
    pml_layer: 10

receivers:
    geophones:
        - { x: 800.0, z: 300.0 }

    das:
        x_aligned:
            - { x: { start: 100, step: 5, end: 900 }, z: 400 }
        z_aligned: []

    snapshots:
        fields: ["vx", "vz"]
        times:  [0.5, 1.0]
```

---

## Stability and discretisation guidelines

The solver automatically checks the **CFL condition** and reduces `dt` if necessary.
For spatial discretisation, a common rule of thumb is:

$$\Delta x \leq \frac{V_\mathrm{min}}{10 \, f_\mathrm{max}}$$

where $V_\mathrm{min}$ is the minimum phase velocity in the model and $f_\mathrm{max}$ is the maximum
frequency content of the wavelet.

The PML thickness should be at least 10 grid cells; thicker layers improve absorption.
For strongly anisotropic media, PML instabilities may occur for certain parameter combinations —
in such cases, increasing the PML thickness often mitigates the issue (see [Method](../method.md)).
