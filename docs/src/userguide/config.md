# Configurations

Configurations are plain Julia dictionaries that follow a fixed schema.
They can be created programmatically using the helper functions [`config_template_2d`](@ref) and [`config_template_3d`](@ref).
These functions require all necessary information as keyword arguments.
Below are two examples of how to use the template functions.

Alternativly, configurations can be written by hand as YAML files.

---

## 2D configuration

```julia
using ElasticFDSG

config = config_template_2d(
    # ── Settings ──────────────────────────────────────────────────────────
    device    = "cpu",          # "cpu" | "cuda" | "metal" | "amd" | "intel"
    precision = "Float32",      # "Float32" | "Float64"
    fd_order  = 4,              # stencil half-width 1–10  
    verbose   = true,           # print simulation summary and progress
    output_file = nothing,      # String path ending in .h5, or nothing

    # ── Time ──────────────────────────────────────────────────────────────
    t_start = 0.0,
    t_end   = 1.0,
    dt      = 0.001,            # will be reduced automatically if CFL is violated

    # ── Source ────────────────────────────────────────────────────────────
    fdom           = 30.0,      # dominant frequency [Hz]
    wavelet        = "ricker",  # "ricker" | "gauss1d"
    wavelet_center = 0.05,      # peak time of the wavelet [s]  (≥ 1.25/fdom)
    seismic_moment = 1e9,       # scalar seismic moment M₀ [N·m]
    src_x = 500.0,              # source x-coordinate [m]
    src_z = 500.0,              # source z-coordinate [m]
    # 2D moment tensor components (symmetric, z is depth axis)
    Mxx = 0.0, Mxz = 1.0, Mzz = 0.0,
    anisotropic = false,        # if true, use anisotropic source radiation

    # ── Boundaries ────────────────────────────────────────────────────────
    # "absorbing" | "else"
    xstart = "absorbing", xend = "absorbing",  
    zstart = "else",      zend = "absorbing",  # reflecting surface at top
    pml_layer = 10,       # number of PML grid cells per absorbing boundary

    # ── Receivers ─────────────────────────────────────────────────────────
    # Geophones — [list of dicts] with x,z locations
    geophones = [
        Dict("x" => 800.0, "z" => 300.0),
        Dict("x" => 900.0, "z" => 300.0),
    ],

    # DAS — axis-aligned strain profiles
    # x_aligned: [list of dicts] fibers running along x at fixed z
    das_x_aligned = [
        Dict("x" => Dict("start"=>100.0, "step"=>5.0, "end"=>900.0), "z"=>400.0),
    ],
    # z_aligned: [list of dicts] fibers running along z at fixed x
    das_z_aligned = [], # if no receiver is needed, pass an empty list 

    # Snapshots
    snapshot_times  = [0.25, 0.5, 0.75, 1.0],          # times [s] to snapshot
    snapshot_fields = ["vx", "vz", "sxx", "sxz", "szz"],
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

config = config_template_3d(
    # ── Settings ──────────────────────────────────────────────────────────
    device    = "cuda",
    precision = "Float32",
    fd_order  = 4,
    verbose   = true,
    output_file = "/path/to/output.h5",

    # ── Time ──────────────────────────────────────────────────────────────
    t_start = 0.0,
    t_end   = 0.8,
    dt      = 0.0005,

    # ── Source ────────────────────────────────────────────────────────────
    fdom           = 40.0,
    wavelet        = "ricker",
    wavelet_center = 0.04,
    seismic_moment = 1e10,
    src_x = 500.0, src_y = 125.0, src_z = 250.0,
    # 3D moment tensor
    Mxx = -1.0, Mxy = 0.0, Mxz = 0.0,
    Myy =  0.0, Myz = 0.0, Mzz = 1.0,
    anisotropic = false,

    # ── Boundaries ────────────────────────────────────────────────────────
    xstart = "absorbing", xend = "absorbing",
    ystart = "absorbing", yend = "absorbing",
    zstart = "free",      zend = "absorbing",
    pml_layer = 10,

    # ── Receivers ─────────────────────────────────────────────────────────
    geophones = [
        Dict("x"=>950.0, "y"=>20.0, "z"=>250.0),
        Dict("x"=>750.0, "y"=>20.0, "z"=>250.0),
        Dict("x"=>550.0, "y"=>20.0, "z"=>250.0),
        # ...
    ],
    das_x_aligned = [
        Dict("x"=>Dict("start"=>0, "step"=>5.0, "end"=>"500.0"), "y"=>50.0, "z"=>500.0)
        # ...
    ],
    das_y_aligned = [],
    das_z_aligned = [
        Dict("x"=>950.0, "y"=>50.0, "z"=>Dict("start"=>0.0, "step"=>5.0, "end"=>500.0)),
        Dict("x"=>250.0, "y"=>50.0, "z"=>Dict("start"=>0.0, "step"=>5.0, "end"=>500.0)),
        # ...
    ],

    # Snapshot planes — one entry per centre point; each centre produces
    # an XY-, XZ-, and YZ-plane snapshot
    snapshot_positions = [
        Dict("x"=>500.0, "y"=>125.0, "z"=>250.0),
        Dict("x"=>250.0, "y"=>250.0, "z"=>250.0),
        # ...
    ],
    snapshot_times  = [0.4, 0.8],
    snapshot_fields = ["vx", "vy", "vz"],
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
        anisotropic: false

boundaries:
    xstart: absorbing
    xend:   absorbing
    zstart: else
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

```julia
runsim("config2d.yaml", "velmod2d.jld2")
```

---

## Stability and discretisation guidelines

The solver automatically checks the **CFL condition** and reduces `dt` if necessary.
For spatial discretisation, a common rule of thumb is:

$$\Delta x \leq \frac{V_\mathrm{min}}{10 \, f_\mathrm{max}}$$

where $V_\mathrm{min}$ is the minimum phase velocity in the model and $f_\mathrm{max}$ is the maximum
frequency content of the wavelet (typically $\approx 2 \, f_\mathrm{dom}$ for a Ricker wavelet).

The PML thickness should be at least 10 grid cells; thicker layers improve absorption.
For strongly anisotropic media, PML instabilities may occur for certain parameter combinations —
in such cases, increasing the PML thickness often mitigates the issue (see [Method](../method.md)).
