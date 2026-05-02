# ElasticFDSG.jl

**ElasticFDSG.jl** is a Julia package for simulating elastic wave propagation in 2D and 3D heterogeneous anisotropic media.
It solves the elastic wave equation in the velocity–stress formulation using a finite-difference staggered-grid (FDSG) scheme.

The package was developed with a focus on a clean, user-friendly workflow: simulations are fully described by a velocity model and a configuration dictionary (or YAML file), and results are returned as a Julia struct or saved to HDF5.

## Features

- 2D and 3D elastic forward modelling on regular grids.
- Vendor-neutral CPU and GPU kernels (CUDA, Metal, AMDGPU, oneAPI) via [KernelAbstractions.jl](https://github.com/JuliaGPU/KernelAbstractions.jl).
- Spatial finite-difference operators of order 2 to 20 (stencil half-width 1–10).
- Second-order leapfrog time integration.
- Elastic isotropic and VTI 2D models parameterised via Thomsen parameters ($\epsilon$, $\delta$).
- Elastic isotropic, VTI, and orthorhombic 3D models parameterised via Tsvankin parameters ($\epsilon_{1,2}$, $\gamma_{1,2}$, $\delta_{1,2,3}$).
- Fully heterogeneous media — stiffness tensors are compressed using a unique-tensor lookup.
- Convolutional Perfectly Matched Layer (C-PML) absorbing boundaries.
- Moment tensor sources with Ricker or Gaussian derivative wavelets.
- Geophone receivers (point particle velocity).
- DAS receivers (axial strain along coordinate-aligned profiles).
- Wavefield snapshots at arbitrary time steps.
- Results serialised to HDF5 or returned directly as a Julia struct.

A step-by-step user guide can be found in the [User Guide](userguide/intro.md).
Working examples are available in the [`examples/`](https://github.com/wtegtow/ElasticFDSG.jl/tree/main/examples) directory of the repository.

## Installation

```julia-repl
julia> using Pkg
julia> Pkg.add(url="https://github.com/wtegtow/ElasticFDSG.jl")
```

## Quick start

```julia
using ElasticFDSG

# Build a minimal 2D velocity model (7 × nx × nz)
nx, nz = 200, 200
h = 10.0   # grid spacing [m]
xc = range(0.0, step=h, length=nx)
zc = range(0.0, step=h, length=nz)
X  = repeat(xc,  1, nz)
Z  = repeat(reshape(zc, 1, :), nx, 1)

velmod = zeros(7, nx, nz)
velmod[1,:,:] .= X;    velmod[2,:,:] .= Z
velmod[3,:,:] .= 3000; velmod[4,:,:] .= 1800   # vp, vs [m/s]
velmod[5,:,:] .= 2500                           # density [kg/m³]
# layers 6 & 7 (Thomsen ε, δ) left at zero → isotropic

# Build a configuration dictionary
config = config_template_2d(
    device       = "cpu",
    precision    = "Float32",
    fd_order     = 4,
    verbose      = true,
    output_file  = nothing,           # return struct instead of saving
    t_start      = 0.0,
    t_end        = 0.5,
    dt           = 0.001,
    fdom         = 30.0,
    wavelet      = "ricker",
    wavelet_center = 0.05,
    seismic_moment = 1e6,
    src_x        = 1000.0,
    src_z        = 500.0,
    Mxx = 0.0, Mxz = 1.0, Mzz = 0.0,
    anisotropic  = false,
    xstart = "absorbing", xend = "absorbing",
    zstart = "absorbing", zend = "absorbing",
    pml_layer    = 10,
    geophones    = [Dict("x"=>1500.0,"z"=>500.0)],
    das_x_aligned = [],
    das_z_aligned = [],
    snapshot_times  = [0.25, 0.5],
    snapshot_fields = ["vx", "vz"],
)

# Run simulation — dimension is auto-detected from the velmod array
fdsg = runsim(config, velmod)
```

## Citing

If you use ElasticFDSG.jl in your research, please cite:

```bibtex
@misc{ElasticFDSG,
  author = {William Tegtow},
  title  = {ElasticFDSG.jl: Simulating elastic wave propagation in 2D and 3D anisotropic media.},
  year   = {2025},
  doi    = {https://doi.org/10.5281/zenodo.14872584}
}
```

!!! note
    This package is under active development. Bug reports and suggestions are very welcome.
