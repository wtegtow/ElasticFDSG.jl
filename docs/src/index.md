# ElasticFDSG.jl

**ElasticFDSG.jl** is a Julia package for simulating elastic wave propagation in 2D and 3D heterogeneous anisotropic media.
It solves the elastic wave equation in the velocity–stress formulation using a finite-difference staggered-grid (FDSG) scheme.

The package was developed with a focus on a clean, user-friendly workflow: simulations are fully described by a velocity model and a configuration dictionary (or YAML file), and results are returned as a Julia struct or saved to HDF5.

## Features

- 2D and 3D elastic forward modelling on regular grids.
- Vendor-neutral CPU and GPU kernels (CUDA, Metal, AMDGPU, oneAPI) via [KernelAbstractions.jl](https://github.com/JuliaGPU/KernelAbstractions.jl).
- Spatial finite-difference operators of order 1-10.
- Second-order leapfrog time integration.
- Heterogeneous isotropic and VTI 2D models (Thomsen parameters).
- Heterogeneous isotropic, VTI, and orthorhombic 3D models (Tsvankin parameters).
- Convolutional Perfectly Matched Layer (C-PML) absorbing boundaries.
- Moment tensor sources.
- Geophone receivers (point particle velocity).
- DAS receivers (axial strain along coordinate-aligned profiles).
- Wavefield snapshots.
- Results saved to HDF5 or returned as a Julia struct.

A step-by-step user guide can be found in the [User Guide](userguide/intro.md).
Working examples are available in the [`examples/`](https://github.com/wtegtow/ElasticFDSG.jl/tree/main/examples) directory of the repository.

## Installation

```julia-repl
julia> using Pkg
julia> Pkg.add(url="https://github.com/wtegtow/ElasticFDSG.jl")
```

## Quick start

```julia
using ElasticFDSG, GLMakie

# Build a minimal 2D velocity model (7 × nx × nz)
h = 10 
x = 0:h:2000
z = 0:h:2000
X  = repeat(xc,  1, nz)
Z  = repeat(reshape(zc, 1, :), nx, 1)

# Set velocities 
vp0 = 4400 
vs0 = 2200
rho0 = 2000 

velmod = zeros(7, nx, nz)
velmod[1,:,:] .= X
velmod[2,:,:] .= Z
velmod[3,:,:] .= vp0
velmod[4,:,:] .= vs0    
velmod[5,:,:] .= rho0                      
# indices 6 & 7 (Thomsen ε, δ) left at zero → isotropic

# Build a configuration dictionary
config = config_template_2d(
    device       = "cpu",
    precision    = "Float64",
    fd_order     = 4,
    verbose      = true,
    output_file  = nothing,   # return struct instead of saving
    t_start      = 0.0,
    t_end        = 0.4,
    dt           = 0.001,
    fdom         = 25.0,
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
    geophones    = [Dict("x"=>1500.0,"z"=>500.0)
                    # ... add more here
    ],
    das_x_aligned = [], 
    das_z_aligned = [
        Dict("x" => 500, "z" => Dict("start"=>500, "step"=>5, "end"=>15000)),
        # ... add more here 
    ],
    snapshot_times  = [0.25, 0.3],
    snapshot_fields = ["vx", "vz"],
)

# Run simulation — dimension is auto-detected from the velmod array
fdsg = runsim(config, velmod); 

# Unpack from fdsg struct 
geophones = fdsg.geophones;
time = fdsg.time.t
geo_data = fdsg.geophones.data 

# Visualize
fig = Figure(size=(800,300))
ax1 = Axis(fig[1,1], title="vx"); ax2 = Axis(fig[1,2], title="vz")
lines!(ax1, time, geo_data[1,1,:], color="black")
lines!(ax2, time, geo_data[1,2,:], color="black")
display(fig)
```


## Citing

If you find ElasticFDSG.jl usefull for your research, consider citing:

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
