# ElasticFDSG.jl

<img src="docs/src/assets/logo.png" alt="ElasticFDSG Logo" width="200"/>

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://wtegtow.github.io/ElasticFDSG.jl/dev/)
[![Build Status](https://github.com/wtegtow/ElasticFDSG.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/wtegtow/ElasticFDSG.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![DOI](https://zenodo.org/badge/923201339.svg)](https://doi.org/10.5281/zenodo.14872584)

ElasticFDSG.jl is a Julia package for solving the elastic wave equation in the velocity–stress formulation using the finite-difference method on a staggered grid.

But why another implementation? In many existing tools, installing dependencies, configuring simulation setups, defining model parameters, and accessing results can be already complicated tasks — particularly for inexperienced users who are looking for a quick and straightforward workflow.
ElasticFDSG was developed to offer a user-friendly experience while also maintaining flexibility to be applied to a wide variety of simulation scenarios.
Users can easily customize their simulations by creating velocity models and configuration files in a straightforward manner, that can be directly passed to the solvers.


## Features 

- 2D and 3D elastic forward modelling on regular grids.
- Vendor neutral CPU and GPU kernel (CPU, CUDA, Metal, AMDGPU, oneAPI) using [KernelAbstractions.jl](https://github.com/JuliaGPU/KernelAbstractions.jl).
- Spatial derivatives of order 1 to 10.
- Second order time marching.
- Elastic isotropic or vertical transversal isotrop (VTI) 2D models using 2D Thomson parameter.
- Elastic isotropic, VTI or orthorhombic (ORT) 3D models using Tsvankin parameter.
- Elastic properties defined on full integer grid points.
- Solver can handle fully heterogeneous media.
- Absorbing boundaries using Convolutional-Perfectly-Matched-Layer.
- Moment tensor sources. 
- Save geophone receiver (velocity point sensors). 
- Save Distributed Acoustic Sensing (DAS) receiver, aligned with model coordinate axis (strain-profiles).
- Save snapshots at specified time steps.
- Easy-to-read source code.

A step by step user guide can be found in the [documentation](https://wtegtow.github.io/ElasticFDSG.jl/dev/).

Basic examples are included in the git-repository `examples/` folder. 

![Demo](docs/src/assets/readme_animation.png)

## Installation

```julia-repl
using Pkg
Pkg.add(ElasticFDSG)
```

## Citing
If you find this package helpful for your research, please consider citing:

```
@misc{ElasticFDSG,
  author       = {William Tegtow},
  title        = {ElasticFDSG.jl: Simulating elastic wave propagation in 2D and 3D anisotropic media.},
  year         = {2025},
  doi          = {https://doi.org/10.5281/zenodo.14872584}
}

```

**Note:**
This package is still in its early stages, and only limited testing has been done so far. Any bug report or suggestion is very welcomed.

**Note:**
Old scripts may need to be revised to be compatible with the latest version.
- With v1.0.2, the structure of the configuration.yaml files has slightly changed. 
- With v1.0.2, 2D velocity models must now be saved as (Nx,Nz) arrays instead of previously (Nz,Nx). 