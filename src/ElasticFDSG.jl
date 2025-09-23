module ElasticFDSG

module dim2 
using CUDA, Metal, AMDGPU, oneAPI
    using KernelAbstractions, GPUArrays
    using JLD2, HDF5, NPZ, YAML
    using StaticArrays, Einsum, ProgressMeter, UnPack
    using LinearAlgebra, Test, Printf

    include(joinpath(@__DIR__, "2d/main.jl"))
    export runsim
end
export dim2

module dim3
    using CUDA, Metal, AMDGPU, oneAPI
    using KernelAbstractions, GPUArrays
    using JLD2, HDF5, NPZ, YAML
    using StaticArrays, Einsum, ProgressMeter, UnPack
    using LinearAlgebra, Test, Printf

    include(joinpath(@__DIR__, "3d/main.jl"))
    export runsim
end
export dim3

# Utils 
using HDF5
include(joinpath(@__DIR__, "shared/io.jl"))
export load_results, print_h5_tree, config_template

function __init__()
@info "Attempted to import GPU backends (CUDA, Metal, AMDGPU, oneAPI). Any failures or warnings are expected and can be ignored if a backend is not available on your system."
end
end
