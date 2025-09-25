module ElasticFDSG

    using HDF5
    include(joinpath(@__DIR__, "shared/io.jl"))
    export load_results, print_h5_tree, config_template

    module dim2 
        using KernelAbstractions, GPUArrays
        using JLD2, HDF5, NPZ, YAML
        using StaticArrays, Einsum, ProgressMeter, UnPack
        using LinearAlgebra, Test, Printf

        include(joinpath(@__DIR__, "2d/main.jl"))
        export runsim
    end
    export dim2

    module dim3
        using KernelAbstractions, GPUArrays
        using JLD2, HDF5, NPZ, YAML
        using StaticArrays, Einsum, ProgressMeter, UnPack
        using LinearAlgebra, Test, Printf

        include(joinpath(@__DIR__, "3d/main.jl"))
        export runsim
    end
    export dim3

end