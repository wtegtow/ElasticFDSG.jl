module ElasticFDSG

    module dim2
        # ----- external dependencies -----
        using Metal, CUDA
        using JLD2, HDF5, NPZ, YAML
        using ProgressMeter

        using Printf, LinearAlgebra
        # ----- scripts -----
        include(joinpath(@__DIR__,"2d/main.jl"))
        export runsim, configtemplate
    end 

    module dim3
        # ----- external dependencies ----- 
        using Metal, CUDA    
        using JLD2, HDF5, NPZ, YAML     
        using ProgressMeter             

        using Printf, LinearAlgebra
        # ----- scripts -----
        include(joinpath(@__DIR__,"3d/main.jl"))
        export runsim, configtemplate
    end 

export dim2, dim3

end
