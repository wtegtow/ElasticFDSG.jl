using ElasticFDSG
using Test

using HDF5, YAML, JLD2, Statistics

@testset "ElasticFDSG.jl" begin

    # 2D tests
    println("Start 2D tests")
    include(joinpath(@__DIR__,"2d/hom_analytic/run.jl"))

    # 3D tests
    println("Start 3D tests")
    include(joinpath(@__DIR__,"3d/hom_analytic/run.jl"))
     
end
