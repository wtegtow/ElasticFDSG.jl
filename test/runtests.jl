using ElasticFDSG
using Test

using HDF5, YAML, Statistics

@testset "ElasticFDSG.jl" begin
    # 2D tests
    include(joinpath(@__DIR__,"2d/hom_analytic/run.jl"))

    # 3D tests
     
end
