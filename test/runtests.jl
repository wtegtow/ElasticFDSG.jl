using ElasticFDSG
using Test

@testset "ElasticFDSG.jl" begin
    
    # runs two small demo simulations. If no error occur, the tests are considered passed.
    include(joinpath(@__DIR__,"2d/demo_small/main.jl"))
    include(joinpath(@__DIR__,"3d/demo_small/main.jl"))
    
    # some more in-depth tests might be included in the future. 
end