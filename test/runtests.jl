using ElasticFDSG
using Test

@testset "ElasticFDSG.jl" begin
    include(joinpath(@__DIR__, "test_dc.jl"))
end;