using ElasticFDSG
using Test

@testset "ElasticFDSG.jl" begin
    include(joinpath(@__DIR__, "validation.jl"))
end;