# This test script runs a small 3d demo simulation. 

function test_small_3d_demo()
    config_path = joinpath(@__DIR__, "config.yaml")
    velmod_path = joinpath(@__DIR__, "velmod.jld2")

    no_error_occured = true 
    try
        ElasticFDSG.dim3.runsim(config_path, velmod_path)
    catch 
        no_error_occured = false 
    end 

    @test no_error_occured
end


test_small_3d_demo()