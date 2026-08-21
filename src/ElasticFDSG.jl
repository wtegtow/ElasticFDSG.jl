module ElasticFDSG

    export runsim
    export load_results
    export config_template_2d, config_template_3d

    using YAML, HDF5, NPZ, JLD2
    using LinearAlgebra, Printf, Dates
    using KernelAbstractions, GPUArrays
    using ProgressMeter

    include(joinpath(@__DIR__, "utils.jl"))
    include(joinpath(@__DIR__, "parser.jl"))
    include(joinpath(@__DIR__, "templates.jl"))
    include(joinpath(@__DIR__, "domain.jl"))
    include(joinpath(@__DIR__, "elastic.jl"))
    include(joinpath(@__DIR__, "fields.jl"))
    include(joinpath(@__DIR__, "time.jl"))
    include(joinpath(@__DIR__, "source.jl"))
    include(joinpath(@__DIR__, "cpml.jl"))
    include(joinpath(@__DIR__, "receiver.jl"))
   
    mutable struct FDSG
        config::Config
        device::Device
        domain::Domain
        elastic::Elastic
        fields::Union{Fields2D, Fields3D}
        time::SimTime
        source::Union{Source2D, Source3D}
        pml::Union{CPML2D, CPML3D}
        geophones::Geophones
        das::DAS
        snapshots::Union{Snapshots2D, Snapshots3D}
    end

    include(joinpath(@__DIR__, "summary.jl"))
    include(joinpath(@__DIR__, "solver.jl"))
    include(joinpath(@__DIR__, "io.jl"))

    """
        runsim(config, velmod) -> FDSG or nothing

    Run an elastic forward simulation.

    # Arguments
    - `config`: simulation configuration — either a `Dict` (built with
      [`config_template_2d`](@ref) / [`config_template_3d`](@ref)) or a
      file path (`String`) to a `.yaml` / `.yml` file.
    - `velmod`: velocity model — either a Julia `AbstractArray` or a file path
      (`String`) to a `.jld2`, `.npy`, or `.npz` file.

    # Returns
    - The populated `FDSG` struct when `config["settings"]["output_file"]` is `nothing`.
    - `nothing` when an output file path is given (results are written to HDF5).

    # Example
    ```julia
    fdsg = runsim(config_dict, velmod_array)
    runsum(config, velmod)                  # from dict and array 
    runsim("config.yaml", "velmod.jld2")    # from file paths
    ```
    """
    function runsim(
        config::Union{String, Dict},
        velmod::Union{String, AbstractArray};
        log_level::Symbol=:warn,
        solve::Bool=true 
    )
        set_log_level!(log_level)
        @logger :debug "Hello from ElasticFDSG"

        config = parse_config(config)
        device = parse_device(config)
        velmod = parse_velmod(velmod)

        domain = init_domain(config, velmod)
        elastic = init_elastic(config, domain, velmod)
        velmod = nothing # free memory 
        GC.gc()
        fields = init_fields(config, domain)
        time = init_time(config, domain, elastic)
        source = init_source(config, domain, elastic, time)
        pml = init_cpml(config, domain, elastic, time, source)
        geophones, das, snapshots = init_receiver(config, domain, elastic, time)

        fdsg = FDSG(config, device, domain, elastic, fields, time, source, pml, geophones, das, snapshots)
        if get(config.dict["settings"], "verbose", true)
            _print_summary(fdsg)
        end

        if solve 
            solve!(fdsg) 
        end

        if isnothing(get(fdsg.config.dict["settings"], "output_file", nothing)) 
            return fdsg # return the FDSG struct if no output file specified
        else
            save_results(fdsg)
            return  
        end
    end;
end