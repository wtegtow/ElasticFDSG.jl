module ElasticFDSG

    export runsim, load_results

    using LinearAlgebra, Printf, Dates
    using KernelAbstractions, GPUArrays
    using YAML, HDF5, NPZ, JLD2, ProgressMeter

    include(joinpath(@__DIR__, "logger.jl"))
    include(joinpath(@__DIR__, "parser.jl"))
    include(joinpath(@__DIR__, "domain.jl"))
    include(joinpath(@__DIR__, "elastic.jl"))
    include(joinpath(@__DIR__, "time.jl"))
    include(joinpath(@__DIR__, "source.jl"))
    include(joinpath(@__DIR__, "fields.jl"))
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
        runsim(config, velmod) -> nothing

    Run an elastic forward simulation and save results to HDF5 as specified in config.

    # Arguments
    - `config`: simulation configuration — either a `Dict` or a file path to a `.yaml` / `.yml` file.
    - `velmod`: velocity model — either a Julia `AbstractArray` or a file path to a `.jld2`, `.npy`, or `.npz` file.

    # Example
    ```julia
    runsim(config, velmod)                  # from dict and array 
    runsim("config.yaml", "velmod.jld2")    # from file paths
    ```
    """
    function runsim(
        config::Union{String, Dict},
        velmod::Union{String, AbstractArray};
        log_level::Symbol=:warn,
        _return::Bool=false)

        set_log_level!(log_level)
        @logger :info "Hello from ElasticFDSG"

        config = parse_config(config)
        device = parse_device(config)
        velmod = parse_velmod(velmod)

        domain                      = init_domain(config, velmod)
        elastic                     = init_elastic(config, domain, velmod)
        time                        = init_time(config, domain, elastic)
        source                      = init_source(config, domain, elastic, time)
        geophones, das, snapshots   = init_receiver(config, domain, elastic, time)
        velmod = nothing; GC.gc()   # free velmod memory before allocating new fields
        fields                      = init_fields(config, domain)
        pml                         = init_cpml(config, domain, elastic, time, source)

        fdsg = FDSG(config, device, domain, elastic, fields, time, source, pml, geophones, das, snapshots)
        config.dict["settings"]["verbose"] && _print_summary(fdsg)
        
        solve!(fdsg) 
        save_results(fdsg)
        _return && return fdsg 
    end
end