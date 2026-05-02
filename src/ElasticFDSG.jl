module ElasticFDSG

    export devmode!
    export runsim
    export load_results
    export config_template_2d, config_template_3d

    using YAML, HDF5, NPZ, JLD2
    using LinearAlgebra, Printf
    using KernelAbstractions, GPUArrays
    using ProgressMeter

    const _DEV = Ref(false)
    devmode!(on::Bool=true) = (_DEV[] = on)
    _log(msg) = _DEV[] && println("  [dev] ", msg)

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

    The simulation dimension (2D / 3D) is detected automatically:
    a `(7, nx, nz)` array triggers a 2D run; a `(13, nx, ny, nz)` array triggers 3D.

    # Returns
    - The populated `FDSG` struct when `settings.output_file` is `nothing`.
    - `nothing` when an output file path is given (results are written to HDF5).

    # Example
    ```julia
    fdsg = runsim(config_dict, velmod_array)
    runsim("config.yaml", "velmod.jld2")   # saves to HDF5
    ```
    """
    function runsim(
        config::Union{String, Dict},
        velmod::Union{String, AbstractArray},
    )
        
        _log("Hello from ElasticFDSG!")

        config = parse_config(config)
        _log("Config parsed successfully")

        device = parse_device(config)
        _log("Running on: $(device.name)")

        velmod = parse_velmod(velmod)
        _log("Velocity model parsed successfully")

        @assert config.dim == velmod.dim "Dimension mismatch between config and velocity model" 

        domain = init_domain(config, velmod)
        _log("Domain initialized with shape $(domain.shape)")

        elastic = init_elastic(config, domain, velmod)
        _log("Elastic parameters initialized")

        velmod = nothing # free memory 
        GC.gc()

        fields = init_fields(config, domain)
        _log("Fields initialized")

        time = init_time(config, domain, elastic)
        _log("Time parameters initialized")

        source = init_source(config, domain, elastic, time)
        _log("Source initialized")

        pml = init_cpml(config, domain, elastic, time, source)
        _log("PML initialized")

        geophones, das, snapshots = init_receiver(config, domain, elastic, time)
        _log("Receivers initialized")

        fdsg = FDSG(config, device, domain, elastic, fields, time, source, pml, geophones, das, snapshots)
        _log("FDSG struct initialized")

        if get(config.dict["settings"], "verbose", true)
            _print_summary(fdsg)
        end

        solve!(fdsg)
        _log("Simulation complete")
    
        if isnothing(get(fdsg.config.dict["settings"], "output_file", nothing)) 
            return fdsg # return the FDSG struct if no output file specified
        else
            save_results(fdsg)
            return  
        end
    end;
end