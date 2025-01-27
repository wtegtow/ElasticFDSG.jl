mutable struct FDSG3D{a,b,c,d,e,f,g,h,i,j}
    settings::a
    domain::b
    elastic::c
    fields::d 
    time::e 
    source::f 
    pml::g
    geophones::h
    das::i 
    snapshots::j
end

include(joinpath(@__DIR__, "../shared/__export__.jl"))
include(joinpath(@__DIR__, "modules/__export__.jl"))

"""
    Elastic_FDSG.dim3.runsim(CONFIGPATH, VELMODPATH)

Run the 3D elastic forward simulation using the specified configuration and velocity model.

# Arguments
- `CONFIGPATH::String`: Path to the `configuration.yaml` file containing simulation settings.
- `VELMODPATH::String`: Path to the velocity model file. Supported formats include `.jld2`, `.npy`, and `.npz`.

# Returns
- `Nothing`: The function runs the simulation and saves the results as specified in the configuration file.
"""
function runsim(CONFIGPATH, VELMODPATH)
    # SETUP
    mlog = MessageLog(String[]);
    settings = init_settings(CONFIGPATH, mlog)
    domain   = init_domain(VELMODPATH, settings)
    elastic  = init_elastic(settings, domain)
    fields   = init_fields(settings, domain)
    time     = init_time(settings, elastic, domain, mlog)
    source   = init_source(settings, domain, elastic, time, mlog)
    pml      = init_pml(settings, domain, elastic, time, source)
    # RECEIVER 
    geophones = init_geophones(settings, domain, time)
    das       = init_das(settings, domain, elastic, time)
    snapshots = init_snapshots(settings, domain, time)
    # INIT
    fdsg3d = FDSG3D{Settings, Domain, Elastic, Fields, Time, Source, Pml, Geophones, DAS, Snapshots}(
                    settings, domain, elastic, fields, time, source, pml, geophones, das, snapshots)
    iwindow(mlog, fdsg3d)
    # SOLVE 
    solver = init_solver(fdsg3d)
    solver(fdsg3d)
    # SAVE RESULTS
    save_results(fdsg3d)
end;