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
ElasticFDSG.dim3.runsim(CONFIGPATH, VELMODPATH)
Run the 3D elastic forward simulation using the specified configuration and velocity model.
# Arguments:
- `CONFIGPATH::String`: Path to the `configuration.yaml` file containing simulation settings.
- `VELMODPATH::String`: Path to the velocity model file. Supported formats: `.jld2`, `.npy`, and `.npz`.
Optional Keyword Arguments:
- `block_size::Tuple=(7,7,7)` : Specify numbers of threads launched on GPU. For CPU it will have no effect.
- `return_::Bool=false`: False -> Results will be saved as .h5 file as specified in the config file. 
                         True  -> Function returns FDSG3D struct.
"""
function runsim(CONFIGPATH::String, VELMODPATH::String; block_size=(7,7,7), return_ = false)
    settings = init_settings(CONFIGPATH; dim=3);
    if settings.showinfo log_progress(" >Log setup...") end
    velmod   = load_velmod(VELMODPATH);
    if settings.showinfo log_progress(" >velocity model initialized") end
    domain   = init_domain(settings, velmod);
    if settings.showinfo log_progress(" >domain initialized") end
    elastic  = init_elastic(settings, domain, velmod);
    if settings.showinfo log_progress(" >elastic parameters initialized") end
    fields   = init_fields(settings, domain);
    if settings.showinfo log_progress(" >fields initialized") end
    velmod = nothing; # free memory 
    time     = init_time(settings, elastic, domain);
    if settings.showinfo log_progress(" >time initialized") end
    source   = init_source(settings, domain, elastic, time);
    if settings.showinfo log_progress(" >source initialized") end
    pml      = init_pml(settings, domain, elastic, time, source);
    if settings.showinfo log_progress(" >pml initialized") end
    geophones, das, snapshots = init_receiver(settings, domain, elastic, time)
    if settings.showinfo log_progress(" >receivers initialized") end

    fdsg3d = FDSG3D{Settings, Domain, Elastic, Fields, Time, Source, Pml, Geophones, DAS, Snapshots}(
                    settings, domain, elastic, fields, time, source, pml, geophones, das, snapshots);
    if settings.showinfo log_progress(" >FDSG3D struct initialized") end
    if settings.showinfo log_progress("") end
    iwindow(fdsg3d)
    solve!(fdsg3d; block_size=block_size)
    if !return_
         save_results(fdsg3d)
         return nothing 
    else return fdsg3d
    end;
end;