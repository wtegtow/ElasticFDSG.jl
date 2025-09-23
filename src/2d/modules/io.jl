function load_velmod(VELMODPATH)
    if endswith(VELMODPATH,"npz")
        velmodfile = npzread(joinpath(VELMODPATH))
        velmod = velmodfile["velmod"];
    elseif endswith(VELMODPATH,"npy")
        velmod = npzread(joinpath(VELMODPATH))
    elseif endswith(VELMODPATH,"jld2")
        velmodfile = jldopen(joinpath(VELMODPATH))
        velmod = velmodfile["velmod"];
        close(velmodfile)
    else
        error("No valid velocity model file found.")
    end;
    return velmod 
end 

function check_configfile2d(config)

    function test_info_block(info_block::Dict, required_info_block)
        for (key, expected) in required_info_block
            @assert haskey(info_block, key) "Configuration-file check failed → missing key: $(key) in block: $(keys(info_block))"
            value = info_block[key]

            # has subblocks -> recursion
            if isa(expected, Vector) && all(x -> x isa Tuple, expected)
                @assert value isa Dict "Configuration-file check failed → key $(key) must be a sub-block (Dict), got $(typeof(value))"
                test_info_block(value, expected)

            elseif isa(expected, Vector)
                @assert value in expected "Configuration-file check failed → value $(value) not in $(expected)"

            else
                @assert value isa expected "Configuration-file check failed → key $(key) must be $(expected), got $(typeof(value))"
            end
        end
    end;

    required_settings = [
        ("device", String),
        ("precision", ["Float32", "Float64"]),
        ("spatial_derivative_order", Int),
        ("show_progress_in_console", Bool),
        ("output_file", String)
    ];

    required_time = [
        ("start", Real),
        ("end", Real),
        ("timestep", Real)
    ]

    required_source = [
        ("dominant_frequency", Real),
        ("wavelet_type", ["ricker", "gauss1d"]),
        ("wavelet_center", Real),
        ("location", [
            ("x", Real),
            ("z", Real)
        ]),
        ("seismic_moment", Real),
        ("moment_tensor", [
            ("Mxx", Real),
            ("Mzz", Real),
            ("Mxz", Real),
            ("anisotropic", Bool)
        ])
    ]

    required_boundaries = [
        ("xstart", String),
        ("xend",   String),
        ("zstart", String),
        ("zend",   String),
        ("pml_layer", Int)
    ]

    required_receivers = [
        ("geophones", Any), # kann leer sein
        ("das", [ 
            ("x_aligned", Any), 
            ("z_aligned", Any)
        ]),
        ("snapshots", [
            ("times", Any),
            ("fields", Any)
        ])
    ]

    test_info_block(config["settings"], required_settings)
    test_info_block(config["time"], required_time)
    test_info_block(config["source"], required_source)
    test_info_block(config["boundaries"], required_boundaries)
    test_info_block(config["receivers"], required_receivers)

end;


function save_results(fdsg2d::FDSG2D)

    file = h5open(fdsg2d.settings.filename, "w")

    g1 = create_group(file, "grid")
    g1["x_coordinates"] = collect(fdsg2d.domain.xcoords[fdsg2d.domain.inner_id[2]])
    g1["z_coordinates"] = collect(fdsg2d.domain.zcoords[fdsg2d.domain.inner_id[1]])

    g2 = create_group(file, "time")
    g2["time"] = collect(fdsg2d.time.t)
    g2["t0"] = fdsg2d.time.t0
    g2["tend"] = fdsg2d.time.tend
    g2["dt"] = fdsg2d.time.dt

    g3 = create_group(file, "source")
    g3["dominant_frequency"] = fdsg2d.source.fdom 
    g3["source_time_function"] = collect(fdsg2d.source.stf)
    g3["Mxx"] = fdsg2d.source.Mxx
    g3["Mxz"] = fdsg2d.source.Mxz
    g3["Mzz"] = fdsg2d.source.Mzz 
    g3["location_xz"] = [fdsg2d.source.x, fdsg2d.source.z]

    if fdsg2d.geophones.n > 0
        g4 = create_group(file, "geophones")
        for geo_idx in 1:fdsg2d.geophones.n
            g4_ = create_group(g4, "geophone_$(geo_idx)")
            g4_["data"] = fdsg2d.geophones.data[geo_idx,:,:]
            g4_["location_xz"] = fdsg2d.geophones.coords[:,geo_idx]
        end 
    end

    # das 
    if fdsg2d.das.x_aligned.n > 0  ||
       fdsg2d.das.z_aligned.n > 0
        g5 = create_group(file, "das")
    end;

    if fdsg2d.das.x_aligned.n > 0
        g5x = create_group(g5, "x_aligned")
        for i in 1:fdsg2d.das.x_aligned.n
        g5x_ = create_group(g5x, "fiber_$(i)")
        g5x_["data"] = fdsg2d.das.x_aligned.data[i]
        v = fdsg2d.das.x_aligned.coords[i]
        g5x_["location_xz"] = hcat(v[1], fill(v[2], length(v[1])))
        end
    end;

    if fdsg2d.das.z_aligned.n > 0
        g5z = create_group(g5, "z_aligned")
        for i in 1:fdsg2d.das.z_aligned.n
        g5z_ = create_group(g5z, "fiber_$(i)")
        g5z_["data"] = fdsg2d.das.z_aligned.data[i]
        v = fdsg2d.das.z_aligned.coords[i]
        g5z_["location_xz"] = hcat(fill(v[1], length(v[2])), v[2])
        end
    end 

    if fdsg2d.snapshots.n > 0
        g6 = create_group(file, "snapshots")
        g6["fields"] = fdsg2d.snapshots.fieldnames
        g6["times"] = fdsg2d.settings.float.(fdsg2d.snapshots.t)
        g6["XZ"] = permutedims(fdsg2d.snapshots.XZ[:,:,fdsg2d.domain.inner_id[1], fdsg2d.domain.inner_id[2]], (1,2,4,3))
    end 
    println("Results saved.")
    close(file)    
end

