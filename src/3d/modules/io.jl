function check_configfile3d(config)

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
            ("y", Real),
            ("z", Real)
        ]),
        ("seismic_moment", Real),
        ("moment_tensor", [
            ("Mxx", Real),
            ("Myy", Real),
            ("Mzz", Real),
            ("Mxy", Real),
            ("Mxz", Real),
            ("Myz", Real),
            ("anisotropic", Bool)
        ])
    ]

    required_boundaries = [
        ("xstart", String),
        ("xend",   String),
        ("ystart", String),
        ("yend",   String),
        ("zstart", String),
        ("zend",   String),
        ("pml_layer", Int)
    ]

    required_receivers = [
        ("geophones", Any), # kann leer sein
        ("das", [ 
            ("x_aligned", Any), 
            ("y_aligned", Any), 
            ("z_aligned", Any)
        ]),
        ("snapshots", [
            ("plane_positions", Any),
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


function save_results(fdsg3d::FDSG3D)

    file = h5open(fdsg3d.settings.filename, "w")

    g1 = create_group(file, "grid")
    g1["x_coordinates"] = collect(fdsg3d.domain.xcoords[fdsg3d.domain.inner_id[3]])
    g1["y_coordinates"] = collect(fdsg3d.domain.ycoords[fdsg3d.domain.inner_id[2]])
    g1["z_coordinates"] = collect(fdsg3d.domain.zcoords[fdsg3d.domain.inner_id[1]])

    g2 = create_group(file, "time")
    g2["time"] = collect(fdsg3d.time.t)
    g2["t0"] = fdsg3d.time.t0
    g2["tend"] = fdsg3d.time.tend
    g2["dt"] = fdsg3d.time.dt

    g3 = create_group(file, "source")
    g3["dominant_frequency"] = fdsg3d.source.fdom 
    g3["source_time_function"] = collect(fdsg3d.source.stf)
    g3["Mxx"] = fdsg3d.source.Mxx
    g3["Mxy"] = fdsg3d.source.Mxy
    g3["Mxz"] = fdsg3d.source.Mxz
    g3["Myy"] = fdsg3d.source.Myy
    g3["Myz"] = fdsg3d.source.Myz 
    g3["Mzz"] = fdsg3d.source.Mzz 
    g3["location_xyz"] = [fdsg3d.source.x, fdsg3d.source.y, fdsg3d.source.z]

    if fdsg3d.geophones.n > 0
        g4 = create_group(file, "geophones")
        for geo_idx in 1:fdsg3d.geophones.n
            g4_ = create_group(g4, "geophone_$(geo_idx)")
            g4_["data"] = fdsg3d.geophones.data[geo_idx,:,:]
            g4_["location_xyz"] = fdsg3d.geophones.coords[:,geo_idx]
        end 
    end

    # das 
    if fdsg3d.das.x_aligned.n > 0  ||
        fdsg3d.das.y_aligned.n > 0 ||
        fdsg3d.das.z_aligned.n > 0

        g5 = create_group(file, "das")
    end;

    if fdsg3d.das.x_aligned.n > 0
        g5x = create_group(g5, "x_aligned")
        for i in 1:fdsg3d.das.x_aligned.n
        g5x_ = create_group(g5x, "fiber_$(i)")
        g5x_["data"] = fdsg3d.das.x_aligned.data[i]
        v = fdsg3d.das.x_aligned.coords[i]
        g5x_["location_xyz"] = hcat(v[1], fill(v[2], length(v[1])), fill(v[3], length(v[1])))
        end
    end;

    if fdsg3d.das.y_aligned.n > 0
        g5y = create_group(g5, "y_aligned")
        for i in 1:fdsg3d.das.y_aligned.n
        g5y_ = create_group(g5y, "fiber_$(i)")
        g5y_["data"] = fdsg3d.das.y_aligned.data[i]
        v = fdsg3d.das.y_aligned.coords[i]
        g5y_["location_xyz"] = hcat(fill(v[1], length(v[2])), v[2], fill(v[3], length(v[2])))
        end
    end 

    if fdsg3d.das.z_aligned.n > 0
        g5z = create_group(g5, "z_aligned")
        for i in 1:fdsg3d.das.z_aligned.n
        g5z_ = create_group(g5z, "fiber_$(i)")
        g5z_["data"] = fdsg3d.das.z_aligned.data[i]
        v = fdsg3d.das.z_aligned.coords[i]
        g5z_["location_xyz"] = hcat(fill(v[1], length(v[3])), fill(v[2], length(v[3])), v[3])
        end
    end 

    if fdsg3d.snapshots.n > 0
        g6 = create_group(file, "snapshots")
        g6["fields"] = fdsg3d.snapshots.fieldnames
        g6["times"] = fdsg3d.snapshots.t
        for i in 1:fdsg3d.snapshots.n
            xidx, yidx, zidx = fdsg3d.snapshots.grid_ids[i]
            g6_ = create_group(g6, "data_$(i)")
            g6_["XY_data"] = permutedims(fdsg3d.snapshots.XY[i,:,:,fdsg3d.domain.inner_id[2], fdsg3d.domain.inner_id[3]],(1, 2, 4, 3))
            g6_["XY_plane_Z"] = fdsg3d.domain.zcoords[zidx]
            g6_["XZ_data"] = permutedims(fdsg3d.snapshots.XZ[i,:,:,fdsg3d.domain.inner_id[1],fdsg3d.domain.inner_id[3]],(1, 2, 4, 3))
            g6_["XZ_plane_Y"] = fdsg3d.domain.ycoords[yidx]
            g6_["YZ_data"] = permutedims(fdsg3d.snapshots.YZ[i,:,:,fdsg3d.domain.inner_id[1],fdsg3d.domain.inner_id[2]],(1, 2, 4, 3))
            g6_["YZ_plane_X"] = fdsg3d.domain.xcoords[xidx]
        end 
    end 
    println("Results saved.")
    close(file)    
end
