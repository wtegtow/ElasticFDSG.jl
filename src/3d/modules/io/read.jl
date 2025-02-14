function check_configfile3d(cf)

    """
    Checks that all required informations are present in the configuration file.
    This function does not perform any type checking. 
    Incorrect types may cause errors later in the application.
    """

    required_keys = [
        "settings.device",
        "settings.spatial_derivative_order",
        "settings.precision",
        
        "settings.save_results",
        "settings.show_summary_in_console",
        "settings.show_progress_in_console",

        "settings.output.destination_folder",
        "settings.output.file_name",
        "settings.output.format",

        "time.start",
        "time.end",
        "time.timestep",

        "source.dominant_frequency",
        "source.wavelet_type",
        "source.wavelet_center",
        "source.amplitude",

        "source.point_source.use",
        "source.point_source.act_on.on_vx",
        "source.point_source.act_on.on_vy",
        "source.point_source.act_on.on_vz",
        "source.point_source.act_on.on_all.phi",
        "source.point_source.act_on.on_all.theta",

        "source.double_couple.use",
        "source.double_couple.dip",
        "source.double_couple.strike",
        "source.double_couple.rake",
        "source.double_couple.anisotropic_moment_tensor",

        "source.location.x",
        "source.location.y",
        "source.location.z",

        "boundaries.xstart",
        "boundaries.xend",
        "boundaries.ystart",
        "boundaries.yend",
        "boundaries.zstart",
        "boundaries.zend",

        "pml.nlayer",
        "pml.reflection_coefficient",

        "receivers.geophones",
        "receivers.das.x_aligned",
        "receivers.das.y_aligned",
        "receivers.das.z_aligned",
        "receivers.snapshots.times",
        "receivers.snapshots.fields",
        "receivers.snapshots.origins"
    ]

    error_name = "Configuration file tests failed."

    function key_exists(d, key_path)
        # Base case: If we have reached the last key in the path
        if length(key_path) == 1
            return haskey(d, key_path[1])
        end
        
        # Recursive case: Go deeper into the dictionary
        if haskey(d, key_path[1]) && isa(d[key_path[1]], Dict)
            return key_exists(d[key_path[1]], key_path[2:end])
        else
            return false
        end
    end

    for key in required_keys
        key_path = split(key, ".")
        if !key_exists(cf, key_path)
            throw(ArgumentError("$error_name '$key' not found"))
        end
    end
end