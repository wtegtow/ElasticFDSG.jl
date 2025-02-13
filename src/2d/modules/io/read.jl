
function check_configfile2d(cf)
   
    required_keys = [
    "settings.device",
    "settings.precision",
    "settings.spatial_derivative_order",
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
    "source.location.x",
    "source.location.y",
    "source.point_source.use",
    "source.point_source.act_on.on_vx",
    "source.point_source.act_on.on_vy",
    "source.point_source.act_on.on_vx_and_vy.force_angle",
    "source.double_couple.use",
    "source.double_couple.strike",
    "boundaries.xstart",
    "boundaries.xend",
    "boundaries.ystart",
    "boundaries.yend",
    "pml.nlayer",
    "pml.reflection_coefficient",
    "receivers.geophones",
    "receivers.das.x_aligned",
    "receivers.das.y_aligned",
    "receivers.snapshots.times",
    "receivers.snapshots.fields"
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
end;