function check_configfile2d(d)

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
        "source.force_angle",
        "source.point_source",
        "source.double_couple",
        "source.location.x",
        "source.location.y",

        "boundaries.left",
        "boundaries.right",
        "boundaries.top",
        "boundaries.bottom",

        "pml.nlayer",
        "pml.reflection_coefficient",

        "receivers.geophones",
        "receivers.das.x_aligned",
        "receivers.das.y_aligned",
        "receivers.snapshots.times",
        "receivers.snapshots.fields"
    ]

    error_name = "Configuration file tests failed."

    for key in required_keys
        keys_hierarchy = split(key, ".")
        current_dict = d
        for k in keys_hierarchy
            if haskey(current_dict, k)
                current_dict = current_dict[k]
            else
                throw(ArgumentError("$error_name '$key' not found"))
            end
        end
    end
end
