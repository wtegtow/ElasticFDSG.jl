struct Settings{T, U}
    config::Dict{Any, Any}
    save::Bool 
    showinfo::Bool
    float::DataType 
    device::String 
    device_name::U 

    filename::String

    N::Int64 
    c::T
end



function init_settings(CONFIGPATH, mlog)

    config = YAML.load_file(CONFIGPATH)
    check_configfile2d(config) # check if the config file is complete

    # some general settings
    float = eval(Symbol(config["settings"]["precision"]))
    device = config["settings"]["device"]
    device_name = check_device(device)
    showinfo = config["settings"]["show_summary_in_console"]
    N = Int64(config["settings"]["spatial_derivative_order"])
    c = float.(diff_coeff(N))

    # output related 
    save = config["settings"]["save_results"]

    output_dir = config["settings"]["output"]["destination_folder"]
    if save == true check_dir(output_dir) end # check if output_dir exists

    format = config["settings"]["output"]["format"]
    if save == true check_format(format) end

    filename = config["settings"]["output"]["file_name"]
    output_path = joinpath(output_dir, (filename*"."*format))
    if save == true check_filename(output_path, mlog) end
    
    # types
    T = typeof(c)
    U = typeof(device_name)

    settings = Settings{T, U}(config, save, showinfo, float, device, device_name, output_path, N, c)
    
    return settings

end;