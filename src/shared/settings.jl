struct Settings{T1,T2,T3} 
    config::T1
    showinfo::Bool
    float::DataType
    filename::String
    N::Int64
    c::T2
    device::String
    device_name::String
    array::T3
end


function init_settings(CONFIGPATH; dim=3)

    # read and check config file
    config = YAML.load_file(CONFIGPATH)
    if dim == 2 
        check_configfile2d(config) 
    elseif dim == 3 
        check_configfile3d(config) 
    else 
        error("dim != {2,3}")
    end
    # set device and corresponding array types
    device = String(config["settings"]["device"]) 
    name_  = nothing
    array_ = nothing 

    # cpu
    if device in ["cpu", "CPU", "Cpu"]
        device = "cpu"
        name_  = "cpu"
        array_ = Array
    
    # nvidia gpu
    elseif device in ["gpu-cuda", "gpu_cuda", "cuda", "Cuda", "CUDA"]

        if !isdefined(Main, :CUDA)
            error("CUDA.jl not loaded in Main namespace.")
        end
        device = "gpu-cuda"
        device_ids = Main.CUDA.devices()
        if isempty(device_ids)
            error("Cuda device not found.")
        else
            name_ = String(Main.CUDA.name(Main.CUDA.device!(0)))
            array_ = Main.CuArray
        end

    # apple gpu
    elseif device in ["gpu-metal", "gpu_metal", "metal", "apple", "Apple", "MacBook", "Mac"]

        if !isdefined(Main, :Metal)
            error("Metal.jl not loaded in Main namespace.")
        end

        device = "gpu-metal"
        device_list = Main.Metal.MTL.devices()
        if isempty(device_list)
            error("Metal device not found.")
        else 
            name_  = String(device_list[1].name)
            array_ = Main.MtlArray
        end

    # AMD (not testet)
    elseif device in ["gpu-amd", "gpu_rocm", "amd", "AMD"]

        if !isdefined(Main, :AMDGPU)
            error("AMDGPU.jl not loaded in Main namespace.")
        end

        device = "gpu-amd"
        device_list = Main.AMDGPU.devices()
        if isempty(device_list)
            error("AMD device not found.")
        end
        name_ = String(Main.AMDGPU.devices()[1])
        array_ = Main.ROCArray

    # oneAPI (not testet)  
    elseif device in ["gpu-intel", "intel", "oneapi", "one_api"]

        if !isdefined(Main, :oneAPI)
            error("oneAPI.jl not in Main namespace.")
        end

        device = "gpu-intel"
        device_list = Main.oneAPI.device()
        if isempty(device_list)
            error("AMD device not found.")
        end
        name_ = String(Main.oneAPI.device()[1])
        array_ = Main.oneArray  

    else 
        error("Device: $device not found. Valid device names are ['cpu', 'cuda', 'metal', 'amd', 'oneapi]")
    end

    # test array operations
    try 
        array_([1,2,3]) + array_([1,2,3])
    catch e
        throw(ErrorException("Array-Addition not working. Please ensure your backend is working correctly. $e"))
    end

    # settings
    showinfo = config["settings"]["show_progress_in_console"] 
    float = eval(Symbol(config["settings"]["precision"]))
    N = Int64(config["settings"]["spatial_derivative_order"])
    c_diff = float.(diff_coeff(N))

    # check paths  
    output_path = config["settings"]["output_file"] * ".h5"
    path_dir = dirname(output_path)
    if !isdir(path_dir)
        error("Output directory $path_dir not found. Please create the directory first.")
    end;

    if isfile(output_path)
    @warn """Output $output_path already exists. 
             Once the calculations are complete, this file will be overwritten."""  _module=nothing _file=nothing _line=nothing
    end;

    # types
    T1 = typeof(config)
    T2 = typeof(c_diff)
    T3 = typeof(array_)

    settings = Settings{T1,T2,T3}(config, showinfo, float, output_path, N, c_diff, device, name_, array_)
    
    return settings
end;
