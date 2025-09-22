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
    zeros::Function
    ones::Function
end


function init_settings(CONFIGPATH)

    # read and check config file
    config = YAML.load_file(CONFIGPATH)
    check_configfile2d(config) 

    # set device with corresponding array types
    device = String(config["settings"]["device"]) 
    name_  = nothing
    array_ = nothing 
    zeros_ = nothing 
    ones_  = nothing

    # cpu
    if device in ["cpu", "CPU","Cpu"]
        device = "cpu"
        name_  = "cpu"
        array_ = Base.Array
        zeros_ = (dims...) -> Base.zeros(float, dims...) 
        ones_  = (dims...) -> Base.ones(float, dims...)
    
    # nvidia gpu
    elseif device in ["gpu-cuda", "gpu_cuda", "cuda", "Cuda"]
        device = "gpu-cuda"
        device_ids = CUDA.devices()
        if isempty(device_ids)
            error("Cuda device not found.")
        else
            name_  = String(CUDA.name(CUDA.device!(0)))
            array_ = CuArray
            zeros_ = (dims...) -> CUDA.zeros(float, dims...) 
            ones_  = (dims...) -> CUDA.ones(float, dims...) 
        end

    # apple gpu
    elseif device in ["gpu-metal", "gpu_metal", "metal", "apple", "Apple", "MacBook", "Mac"]
        device = "gpu-metal"
        device_list = Metal.MTL.devices()

        if isempty(device_list)
            error("Metal device not found.")
        else 
            name_  = String(device_list[1].name)
            array_ = MtlArray
            zeros_ = (dims...) -> Metal.zeros(float, dims...) 
            ones_  = (dims...) -> Metal.ones(float, dims...) 
        end

    # AMD not testet 
    elseif device in ["gpu-amd", "gpu_rocm", "amd", "AMD"]
        device = "gpu-amd"
        try
            AMDGPU.devices()
        catch err
            error("AMDGPU device not found: $err")
        end
        name_ = AMDGPU.devices()[1]
        array_ = ROCArray
        zeros_ = (dims...) -> AMDGPU.zeros(float, dims...)
        ones_  = (dims...) -> AMDGPU.ones(float, dims...)

    # oneAPI not testet yet  
    elseif device in ["gpu-intel", "intel", "oneapi", "one_api"]
        device = "gpu-intel"
        try
            oneAPI.device()
        catch err
            error("oneAPI device not found: $err")
        end
        name_ = oneAPI.device()[1]  
        array_ = oneArray  
        zeros_ = (dims...) -> oneAPI.zeros(float, dims...)
        ones_  = (dims...) -> oneAPI.ones(float, dims...)

    else 
        error("Device: $device not found. Valid device names are ['cpu', 'gpu_metal', 'gpu_cuda']")
    end

    # settings
    showinfo = config["settings"]["show_progress_in_console"] 
    float = eval(Symbol(config["settings"]["precision"]))
    N = Int64(config["settings"]["spatial_derivative_order"])
    c = float.(diff_coeff(N))

    # paths  
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
    T2 = typeof(c)
    T3 = typeof(array_)

    settings = Settings{T1,T2,T3}(config, showinfo, float, output_path, N, c, device, name_, array_, zeros_, ones_)
    
    return settings
end;
