# ============================================================
# Config parser 
# ============================================================
struct Config
    dict::Dict{String, Any}
    dim::Int
end

const _REQUIRED_SETTINGS = [
    ("device",                   String),
    ("precision",                ["Float32", "Float64"]),
    ("spatial_derivative_order", Int),
    ("verbose", Bool),
    ("output_file",              Union{String, Nothing}),
]

const _REQUIRED_TIME = [
    ("start",    Real),
    ("end",      Real),
    ("timestep", Real),
]

const _REQUIRED_SOURCE_2D = [
    ("dominant_frequency", Real),
    ("wavelet_type",       ["ricker", "gauss1d"]),
    ("wavelet_center",     Real),
    ("seismic_moment",     Real),
    ("location", [
        ("x", Real),
        ("z", Real),
    ]),
    ("moment_tensor", [
        ("Mxx",        Real),
        ("Mxz",        Real),
        ("Mzz",        Real),
        ("anisotropic", Bool),
    ]),
]

const _REQUIRED_SOURCE_3D = [
    ("dominant_frequency", Real),
    ("wavelet_type",       ["ricker", "gauss1d"]),
    ("wavelet_center",     Real),
    ("seismic_moment",     Real),
    ("location", [
        ("x", Real),
        ("y", Real),
        ("z", Real),
    ]),
    ("moment_tensor", [
        ("Mxx",        Real),
        ("Mxy",        Real),
        ("Mxz",        Real),
        ("Myy",        Real),
        ("Myz",        Real),
        ("Mzz",        Real),
        ("anisotropic", Bool),
    ]),
]

const _REQUIRED_BOUNDARIES_2D = [
    ("xstart",    String),
    ("xend",      String),
    ("zstart",    String),
    ("zend",      String),
    ("pml_layer", Int),
]

const _REQUIRED_BOUNDARIES_3D = [
    ("xstart",    String),
    ("xend",      String),
    ("ystart",    String),
    ("yend",      String),
    ("zstart",    String),
    ("zend",      String),
    ("pml_layer", Int),
]

const _REQUIRED_RECEIVERS_2D = [
    ("geophones", Any),
    ("das", [
        ("x_aligned", Any),
        ("z_aligned", Any),
    ]),
    ("snapshots", [
        ("times",  Any),
        ("fields", Any),
    ]),
]

const _REQUIRED_RECEIVERS_3D = [
    ("geophones", Any),
    ("das", [
        ("x_aligned", Any),
        ("y_aligned", Any),
        ("z_aligned", Any),
    ]),
    ("snapshots", [
        ("plane_positions", Any),
        ("times",           Any),
        ("fields",          Any),
    ]),
]


function _check_block(block::Dict, schema, context::String="")
    for (key, expected) in schema
        if !haskey(block, key)
            error("Config check failed: missing key \"$key\"$(isempty(context) ? "" : " in \"$context\"")")
        end
        value = block[key]

        if isa(expected, Vector) && all(x -> x isa Tuple, expected)
            # sub-block: recurse
            if !(value isa Dict)
                error("Config check failed: \"$key\" must be a sub-block (Dict), got $(typeof(value))")
            end
            _check_block(value, expected, key)

        elseif isa(expected, Vector)
            # enum check
            if !(value in expected)
                error("Config check failed: value \"$value\" for key \"$key\" not in allowed values $expected")
            end

        elseif expected === Any
            nothing  # no constraint

        else
            # type check
            if !(value isa expected)
                error("Config check failed: \"$key\" must be a $(expected), got $(typeof(value))")
            end
        end
    end
end


function _check_config(dict::Dict, dim::Int)
    for section in ("settings", "time", "source", "boundaries", "receivers")
        if !haskey(dict, section)
            error("Config check failed: missing top-level key \"$section\"")
        end
    end

    _check_block(dict["settings"],   _REQUIRED_SETTINGS,               "settings")
    _check_block(dict["time"],       _REQUIRED_TIME,                   "time")
    _check_block(dict["source"],     dim == 2 ? _REQUIRED_SOURCE_2D     : _REQUIRED_SOURCE_3D,     "source")
    _check_block(dict["boundaries"], dim == 2 ? _REQUIRED_BOUNDARIES_2D : _REQUIRED_BOUNDARIES_3D, "boundaries")
    _check_block(dict["receivers"],  dim == 2 ? _REQUIRED_RECEIVERS_2D  : _REQUIRED_RECEIVERS_3D,  "receivers")
end


function _detect_dim(dict::Dict)::Int
    try
        return haskey(dict["boundaries"], "ystart") ? 3 : 2
    catch
        return 2
    end
end

function _load_dict(input::String)
    @assert ispath(input) "Path to config file does not exist, got: $input"
    supported_formats = [".yaml", ".yml"]

    if endswith(lowercase(input), ".yaml") || endswith(lowercase(input), ".yml")
        return YAML.load_file(input)
    else
        error("Unsupported config file format. Only $(join(supported_formats, ", ")) files are supported.")
    end
end

function parse_config(input::Union{String, Dict})::Config
    dict = if input isa String
        _load_dict(input)
    else
        input
    end
    dim = _detect_dim(dict)
    _check_config(dict, dim)

    # additional checks
    path = get(dict["settings"], "output_file", nothing)

    if !(isnothing(path) || path == "null") 
        output_file = String(path)
        output_dir = dirname(output_file)

        if !isdir(output_dir)
            # we dont allow creating new directories
            error("Output directory does not exist: $output_dir")
        end

        if !endswith(lowercase(output_file), ".h5")
            error("Unsupported output file format: $output_file. Only HDF5 (.h5) is supported.")
        end

        if isfile(output_file)            
            @warn "Output file: $output_file already exists and will be overwritten"    
        end
    end



    return Config(dict, dim)
end


# ============================================================
# Velocity model parser
# ============================================================

const _VELMOD_NFIELDS_2D = 7
const _VELMOD_NFIELDS_3D = 13

struct VelocityModel2D{T<:AbstractArray}
    X::T        # x-coordinate grid  (nx, nz)
    Z::T        # z-coordinate grid  (nx, nz)
    vp::T
    vs::T
    rho::T
    eps::T      # Thomsen epsilon
    del::T      # Thomsen delta
end

struct VelocityModel3D{T<:AbstractArray}
    X::T        # (nx, ny, nz)
    Y::T
    Z::T
    vp::T
    vs::T
    rho::T
    eps1::T
    eps2::T
    gam1::T
    gam2::T
    del1::T
    del2::T
    del3::T
end

struct VelocityModel 
    fields::Union{VelocityModel2D, VelocityModel3D}
    dim::Int
end

function _load_velmod_file(path::String)::AbstractArray
    @assert ispath(path) "Velocity model file not found: $path"

    supported_formats = [".npy", ".npz", ".jld2"]

    if endswith(path, ".npz")
        file = npzread(path)
        @assert haskey(file, "velmod") "NPZ file does not contain key \"velmod\""
        return file["velmod"]

    elseif endswith(path, ".npy")
        return npzread(path)

    elseif endswith(path, ".jld2")
        file = jldopen(path)
        @assert haskey(file, "velmod") "JLD2 file does not contain key \"velmod\""
        arr = file["velmod"]
        close(file)
        return arr

    else
        error("Unsupported velocity model format. Supported: $(join(supported_formats, ", "))")
    end
end

function _check_velmod(arr::AbstractArray)
    nd = ndims(arr)

    if nd == 3
        nf = size(arr, 1)
        if nf != _VELMOD_NFIELDS_2D
            error("Velmod check failed: 3-D array (2D sim) must have $_VELMOD_NFIELDS_2D fields in dim 1, got $nf")
        end

    elseif nd == 4
        nf = size(arr, 1)
        if nf != _VELMOD_NFIELDS_3D
            error("Velmod check failed: 4-D array (3D sim) must have $_VELMOD_NFIELDS_3D fields in dim 1, got $nf")
        end
    else
        error("Velmod check failed: expected a 3-D (2D sim) or 4-D (3D sim) array, got $(nd)-D")
    end

    if !all(isfinite, arr)
        error("Velmod check failed: array contains NaN or Inf values")
    end
end

function _build_velmod(arr::AbstractArray)
    nd = ndims(arr)
    if nd == 3
        return VelocityModel(
            VelocityModel2D(
                arr[1, :, :],   # X
                arr[2, :, :],   # Z
                arr[3, :, :],   # vp
                arr[4, :, :],   # vs
                arr[5, :, :],   # rho
                arr[6, :, :],   # eps
                arr[7, :, :],   # del
            ),
            2
        )
    else
        return VelocityModel(
            VelocityModel3D(
                arr[1,  :, :, :],   # X
                arr[2,  :, :, :],   # Y
                arr[3,  :, :, :],   # Z
                arr[4,  :, :, :],   # vp
                arr[5,  :, :, :],   # vs
                arr[6,  :, :, :],   # rho
                arr[7,  :, :, :],   # eps1
                arr[8,  :, :, :],   # eps2
                arr[9,  :, :, :],   # gam1
                arr[10, :, :, :],   # gam2
                arr[11, :, :, :],   # del1
                arr[12, :, :, :],   # del2
                arr[13, :, :, :],   # del3
            ),
            3
        )
    end
end

function parse_velmod(input::Union{String, AbstractArray})
    arr = if input isa String
        _load_velmod_file(input)
    else
        input
    end
    _check_velmod(arr)
    return _build_velmod(arr)
end


# ============================================================
# Device parser 
# ============================================================

struct Device
    name::String        # human-readable device name
    backend             # KernelAbstractions backend instance
    array               # Array constructor (Array, CuArray, MtlArray, …)
end

function parse_device(config::Config)::Device
    device_str = String(config.dict["settings"]["device"])

    # CPU 
    if device_str in ["cpu", "CPU", "Cpu"]
        return Device("cpu", KernelAbstractions.CPU(), Array)    
        
    # NVIDIA CUDA
    elseif device_str in ["cuda", "gpu-cuda", "gpu_cuda", "CUDA", "Cuda"]
        isdefined(Main, :CUDA) || error("CUDA.jl not loaded in Main. Run `using CUDA` first.")
        devs = Main.CUDA.devices()
        isempty(devs) && error("No CUDA device found.")
        name = String(Main.CUDA.name(Main.CUDA.device!(0)))
        return Device(name, Main.CUDA.CUDABackend(), Main.CuArray)

    # Apple Metal
    elseif device_str in ["metal", "gpu-metal", "gpu_metal", "Metal", "apple", "Apple", "Mac"]
        isdefined(Main, :Metal) || error("Metal.jl not loaded in Main. Run `using Metal` first.")
        devs = Main.Metal.MTL.devices()
        isempty(devs) && error("No Metal device found.")
        name = String(devs[1].name)
        return Device(name, Main.Metal.MetalBackend(), Main.MtlArray)

    # AMD ROCm (not tested, but added for completeness)
    elseif device_str in ["amd", "rocm", "gpu-amd", "gpu_rocm", "AMD"]
        isdefined(Main, :AMDGPU) || error("AMDGPU.jl not loaded in Main. Run `using AMDGPU` first.")
        devs = Main.AMDGPU.devices()
        isempty(devs) && error("No AMD device found.")
        name = String(devs[1])
        return Device(name, Main.AMDGPU.ROCBackend(), Main.ROCArray)

    # Intel oneAPI (not tested, but added for completeness)
    elseif device_str in ["intel", "oneapi", "one_api", "gpu-intel"]
        isdefined(Main, :oneAPI) || error("oneAPI.jl not loaded in Main. Run `using oneAPI` first.")
        dev = Main.oneAPI.device()
        isnothing(dev) && error("No oneAPI device found.")
        name = String(dev)
        return Device(name, Main.oneAPI.oneAPIBackend(), Main.oneArray)

    else
        error("Unknown device: \"$device_str\". Valid options: cpu, cuda, metal, amd, oneapi")
    end
end