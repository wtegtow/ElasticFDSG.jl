function _nearest_id(val, coords)
    argmin(abs.(coords .- val))
end

function _check_in_domain(locs::Matrix, domain::Domain)
    for k in axes(locs, 1)
        c = domain.coordinates[k]
        if any(locs[k,:] .< first(c)) || any(locs[k,:] .> last(c))
            error("Some receivers are outside the domain along axis $k.")
        end
    end
end

function _axis_keys(N)
    N == 2 ? ("x", "z") : ("x", "y", "z")
end

function _field_dict(fields::Fields2D)
    Dict("vx"  => fields.vx,  "vz"  => fields.vz,
         "sxx" => fields.sxx, "sxz" => fields.sxz, "szz" => fields.szz)
end

function _field_dict(fields::Fields3D)
    Dict("vx"  => fields.vx,  "vy"  => fields.vy,  "vz"  => fields.vz,
         "sxx" => fields.sxx, "sxy" => fields.sxy, "sxz" => fields.sxz,
         "syy" => fields.syy, "syz" => fields.syz, "szz" => fields.szz)
end

# ============================================================
# Geophones
# ============================================================

mutable struct Geophones{T<:AbstractFloat}
    n::Int
    ids::Matrix{Int}     # (ndim, ngeo)
    data::Array{T,3}     # (ngeo, ncomp, nt)   ncomp = ndim
    coords::Matrix{T}    # (ndim, ngeo)  — nearest neighbor grid coordinates
end

function _build_locs(gcfg, N, fp)
    keys = _axis_keys(N)
    ngeo = length(gcfg)
    return [fp(gcfg[i][keys[k]]) for k in 1:N, i in 1:ngeo]
end

function init_geophones(config::Config, domain::Domain{N}, time::SimTime) where N
    fp   = eval(Symbol(config.dict["settings"]["precision"]))
    gcfg = something(get(get(config.dict, "receivers", Dict()), "geophones", nothing), [])

    if isempty(gcfg)
        return Geophones{fp}(0, zeros(Int, N, 0), zeros(fp, 0, N, 0), zeros(fp, N, 0))
    end

    ngeo = length(gcfg)
    locs = _build_locs(gcfg, N, fp)
    _check_in_domain(locs, domain)

    ids    = [_nearest_id(locs[k, i], domain.coordinates[k]) for k in 1:N, i in 1:ngeo]
    coords = [domain.coordinates[k][ids[k, i]] for k in 1:N, i in 1:ngeo]
    data   = zeros(fp, ngeo, N, time.nt)
    return Geophones(ngeo, ids, data, coords)
end

function save_geophones!(geo::Geophones, fields::Fields2D, ti)
    geo.n == 0 && return
    for n in 1:geo.n, (c, v) in enumerate((fields.vx, fields.vz))
        geo.data[n, c, ti] = v[geo.ids[1,n], geo.ids[2,n]]
    end
end

function save_geophones!(geo::Geophones, fields::Fields3D, ti)
    geo.n == 0 && return
    for n in 1:geo.n, (c, v) in enumerate((fields.vx, fields.vy, fields.vz))
        geo.data[n, c, ti] = v[geo.ids[1,n], geo.ids[2,n], geo.ids[3,n]]
    end
end


# ============================================================
# DAS
# ============================================================

mutable struct Fiber{T<:AbstractFloat}       
    axis::String        # "axis"_aligned
    coords::Matrix{T}   # (ndim, nch) 
    ids::Matrix{Int}    # (ndim, nch)
    data::Matrix{T}     # (nch, nt)
end

mutable struct DAS
    fibers::Union{Vector{Fiber}, Nothing}
end

function init_das(config::Config, domain::Domain, time::SimTime)

    rcv_cfg = config.dict["receivers"]
    das_cfg = get(rcv_cfg, "das", nothing)
    isnothing(das_cfg) && return DAS(nothing)

    fp = eval(Symbol(config.dict["settings"]["precision"]))
    dim = length(domain.coordinates)

    axis = dim == 2 ? ("x", "z") : ("x", "y", "z")
    axis_keys = dim == 2 ? ("x_aligned", "z_aligned") : 
                           ("x_aligned", "y_aligned", "z_aligned")

    fibers = Fiber[]
    for key in axis_keys
        !haskey(das_cfg, key) && continue 

        for fbr in das_cfg[key]
            
            # this loop only determines nchannel and the aligned axis
            pts = nothing 
            axs = nothing
            for a in axis
                if fbr[a] isa Dict
                    pts = fp.(collect(fbr[a]["start"]:fbr[a]["step"]:fbr[a]["end"]))
                    axs = a 
                    break
                end 
            end

            nchannel = length(pts)
            coords = zeros(fp, dim, nchannel)
            ids    = zeros(Int, dim, nchannel)
            data   = zeros(fp, nchannel, time.nt)

            # locations 
            for (i, a) in enumerate(axis)
                if a == axs 
                    coords[i, :] .= pts
                else
                    coords[i, :] .= fp.(repeat([fbr[a]], nchannel))
                end
            end
            _check_in_domain(coords, domain)

            # location indices
            for i in 1:dim , ch in 1:nchannel
                ids[i, ch] = _nearest_id(coords[i, ch], domain.coordinates[i])
            end

            push!(fibers, Fiber(key, coords, ids, data))

        end
    end
    @logger :debug "$(length(fibers)) Fibers registered"
    das = DAS(fibers)
    return das
end

function save_das!(das::DAS, fields::Fields2D, domain::Domain, N, ti)
    
    # strain rate       =      1/2 * (∇v + ∇vᵀ)
    # axial strain rate = nᵀ · 1/2 * (∇v + ∇vᵀ) · n
    # x-aligned: n=(1,0) -> ∂vx/∂x
    # z-aligned: n=(0,1) -> ∂vz/∂z

    isnothing(das.fibers) && return 

    vx = Array(fields.vx) # copy to cpu to avoid scalar indexing issue
    vz = Array(fields.vz)
    dx = step(domain.coordinates[1])
    dz = step(domain.coordinates[2])
    T = eltype(vx)
    c_fd = diff_coeff(N)
   
    for fiber in das.fibers

        for (ic, (x, z)) in enumerate(eachcol(fiber.ids))

            if fiber.axis == "x_aligned"
                vx_x = zero(T)
                for i in 1:N
                    vx_x += c_fd[i]/dx * (vx[x+i, z] - vx[x-(i-1), z])
                end
                fiber.data[ic, ti] = vx_x

            elseif fiber.axis == "z_aligned"
                vz_z = zero(T)
                for i in 1:N
                    vz_z += c_fd[i]/dz * (vz[x, z+(i-1)] - vz[x, z-i])
                end
                fiber.data[ic, ti] = vz_z 
            end 
        end
    end
end

function save_das!(das::DAS, fields::Fields3D, domain::Domain, N, ti)
    
    # strain rate       =      1/2 * (∇v + ∇vᵀ)
    # axial strain rate = nᵀ · 1/2 * (∇v + ∇vᵀ) · n
    # x-aligned: n=(1,0,0) -> ∂vx/∂x
    # y-aligned: n=(0,1,0) -> ∂vy/∂y
    # z-aligned: n=(0,0,1) -> ∂vz/∂z

    isnothing(das.fibers) && return 

    vx = Array(fields.vx) # copy to cpu to avoid scalar indexing issue
    vy = Array(fields.vy)
    vz = Array(fields.vz)
    dx = step(domain.coordinates[1])
    dy = step(domain.coordinates[2])
    dz = step(domain.coordinates[3])
    T = eltype(vx)
    c_fd = diff_coeff(N)
   
    for fiber in das.fibers

        for (ic, (x, y, z)) in enumerate(eachcol(fiber.ids))

            if fiber.axis == "x_aligned"
                vx_x = zero(T)
                for i in 1:N
                    vx_x += c_fd[i]/dx * (vx[x+i, y, z] - vx[x-(i-1), y, z])
                end
                fiber.data[ic, ti] = vx_x

            elseif fiber.axis == "y_aligned"
                vy_y = zero(T)
                for i in 1:N
                    vy_y += c_fd[i]/dy * (vy[x, y+(i-1), z] - vy[x, y-i, z])
                end
                fiber.data[ic, ti] = vy_y

            elseif fiber.axis == "z_aligned"
                vz_z = zero(T)
                for i in 1:N
                    vz_z += c_fd[i]/dz * (vz[x, y, z+(i-1)] - vz[x, y, z-i])
                end
                fiber.data[ic, ti] = vz_z 
            end 
        end
    end
end

# ============================================================
# Snapshots
# ============================================================

mutable struct Snapshots2D{T<:AbstractFloat}
    n::Int
    data::Array{T,4}           # (ntime, nfields, nx, nz)
    fieldnames::Vector{String}
    tid_map::Dict{Int,Int}
end

function init_snapshots(config::Config, domain::Domain{2}, time::SimTime)
    fp    = eval(Symbol(config.dict["settings"]["precision"]))

    scfg = get(config.dict["receivers"], "snapshots", nothing)
    if isnothing(scfg) # default case if empty
        return Snapshots2D{fp}(0, zeros(fp, 0, 0, 0, 0), String[], Dict{Int,Int}())
    end
    tsnap = scfg["times"]
    flds  = scfg["fields"]

    nx, nz  = domain.shape
    ntime   = length(tsnap)
    nfields = length(flds)
    tids    = [_nearest_id(t, time.t) for t in tsnap]
    tid_map = Dict(tids[i] => i for i in eachindex(tids))
    data    = zeros(fp, ntime, nfields, nx, nz)

    return Snapshots2D(ntime * nfields, data, String.(flds), tid_map)
end

function save_snapshots!(snap::Snapshots2D, fields::Fields2D, ti)
    
    # new guard
    snap.n == 0 && return
    haskey(snap.tid_map, ti) || return

    tid = snap.tid_map[ti]
    fd  = _field_dict(fields)
    for (fi, name) in enumerate(snap.fieldnames)
        snap.data[tid, fi, :, :] .= Array(fd[name]) # Array() to avoid GPU-host issues if fields on device
    end
end

mutable struct Snapshots3D{T<:AbstractFloat}
    n::Int
    XY::Array{T,5}             # (nsnaps, ntime, nfields, nx, ny)
    XZ::Array{T,5}             # (nsnaps, ntime, nfields, nx, nz)
    YZ::Array{T,5}             # (nsnaps, ntime, nfields, ny, nz)
    fieldnames::Vector{String}
    tid_map::Dict{Int,Int}
    grid_ids::Vector{Vector{Int}}  # per-snap: [ix, iy, iz]
end

function init_snapshots(config::Config, domain::Domain{3}, time::SimTime)
    fp     = eval(Symbol(config.dict["settings"]["precision"]))

    scfg = get(config.dict["receivers"], "snapshots", nothing)
    if isnothing(scfg) # default case if empty
        empty5 = zeros(fp, 0, 0, 0, 0, 0)
        return Snapshots3D{fp}(0, empty5, empty5, empty5, String[], Dict{Int,Int}(), Vector{Int}[])
    end
    tsnap   = scfg["times"]
    flds    = scfg["fields"]
    planes  = scfg["plane_positions"]

    nx, ny, nz = domain.shape
    nsnaps  = length(planes)
    ntime   = length(tsnap)
    nfields = length(flds)
    tids    = [_nearest_id(t, time.t) for t in tsnap]
    tid_map = Dict(tids[i] => i for i in eachindex(tids))
    grid_ids = [[_nearest_id(fp(p["x"]), domain.coordinates[1]),
                 _nearest_id(fp(p["y"]), domain.coordinates[2]),
                 _nearest_id(fp(p["z"]), domain.coordinates[3])] for p in planes]

    XY = zeros(fp, nsnaps, ntime, nfields, nx, ny)
    XZ = zeros(fp, nsnaps, ntime, nfields, nx, nz)
    YZ = zeros(fp, nsnaps, ntime, nfields, ny, nz)

    return Snapshots3D(nsnaps * ntime * nfields, XY, XZ, YZ, String.(flds), tid_map, grid_ids)
end

function save_snapshots!(snap::Snapshots3D, fields::Fields3D, ti)
    
    # new guard
    snap.n == 0 && return
    haskey(snap.tid_map, ti) || return

    tid = snap.tid_map[ti]
    fd  = _field_dict(fields)
    for n in 1:length(snap.grid_ids)
        ix, iy, iz = snap.grid_ids[n]
        for (fi, name) in enumerate(snap.fieldnames)
            f = fd[name]
            snap.XY[n, tid, fi, :, :] .= Array(f[:, :, iz]) # Array() to ensure GPU compability
            snap.XZ[n, tid, fi, :, :] .= Array(f[:, iy, :])
            snap.YZ[n, tid, fi, :, :] .= Array(f[ix, :, :])
        end
    end
end

function init_receiver(config::Config, domain::Domain, elastic::Elastic, time::SimTime)

    geophones = init_geophones(config, domain, time)
    das       = init_das(config, domain, time)
    dasn = isnothing(das.fibers) ? 0 : length(das.fibers)
    snapshots = init_snapshots(config, domain, time)

    nrec = geophones.n + snapshots.n + dasn
    if nrec == 0
        @logger :warn "Receiver list is empty. No data will be saved."
    end
    @logger :debug "All receivers registered"
    return geophones, das, snapshots
end
