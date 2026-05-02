# ============================================================
# Shared helpers
# ============================================================

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

function _axis_keys(N)
    N == 2 ? ("x", "z") : ("x", "y", "z")
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
    GPUArrays.@allowscalar for n in 1:geo.n, (c, v) in enumerate((fields.vx, fields.vz))
        geo.data[n, c, ti] = v[geo.ids[1,n], geo.ids[2,n]]
    end
end

function save_geophones!(geo::Geophones, fields::Fields3D, ti)
    geo.n == 0 && return
    GPUArrays.@allowscalar for n in 1:geo.n, (c, v) in enumerate((fields.vx, fields.vy, fields.vz))
        geo.data[n, c, ti] = v[geo.ids[1,n], geo.ids[2,n], geo.ids[3,n]]
    end
end


# ============================================================
# DAS
# ============================================================

mutable struct Fibers{T<:AbstractFloat}
    n::Int
    axis::Int                    # varying spatial axis (1=x, 2=y or z in 2D, 3=z in 3D)
    ids::Vector{Any}             # per-fiber: Vector{Any}[ndim], axis k: Vector{Int} or Int
    cinv::Vector{Array{T,3}}     # per-fiber: (nch, ndim, ndim) — upper-left normal-stress compliance
    data::Vector{Matrix{T}}      # per-fiber: (nch, nt)
    coords::Vector{Matrix{T}}    # per-fiber: (ndim, nch) — nearest neighbor grid coordinates
end

mutable struct DAS
    fibers::Vector{Fibers}       # one per spatial axis
end

function _cinv_from_stiffness(s::Stiffness2D, fp)
    s.c44 == 0 && return fp.(diagm([1/s.c11, 1/s.c33]))
    return fp.(inv([s.c11 s.c13; s.c13 s.c33]))
end

function _cinv_from_stiffness(s::Stiffness3D, fp)
    (s.c55 == 0 || s.c44 == 0) && return fp.(diagm([1/s.c11, 1/s.c22, 1/s.c33]))
    return fp.(inv([s.c11 s.c12 s.c13;
                    s.c12 s.c22 s.c23;
                    s.c13 s.c23 s.c33]))
end

function _register_fibers(fp, axis, fibers_cfg, axis_keys, domain, elastic, nt)
    ndim = length(axis_keys)
    if isnothing(fibers_cfg) || isempty(fibers_cfg)
        return Fibers{fp}(0, axis, Any[], Array{fp,3}[], Matrix{fp}[], Matrix{fp}[])
    end

    nfibers    = length(fibers_cfg)
    all_ids    = Vector{Any}(undef, nfibers)
    all_cinv   = Vector{Array{fp,3}}(undef, nfibers)
    all_data   = Vector{Matrix{fp}}(undef, nfibers)
    all_coords = Vector{Matrix{fp}}(undef, nfibers)

    for i in 1:nfibers
        fiber = fibers_cfg[i]

        # per-axis ids: Vector{Int} for varying axis, Int for fixed axes
        ids_i = Vector{Any}(undef, ndim)
        for (k, key) in enumerate(axis_keys)
            raw = fiber[key]
            if k == axis
                pts      = fp.(collect(raw["start"]:raw["step"]:raw["end"]))
                ids_i[k] = [_nearest_id(v, domain.coordinates[k]) for v in pts]
            else
                ids_i[k] = _nearest_id(fp(raw), domain.coordinates[k])
            end
        end
        all_ids[i] = ids_i

        # build Cinv for each channel
        # iterate axes in (x,y,z) order; fixed axes wrapped in 1-element vector
        nch    = length(ids_i[axis])
        cinv_i = zeros(fp, nch, ndim, ndim)
        iter   = [k == axis ? ids_i[k] : [ids_i[k]] for k in 1:ndim]
        for (ch, idx) in enumerate(Iterators.product(iter...))
            s = elastic.c_tensors[elastic.c_lookup[idx...]].fields
            cinv_i[ch, :, :] .= _cinv_from_stiffness(s, fp)
        end
        all_cinv[i] = cinv_i
        all_data[i] = zeros(fp, nch, nt)

        coords_i = zeros(fp, ndim, nch)
        for ch in 1:nch
            for k in 1:ndim
                if k == axis
                    coords_i[k, ch] = domain.coordinates[k][ids_i[k][ch]]
                else
                    coords_i[k, ch] = domain.coordinates[k][ids_i[k]]
                end
            end
        end
        all_coords[i] = coords_i
    end

    return Fibers(nfibers, axis, all_ids, all_cinv, all_data, all_coords)
end

function init_das(config::Config, domain::Domain{2}, elastic::Elastic, time::SimTime)
    fp           = eval(Symbol(config.dict["settings"]["precision"]))
    das_raw      = get(get(config.dict, "receivers", Dict()), "das", nothing)
    das_cfg      = something(das_raw, Dict())
    axis_keys    = _axis_keys(2)
    orientations = ("x_aligned", "z_aligned")
    fibers = [_register_fibers(fp, axis, get(das_cfg, orient, nothing),
                               axis_keys, domain, elastic, time.nt)
              for (axis, orient) in enumerate(orientations)]
    return DAS(fibers)
end

function init_das(config::Config, domain::Domain{3}, elastic::Elastic, time::SimTime)
    fp           = eval(Symbol(config.dict["settings"]["precision"]))
    das_raw      = get(get(config.dict, "receivers", Dict()), "das", nothing)
    das_cfg      = something(das_raw, Dict())
    axis_keys    = _axis_keys(3)
    orientations = ("x_aligned", "y_aligned", "z_aligned")
    fibers = [_register_fibers(fp, axis, get(das_cfg, orient, nothing),
                               axis_keys, domain, elastic, time.nt)
              for (axis, orient) in enumerate(orientations)]
    return DAS(fibers)
end

function save_das!(das::DAS, fields::Fields2D, ti)
    for fg in das.fibers
        fg.n == 0 && continue
        for i in 1:fg.n
            ids   = fg.ids[i]
            sxx_v = vec(Array(fields.sxx[ids[1], ids[2]]))
            szz_v = vec(Array(fields.szz[ids[1], ids[2]]))
            for ch in 1:size(fg.data[i], 1)
                fg.data[i][ch, ti] = dot(fg.cinv[i][ch, fg.axis, :], [sxx_v[ch], szz_v[ch]])
            end
        end
    end
end

function save_das!(das::DAS, fields::Fields3D, ti)
    for fg in das.fibers
        fg.n == 0 && continue
        for i in 1:fg.n
            ids   = fg.ids[i]
            sxx_v = vec(Array(fields.sxx[ids[1], ids[2], ids[3]]))
            syy_v = vec(Array(fields.syy[ids[1], ids[2], ids[3]]))
            szz_v = vec(Array(fields.szz[ids[1], ids[2], ids[3]]))
            for ch in 1:size(fg.data[i], 1)
                fg.data[i][ch, ti] = dot(fg.cinv[i][ch, fg.axis, :], [sxx_v[ch], syy_v[ch], szz_v[ch]])
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

mutable struct Snapshots3D{T<:AbstractFloat}
    n::Int
    XY::Array{T,5}             # (nsnaps, ntime, nfields, nx, ny)
    XZ::Array{T,5}             # (nsnaps, ntime, nfields, nx, nz)
    YZ::Array{T,5}             # (nsnaps, ntime, nfields, ny, nz)
    fieldnames::Vector{String}
    tid_map::Dict{Int,Int}
    grid_ids::Vector{Vector{Int}}  # per-snap: [ix, iy, iz]
end

function init_snapshots(config::Config, domain::Domain{2}, time::SimTime)
    fp    = eval(Symbol(config.dict["settings"]["precision"]))
    scfg  = something(get(get(config.dict, "receivers", Dict()), "snapshots", nothing), Dict())
    tsnap = get(scfg, "times",  nothing)
    flds  = get(scfg, "fields", nothing)

    if isnothing(tsnap) || isnothing(flds)
        return Snapshots2D{fp}(0, zeros(fp, 0, 0, 0, 0), String[], Dict{Int,Int}())
    end

    nx, nz  = domain.shape
    ntime   = length(tsnap)
    nfields = length(flds)
    tids    = [_nearest_id(t, time.t) for t in tsnap]
    tid_map = Dict(tids[i] => i for i in eachindex(tids))
    data    = zeros(fp, ntime, nfields, nx, nz)

    return Snapshots2D(ntime * nfields, data, String.(flds), tid_map)
end

function init_snapshots(config::Config, domain::Domain{3}, time::SimTime)
    fp     = eval(Symbol(config.dict["settings"]["precision"]))
    scfg   = something(get(get(config.dict, "receivers", Dict()), "snapshots", nothing), Dict())
    tsnap  = get(scfg, "times",          nothing)
    flds   = get(scfg, "fields",         nothing)
    planes = get(scfg, "plane_positions", nothing)

    if isnothing(tsnap) || isnothing(flds) || isnothing(planes)
        empty5 = zeros(fp, 0, 0, 0, 0, 0)
        return Snapshots3D{fp}(0, empty5, empty5, empty5, String[], Dict{Int,Int}(), Vector{Int}[])
    end

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

function save_snapshots!(snap::Snapshots2D, fields::Fields2D, ti)
    snap.n == 0 || !haskey(snap.tid_map, ti) && return
    tid = snap.tid_map[ti]
    fd  = _field_dict(fields)
    for (fi, name) in enumerate(snap.fieldnames)
        snap.data[tid, fi, :, :] .= Array(fd[name]) # Array() to avoid GPU-host issues if fields on device
    end
end

function save_snapshots!(snap::Snapshots3D, fields::Fields3D, ti)
    snap.n == 0 || !haskey(snap.tid_map, ti) && return
    tid = snap.tid_map[ti]
    fd  = _field_dict(fields)
    for n in 1:length(snap.grid_ids)
        ix, iy, iz = snap.grid_ids[n]
        for (fi, name) in enumerate(snap.fieldnames)
            f = fd[name]
            snap.XY[n, tid, fi, :, :] .= Array(f[:, :, iz]) # Array() to avoid GPU-host issues if fields on device
            snap.XZ[n, tid, fi, :, :] .= Array(f[:, iy, :])
            snap.YZ[n, tid, fi, :, :] .= Array(f[ix, :, :])
        end
    end
end


# ============================================================
# Top-level init
# ============================================================

function init_receiver(config::Config, domain::Domain, elastic::Elastic, time::SimTime)
    geophones = init_geophones(config, domain, time)
    das       = init_das(config, domain, elastic, time)
    snapshots = init_snapshots(config, domain, time)

    nrec = geophones.n + sum(f.n for f in das.fibers) + snapshots.n
    if nrec == 0
        @warn "Receiver list is empty. No data will be saved." _module=nothing _file=nothing _line=nothing
    end

    return geophones, das, snapshots
end
