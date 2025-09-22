# =======================================
# geophones
# =======================================
mutable struct Geophones{T1,T2,T3}
    data::T1
    n::Int64
    coords::T2
    xids::T3
    zids::T3
end

function init_geophones(settings::Settings, 
                        domain::Domain,
                        time::Time)

    geophones = nothing

    # if empty
    if isnothing(settings.config["receivers"]["geophones"])
        T = typeof(nothing)
        geophones = Geophones{T,T,T}(nothing, 0, nothing, nothing, nothing)
    else
        ngeo = length(settings.config["receivers"]["geophones"])
        geo_coords = zeros(settings.float, 2, ngeo)
    
        for i in 1:ngeo
            geo_coords[1, i] = settings.config["receivers"]["geophones"][i]["x"]
            geo_coords[2, i] = settings.config["receivers"]["geophones"][i]["z"]
        end 

        geo_x = geo_coords[1,:]
        geo_z = geo_coords[2,:]

        x0 = domain.xcoords[begin]
        xend = domain.xcoords[end]
        z0 = domain.zcoords[begin]
        zend = domain.zcoords[end]

        # are receivers within the domain?
        if all(geo_x .>= x0)   && 
        all(geo_x .<= xend) && 
        all(geo_z .>= z0)   &&
        all(geo_z .<= zend)
            ;
        else
            error("Some geophone receivers are defined outsite the domain.")
        end

        geo_x_ids = [argmin(abs.(domain.xcoords .- geo_xi))  for geo_xi in geo_x]
        geo_z_ids = [argmin(abs.(domain.zcoords .- geo_zi))  for geo_zi in geo_z]

        ngeo = length(geo_x_ids)
        data = zeros(settings.float, ngeo, 2, time.nt)

        T1 = typeof(data)
        T2 = typeof(geo_coords)
        T3 = typeof(geo_x_ids)
        geophones = Geophones{T1, T2, T3}(data, ngeo, geo_coords, geo_x_ids, geo_z_ids)
    end;
    return geophones

end;

function save_geophones!(geophones::Geophones, vx, vz, ti)
    for n in 1:geophones.n  
        geophones.data[n,1,ti] = vx[geophones.zids[n],geophones.xids[n]]
        geophones.data[n,2,ti] = vz[geophones.zids[n],geophones.xids[n]]
    end
end;

# =======================================
# das
# =======================================

mutable struct Fibers{T1,T2,T3,T4}
    n::Int64
    orientation::String
    data::T1
    coords::T2
    ids::T3
    cinv::T4
end

mutable struct DAS
    x_aligned::Fibers 
    z_aligned::Fibers
end

function register_fibers(settings::Settings, domain::Domain, elastic::Elastic, time::Time, fiber_orientation)
 
    fiber_coords = []
    fiber_ids = []
    fiber_cinv = []
    fiber_data = []
    nfibers = 0

    fibers = settings.config["receivers"]["das"][fiber_orientation]
    if isnothing(fibers)
        return 0,fiber_orientation,nothing,nothing,nothing,nothing
    else
        nfibers = length(fibers)

        # coordinates of fibers
        for i in 1:nfibers

            fiber = fibers[i]
            if fiber_orientation == "x_aligned"
                x = fiber["x"] 
                fiber_x = settings.float.(collect(x["start"]:x["step"]:x["end"]))
                fiber_z = fiber["z"]
                push!(fiber_coords, [fiber_x, fiber_z])
            
            elseif fiber_orientation == "z_aligned"
                fiber_x = fiber["x"]
                z = fiber["z"]   
                fiber_z = settings.float.(collect(z["start"]:z["step"]:z["end"]))
                push!(fiber_coords, [fiber_x, fiber_z])
            end
        end

        # indices of fibers
        for i in 1:nfibers
            fiber_x, fiber_z = fiber_coords[i]
            if fiber_orientation == "x_aligned"
                fiber_x_ids = [argmin(abs.(domain.xcoords .- v)) for v in fiber_x]
                fiber_z_ids =  argmin(abs.(domain.zcoords .- fiber_z))
                push!(fiber_ids, [fiber_x_ids, fiber_z_ids])
            
            elseif fiber_orientation == "z_aligned"
                fiber_x_ids =  argmin(abs.(domain.xcoords .- fiber_x))
                fiber_z_ids = [argmin(abs.(domain.zcoords .- v)) for v in fiber_z]
                push!(fiber_ids, [fiber_x_ids, fiber_z_ids])
            end
        end

        # inverted C-Tensors along fiber coordinates
        for k in 1:nfibers
            fiber_x_ids, fiber_z_ids = fiber_ids[k]
            fiber_dim = length(fiber_x_ids) * length(fiber_z_ids)
            Cinv = zeros(settings.float, fiber_dim , 2, 2)
            channel_idx = 1
            for zid in fiber_z_ids, xid in fiber_x_ids 

                c11, c13, c33, c55, rho = elastic.c_tensors[elastic.c_lookup[zid, xid], :]
                is_liquid = c55 == 0 
                if is_liquid
                    Cinv[channel_idx, :, :] .= 0.0   
                    Cinv[channel_idx,1,1] = 1/c11       # exx = sxx / c11
                    Cinv[channel_idx,2,2] = 1/c33       # ezz = szz / c33
                else
                C = [c11 c13;
                     c13 c33]
                Cinv[channel_idx,:,:] .= inv(C)
                end
                channel_idx += 1
            end
            push!(fiber_cinv, Cinv)
        end

        # data arrays 
        for i in 1:nfibers
            fiber_i_dim = length(fiber_ids[i][1]) * length(fiber_ids[i][2])
            data_i = zeros(settings.float, fiber_i_dim, time.nt)
            push!(fiber_data, data_i)
        end

    end
    
    return nfibers, fiber_coords, fiber_ids, fiber_cinv, fiber_data
end;

function init_das(settings::Settings, domain::Domain, elastic::Elastic, time::Time)
    fibers = []
    for orientation in ["x_aligned", "z_aligned"]
        nfibers, fiber_coords, fiber_ids, fiber_cinv, fiber_data = register_fibers(settings, domain, elastic, time, orientation)
        # types 
        T1 = typeof(fiber_data)
        T2 = typeof(fiber_coords)
        T3 = typeof(fiber_ids)
        T4 = typeof(fiber_cinv)
        push!(fibers, Fibers{T1,T2,T3,T4}(nfibers, orientation, fiber_data, fiber_coords, fiber_ids, fiber_cinv))
    end
    # gather all fibers in DAS struct
    das = DAS(fibers...)
    return das 
end;

function save_das!(das::DAS, sxx, szz, ti)
    fibers = [das.x_aligned, das.z_aligned] 
    for (component, fiber) in enumerate(fibers)
        for i in 1:fiber.n 
            coord = (fiber.ids[i][2], fiber.ids[i][1])
            sxx_along_fiber = Array{eltype(fiber.data[i])}(sxx[coord...])
            szz_along_fiber = Array{eltype(fiber.data[i])}(szz[coord...])
            for channel in eachindex(fiber.ids[i][component])
                sigma = [sxx_along_fiber[channel] szz_along_fiber[channel]]
                fiber.data[i][channel, ti] = (sigma * fiber.cinv[i][channel,component,:])[1]
            end
        end 
    end
end

# =======================================
# snapshots
# =======================================

mutable struct Snapshots{T1,T2,T3,T4,T5}
    n::Int64
    XZ::T1
    fieldnames::T2
    t::T3
    tid::T4
    tid_map::T5
end 

function init_snapshots(settings::Settings, domain::Domain, time::Time)

    tsnap = settings.config["receivers"]["snapshots"]["times"]
    fields = settings.config["receivers"]["snapshots"]["fields"]

    # if empty
    if isnothing(tsnap) || isnothing(fields) 
        T = typeof(nothing)
        snapshots = Snapshots{T,T,T,T,T}(0, repeat([nothing],6)...)
    else 

    nfields = length(fields)
    ntime = length(tsnap)
    
    # data arrays for snapshots
    XZ = zeros(settings.float, ntime, nfields, domain.nz, domain.nx);

    # time indices of snapshots
    tsnap_ids = [argmin(abs.(time.t .- tsnap_i)) for tsnap_i in tsnap]
    tid_map = Dict(tsnap_ids[i] => i for i in eachindex(tsnap_ids))

    # types
    T1 = typeof(XZ)
    T2 = typeof(fields)
    T3 = typeof(tsnap)
    T4 = typeof(tsnap_ids)
    T5 = typeof(tid_map)
    snapshots = Snapshots{T1,T2,T3,T4,T5}(nfields*ntime, XZ, fields, tsnap, tsnap_ids, tid_map);
    end 

    return snapshots
end;

function save_snapshots!(snapshots::Snapshots, vx, vz, sxx, sxz, szz, ti)
    
    if haskey(snapshots.tid_map, ti)
        tid = snapshots.tid_map[ti]
        fieldnum = 1
        for field in snapshots.fieldnames
            if field == "vx"
                field_ = vx
            elseif field == "vz"
                field_ = vz 
            elseif field == "sxx"
                field_ = sxx
            elseif field == "sxz"
                field_ = sxz
            elseif field == "szz"
                field_ = szz
            end
            snapshots.XZ[tid, fieldnum, :, :] .= Array{eltype(snapshots.XZ)}(field_)
            fieldnum += 1
        end
    end 
end;

# =======================================
# receiver
# =======================================

function init_receiver(settings::Settings, domain::Domain, elastic::Elastic, time::Time)

    geophones = init_geophones(settings, domain, time); 
    das = init_das(settings, domain, elastic, time);
    snapshots = init_snapshots(settings, domain, time);

    nrec = snapshots.n  +
        geophones.n     +
        das.x_aligned.n +
        das.z_aligned.n

    if nrec == 0
        @warn "Receiver list is empty. No data will be saved." _module=nothing _file=nothing _line=nothing
    end;
    return geophones, das, snapshots
end;