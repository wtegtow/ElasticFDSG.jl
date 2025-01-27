struct Source{U, V}
    x::U
    y::U
    z::U
    xid::Int64
    yid::Int64
    zid::Int64
    fdom::U
    rhosrc::U
    STF::V
end

function init_source(settings::Settings,
                    domain::Domain,
                    elastic::Elastic,
                    time::Time,
                    Log)

    x0 = domain.xcoords[begin]
    xend = domain.xcoords[end]
    y0 = domain.ycoords[begin]
    yend = domain.ycoords[end]
    z0 = domain.zcoords[begin]
    zend = domain.zcoords[end]

    fdom = settings.config["source"]["dominant_frequency"] 
    ts = settings.config["source"]["wavelet_center"]   
    f0 = settings.config["source"]["amplitude"] 

    source_x = settings.float(settings.config["source"]["location"]["x"])
    source_y = settings.float(settings.config["source"]["location"]["y"])
    source_z = settings.float(settings.config["source"]["location"]["z"])

    # is source within the domain?
    if source_x >= x0   && 
       source_x <= xend && 
       source_y >= y0   && 
       source_y <= yend &&
       source_z >= z0   && 
       source_z <= zend 
    ;
    else
        throw(DimensionMismatch("Source location is not within the defined domain."))
    end

    # is the wavelet truncated?
    if ts < (1 / fdom)
        add_message!(Log, "Warning: Wavelet is not fully included in the source-time function (ts !< 1/fdom).")
    end

    # indices of source
    sx_id = argmin(abs.(domain.xcoords .- source_x))
    sy_id = argmin(abs.(domain.ycoords .- source_y))
    sz_id = argmin(abs.(domain.zcoords .- source_z))

    # density at source location
    rhosrc = elastic.rho[sz_id, sy_id, sx_id]
    rhosrc = settings.float(rhosrc)
    
    # wavelet
    wavelet_type = settings.config["source"]["wavelet_type"]
    # source time function 
    if wavelet_type == "ricker"
        STF = ricker(time.t, ts, fdom) .* f0
    elseif wavelet_type == "gauss1d"
        STF = gauss1d(time.t, ts, fdom) .* f0
    end; 
    STF = settings.float.(STF)

    # types 
    U = typeof(source_x)
    V = typeof(STF)
    W = typeof(rhosrc)
    source = Source{U, V}(source_x, source_y, source_z,
                          sx_id, sy_id, sz_id,
                          fdom, rhosrc, STF)

    return source 
    
end;