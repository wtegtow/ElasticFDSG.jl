struct Source{U, V}
    x::U
    y::U
    fdom::U
    xid::Int64
    yid::Int64
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

    fdom = settings.config["source"]["dominant_frequency"] 
    ts = settings.config["source"]["wavelet_center"]   
    f0 = settings.config["source"]["amplitude"] 

    source_x = settings.float(settings.config["source"]["location"]["x"])
    source_y = settings.float(settings.config["source"]["location"]["y"])

    # is source within the domain?
    if source_x >= x0   && 
       source_x <= xend && 
       source_y >= y0   && 
       source_y <= yend
    ;
    else
        throw(DimensionMismatch("Source location is not within the defined domain."))
    end

    # is wavelet truncated?
    if ts < (1 / fdom)
        add_message!(Log, "Warning: Wavelet is not fully included in the source-time function (ts !< 1/fdom).")
    end

    # indices of source
    sx_id = argmin(abs.(domain.xcoords .- source_x))
    sy_id = argmin(abs.(domain.ycoords .- source_y))
    
    # wavelet
    wavelet_type = settings.config["source"]["wavelet_type"]
    if wavelet_type == "ricker"
        STF = ricker(time.t, ts, fdom) .* f0
    elseif wavelet_type == "gauss1d"
        STF = gauss1d(time.t, ts, fdom) .* f0
    end; 

    rhosrc = elastic.rho[sy_id,sx_id]

    # types 
    U = typeof(source_x)
    V = typeof(STF)
    W = typeof(rhosrc)
    source = Source{U, V}(source_x, source_y, fdom, sx_id, sy_id, rhosrc, STF)

    return source 
end;
