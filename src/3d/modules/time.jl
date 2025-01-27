struct Time{U, V}
    t0::U
    tend::U
    dt::U
    nt::Int64
    t::V
end

function init_time(settings::Settings,
                   elastic::Elastic,
                   domain::Domain,
                   mlog)

    # read 
    t0 = settings.config["time"]["start"]
    tend = settings.config["time"]["end"]
    dt = settings.config["time"]["timestep"]

    t0 = settings.float(t0)
    tend = settings.float(tend)

    # check stability criterion
    courant = 0.65
    dt_stable = courant/(elastic.vpmax * sqrt(1/domain.dx^2 + 1.0/domain.dy^2 + 1.0/domain.dz^2))

    if dt > dt_stable
        add_message!(mlog, "Warning: Given dt does not satisfy the stability criterion. dt changed to $dt_stable")
        dt = dt_stable
    end

    # time vector
    dt = settings.float(dt)
    t = t0:dt:tend
    nt = length(t)

    # types 
    U = typeof(t0)
    V = typeof(t)
    time = Time{U, V}(t0,tend,dt,nt,t) 

    return time
end;
