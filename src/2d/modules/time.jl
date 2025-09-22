struct Time{T1,T2}
    t0::T1
    tend::T1
    dt::T1
    nt::Int64
    t::T2
end

function init_time(settings::Settings,
                   elastic::Elastic,
                   domain::Domain)

    # read 
    t0 = settings.config["time"]["start"]
    tend = settings.config["time"]["end"]
    dt = settings.config["time"]["timestep"]

    t0 = settings.float(t0)
    tend = settings.float(tend)

    # check CFL stability criterion
    courant = 0.55
    dt_stable = courant/(elastic.vmax * sqrt(1/domain.dx^2 + 1.0/domain.dz^2))

    if dt > dt_stable
        @warn "Given Δt does not satisfy the CFL-stability criterion. Δt changed to $dt_stable."  _module=nothing _file=nothing _line=nothing
        dt = dt_stable
    end

    # time vector
    dt = settings.float(dt)
    t = t0:dt:tend
    nt = length(t)
    tend = t[end] 

    # types 
    T1 = typeof(t0)
    T2 = typeof(t)
    time = Time{T1,T2}(t0, tend, dt, nt, t) 

    return time
end;
