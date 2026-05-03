struct SimTime{T}
    t0::T
    tend::T
    dt::T
    nt::Int
    t::StepRangeLen{T}
end

function init_time(config::Config, domain::Domain, elastic::Elastic)
    fp   = eval(Symbol(config.dict["settings"]["precision"]))
    tcfg = config.dict["time"]

    t0   = fp(tcfg["start"])
    tend = fp(tcfg["end"])
    dt   = fp(tcfg["timestep"])

    courant    = fp(0.6)
    spacings   = map(c -> fp(step(c)), domain.coordinates)
    dt_stable  = courant / (elastic.vmax * sqrt(sum(1/d^2 for d in spacings)))

    if dt > dt_stable
        @warn "Given Δt=$dt does not satisfy the CFL criterion. Δt changed to $dt_stable." _module=nothing _file=nothing _line=nothing
        dt = dt_stable
    end

    t    = t0:dt:tend+(fp(0.51)*dt)
    tend = last(t)
    nt   = length(t)

    return SimTime(t0, tend, dt, nt, t)
end
