mutable struct Fields2D{T<:AbstractArray}
    vx::T;  vz::T
    sxx::T; sxz::T; szz::T
end

mutable struct Fields3D{T<:AbstractArray}
    vx::T;  vy::T;  vz::T
    sxx::T; sxy::T; sxz::T
    syy::T; syz::T; szz::T
end

function init_fields(config::Config, domain::Domain{2})
    fp = eval(Symbol(config.dict["settings"]["precision"]))
    z  = () -> zeros(fp, domain.shape)
    return Fields2D(z(), z(), z(), z(), z())
end

function init_fields(config::Config, domain::Domain{3})
    fp = eval(Symbol(config.dict["settings"]["precision"]))
    z  = () -> zeros(fp, domain.shape)
    return Fields3D(z(), z(), z(), z(), z(), z(), z(), z(), z())
end