mutable struct Fields{T1}
    vx::T1
    vy::T1  
    vz::T1
    sxx::T1
    sxy::T1
    sxz::T1
    syy::T1
    syz::T1
    szz::T1
end

function init_fields(settings::Settings, domain::Domain)
    # allocate primary fields
    vx = zeros(settings.float,domain.dim);
    vy = zeros(settings.float,domain.dim);
    vz = zeros(settings.float,domain.dim);
    sxx = zeros(settings.float,domain.dim);
    sxy = zeros(settings.float,domain.dim);
    sxz = zeros(settings.float,domain.dim);
    syy = zeros(settings.float,domain.dim);
    syz = zeros(settings.float,domain.dim);
    szz = zeros(settings.float,domain.dim);
    # type
    T1 = typeof(vx)
    fields = Fields{T1}(vx, vy, vz, sxx, sxy, sxz, syy, syz, szz)
    return fields 
end