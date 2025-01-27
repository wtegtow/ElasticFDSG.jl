mutable struct Fields{T}
    vx::T 
    vy::T  
    vz::T
    sxx::T 
    sxy::T
    sxz::T 
    syy::T
    syz::T 
    szz::T
end

function init_fields(settings::Settings, domain::Domain)
    
    # allocate field quantities
    vx = zeros(settings.float,domain.dim);
    vy = zeros(settings.float,domain.dim);
    vz = zeros(settings.float,domain.dim);
    sxx = zeros(settings.float,domain.dim);
    sxy = zeros(settings.float,domain.dim);
    sxz = zeros(settings.float,domain.dim);
    syy = zeros(settings.float,domain.dim);
    syz = zeros(settings.float,domain.dim);
    szz = zeros(settings.float,domain.dim);

    T = typeof(vx)
    fields = Fields{T}(vx, vy, vz, sxx, sxy, sxz, syy, syz, szz)
    return fields 
end