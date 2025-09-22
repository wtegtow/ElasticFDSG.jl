mutable struct Fields{T1}
    vx::T1
    vz::T1
    sxx::T1
    sxz::T1
    szz::T1
end

function init_fields(settings::Settings, domain::Domain)
    # allocate primary fields
    vx = zeros(settings.float,domain.dim);
    vz = zeros(settings.float,domain.dim);
    sxx = zeros(settings.float,domain.dim);
    sxz = zeros(settings.float,domain.dim);
    szz = zeros(settings.float,domain.dim);
    # type
    T1 = typeof(vx)
    fields = Fields{T1}(vx, vz, sxx, sxz, szz)
    return fields 
end