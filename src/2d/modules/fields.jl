
mutable struct Fields{T}
    vx::T 
    vy::T  
    sxx::T 
    syy::T 
    sxy::T
end


function init_fields(settings::Settings, domain::Domain)
    
    # allocate field quantities
    vx = zeros(settings.float,domain.dim);
    vy = zeros(settings.float,domain.dim);
    sxx = zeros(settings.float,domain.dim);
    syy = zeros(settings.float,domain.dim);
    sxy = zeros(settings.float,domain.dim);

    array_type = typeof(vx)

    fields = Fields{array_type}(vx, vy, sxx, syy, sxy)
    return fields 

end
