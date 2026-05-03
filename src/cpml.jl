function cmpl(N, npoints_pml, use_pml_start, use_pml_end, domain_, 
              vmax, fdom, dt, rcoef, FLOAT)

    #= This function contains recycled code from: 
     https://github.com/geodynamics/seismic_cpml
     Roland Martin and Dimitri Komatitsch and Stephen D. Gedney,
     A variational formulation of a stabilized unsplit convolutional perfectly
     matched layer for the isotropic or anisotropic seismic wave equation.
    =#

    # params
    K_max = 1
    α_max = π * fdom  
    NPOW = 3

    n = length(domain_)
    h = abs(domain_[2]-domain_[1])
    small_float = h / 1e3

    # output arrays
    K_evn  = ones(FLOAT,n);
    K_odd  = ones(FLOAT,n);
    a_evn  = zeros(FLOAT,n);
    a_odd  = zeros(FLOAT,n);
    b_evn  = zeros(FLOAT,n);
    b_odd  = zeros(FLOAT,n);
    
    # intermediate arrays
    α_evn = zeros(FLOAT,n)
    α_odd = zeros(FLOAT,n);
    d_evn = zeros(FLOAT,n);
    d_odd = zeros(FLOAT,n);

    # some usefull variables
    pml_width = npoints_pml * h
    ghost_width = N*h
    d0 = -(NPOW-1) * vmax * log(rcoef) / (2 * pml_width)

    # pml region
    pml_start = domain_[begin] + ghost_width : h : domain_[begin] + ghost_width + pml_width - h
    pml_end = domain_[end] - ghost_width - pml_width + h : h : domain_[end] - ghost_width 
   
    # even grid points 
    for i in eachindex(domain_)

        val_i = domain_[i]
        val_in_pml = nothing
        in_pml = false

        if any(isapprox.(val_i, pml_start, atol=small_float)) && use_pml_start == true
            val_in_pml = pml_start[end] - val_i
            in_pml = true

        elseif any(isapprox.(val_i, pml_end, atol=small_float)) && use_pml_end == true
            val_in_pml = val_i - pml_end[1]
            in_pml = true
        end

        if  in_pml == true
            val_in_pml /= pml_width
            d_evn[i] = d0 * val_in_pml^NPOW
            K_evn[i] = 1 + (K_max - 1) * val_in_pml^NPOW
            α_evn[i] = α_max * (1 - val_in_pml)
        
        if α_evn[i] < 0
           α_evn[i] = 0
        end
        end
    end

    # odd grid points
    for i in eachindex(domain_)
        val_i = domain_[i]
        val_in_pml = nothing
        in_pml = false

        if any(isapprox.(val_i, pml_start, atol=small_float)) && use_pml_start == true
            val_in_pml = pml_start[end] - (val_i + h/2)
            in_pml = true

        elseif any(isapprox.(val_i, pml_end, atol=small_float)) && use_pml_end == true
            val_in_pml = (val_i + h/2) - pml_end[1]
            in_pml = true
        end

        if in_pml == true
            val_in_pml /= pml_width 
            d_odd[i] = d0 * val_in_pml^NPOW
            K_odd[i] = 1 + (K_max - 1) * val_in_pml^NPOW
            α_odd[i] = α_max * (1 - val_in_pml)

        if α_odd[i] < 0
            α_odd[i] = 0
        end
        end
    end

    # fill values 
    for i in eachindex(domain_)
        b_evn[i] = exp(- (d_evn[i] / K_evn[i] + α_evn[i]) * dt)
        b_odd[i] = exp(- (d_odd[i] / K_odd[i] + α_odd[i]) * dt)

        if d_evn[i] != 0
            a_evn[i] = d_evn[i] * (b_evn[i] - 1) / (K_evn[i] * (d_evn[i] + K_evn[i] * α_evn[i]))
        end

        if d_odd[i] != 0
            a_odd[i] = d_odd[i] * (b_odd[i] - 1) / (K_odd[i] * (d_odd[i] + K_odd[i] * α_odd[i]))
        end
    end    

    return K_evn, K_odd, a_evn, a_odd, b_evn , b_odd
end


mutable struct CPML2D{T<:AbstractVector, M<:AbstractMatrix}
    # memory variables (sized pml_dim)
    vx_x::T;   vz_z::T
    vz_x::T;   vx_z::T
    sxx_x::T;  sxz_z::T
    sxz_x::T;  szz_z::T
    # damping tables [a; b; K] per axis, evn/odd grids
    x_evn::M;  x_odd::M
    z_evn::M;  z_odd::M
end

mutable struct CPML3D{T<:AbstractVector, M<:AbstractMatrix}
    # memory variables (sized pml_dim)
    vx_x::T;   vy_y::T;   vz_z::T
    vy_x::T;   vx_y::T
    vz_x::T;   vx_z::T
    vz_y::T;   vy_z::T
    sxx_x::T;  sxy_y::T;  sxz_z::T
    sxy_x::T;  syy_y::T;  syz_z::T
    sxz_x::T;  syz_y::T;  szz_z::T
    # damping tables [a; b; K] per axis, evn/odd grids
    x_evn::M;  x_odd::M
    y_evn::M;  y_odd::M
    z_evn::M;  z_odd::M
end

function _cpml_coeff_table(fp, N, npoints_pml, use_s, use_e, coords, vmax, fdom, dt)
    K_evn, K_odd, a_evn, a_odd, b_evn, b_odd = cmpl(
        N, npoints_pml, use_s, use_e, collect(coords), vmax, fdom, dt, fp(1e-10), fp)
    n   = length(coords)
    evn = zeros(fp, 3, n);  evn[1,:] .= a_evn;  evn[2,:] .= b_evn;  evn[3,:] .= K_evn
    odd = zeros(fp, 3, n);  odd[1,:] .= a_odd;  odd[2,:] .= b_odd;  odd[3,:] .= K_odd
    return evn, odd
end


function init_cpml(config::Config, domain::Domain{2}, elastic::Elastic, time::SimTime, source)
    fp          = eval(Symbol(config.dict["settings"]["precision"]))
    N           = Int(config.dict["settings"]["spatial_derivative_order"])
    npoints_pml = Int(config.dict["boundaries"]["pml_layer"])
    bounds      = config.dict["boundaries"]
    pml_dim     = count(domain.pml_lookup .> 0)

    z    = () -> zeros(fp, pml_dim)
    xs, xz = get(bounds,"xstart","none")=="absorbing", get(bounds,"xend","none")=="absorbing"
    zs, ze = get(bounds,"zstart","none")=="absorbing", get(bounds,"zend","none")=="absorbing"

    x_evn, x_odd = _cpml_coeff_table(fp, N, npoints_pml, xs, xz, domain.coordinates[1], elastic.vmax, source.fdom, time.dt)
    z_evn, z_odd = _cpml_coeff_table(fp, N, npoints_pml, zs, ze, domain.coordinates[2], elastic.vmax, source.fdom, time.dt)

    return CPML2D(z(),z(), z(),z(), z(),z(), z(),z(), x_evn, x_odd, z_evn, z_odd)
end


function init_cpml(config::Config, domain::Domain{3}, elastic::Elastic, time::SimTime, source)
    fp          = eval(Symbol(config.dict["settings"]["precision"]))
    N           = Int(config.dict["settings"]["spatial_derivative_order"])
    npoints_pml = Int(config.dict["boundaries"]["pml_layer"])
    bounds      = config.dict["boundaries"]
    pml_dim     = count(domain.pml_lookup .> 0)

    z    = () -> zeros(fp, pml_dim)
    xs, xe = get(bounds,"xstart","none")=="absorbing", get(bounds,"xend","none")=="absorbing"
    ys, ye = get(bounds,"ystart","none")=="absorbing", get(bounds,"yend","none")=="absorbing"
    zs, ze = get(bounds,"zstart","none")=="absorbing", get(bounds,"zend","none")=="absorbing"

    x_evn, x_odd = _cpml_coeff_table(fp, N, npoints_pml, xs, xe, domain.coordinates[1], elastic.vmax, source.fdom, time.dt)
    y_evn, y_odd = _cpml_coeff_table(fp, N, npoints_pml, ys, ye, domain.coordinates[2], elastic.vmax, source.fdom, time.dt)
    z_evn, z_odd = _cpml_coeff_table(fp, N, npoints_pml, zs, ze, domain.coordinates[3], elastic.vmax, source.fdom, time.dt)

    return CPML3D(z(),z(),z(), z(),z(), z(),z(), z(),z(),
                  z(),z(),z(), z(),z(),z(), z(),z(),z(),
                  x_evn, x_odd, y_evn, y_odd, z_evn, z_odd)
end