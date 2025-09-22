mutable struct Pml{T1, T2}
    # fields  
    vx_x_old::T1
    vy_y_old::T1
    vz_z_old::T1
    vy_x_old::T1
    vx_y_old::T1
    vz_x_old::T1
    vx_z_old::T1
    vz_y_old::T1
    vy_z_old::T1
    sxx_x_old::T1
    sxy_y_old::T1
    sxz_z_old::T1
    sxy_x_old::T1
    syy_y_old::T1
    syz_z_old::T1
    sxz_x_old::T1
    syz_y_old::T1
    szz_z_old::T1
    # damping 
    x_evn::T2
    x_odd::T2
    y_evn::T2
    y_odd::T2
    z_evn::T2
    z_odd::T2
end


function init_pml(settings::Settings,
                  domain::Domain,
                  elastic::Elastic,
                  time::Time,
                  source::Source)

    N = settings.N
    npoints_pml = settings.config["boundaries"]["pml_layer"]
    rcoeff = 1e-10
    
    # use pml? 
    xstart = settings.config["boundaries"]["xstart"] == "absorbing" ? true : false
    xend   = settings.config["boundaries"]["xend"]   == "absorbing" ? true : false
    ystart = settings.config["boundaries"]["ystart"] == "absorbing" ? true : false
    yend   = settings.config["boundaries"]["yend"]   == "absorbing" ? true : false
    zstart = settings.config["boundaries"]["zstart"] == "absorbing" ? true : false
    zend   = settings.config["boundaries"]["zend"]   == "absorbing" ? true : false

    K_x_evn, K_x_odd, 
    a_x_evn, a_x_odd, 
    b_x_evn, b_x_odd = cmpl(N, npoints_pml, xstart, xend, domain.xcoords, elastic.vmax, source.fdom, time.dt, rcoeff, settings.float);

    K_y_evn, K_y_odd, 
    a_y_evn, a_y_odd, 
    b_y_evn, b_y_odd = cmpl(N, npoints_pml, ystart, yend, domain.ycoords, elastic.vmax, source.fdom, time.dt, rcoeff, settings.float);

    K_z_evn, K_z_odd, 
    a_z_evn, a_z_odd, 
    b_z_evn, b_z_odd = cmpl(N, npoints_pml, zstart, zend, domain.zcoords, elastic.vmax, source.fdom, time.dt, rcoeff, settings.float);
 
    pml_dim = length(domain.pml_lookup)

    vx_x_old = zeros(settings.float, pml_dim)
    vy_y_old = zeros(settings.float, pml_dim)
    vz_z_old = zeros(settings.float, pml_dim)
    vy_x_old = zeros(settings.float, pml_dim)
    vx_y_old = zeros(settings.float, pml_dim)
    vz_x_old = zeros(settings.float, pml_dim)
    vx_z_old = zeros(settings.float, pml_dim)
    vz_y_old = zeros(settings.float, pml_dim)
    vy_z_old = zeros(settings.float, pml_dim)
    sxx_x_old = zeros(settings.float, pml_dim)
    sxy_y_old = zeros(settings.float, pml_dim)
    sxz_z_old = zeros(settings.float, pml_dim)
    sxy_x_old = zeros(settings.float, pml_dim)
    syy_y_old = zeros(settings.float, pml_dim)
    syz_z_old = zeros(settings.float, pml_dim)
    sxz_x_old = zeros(settings.float, pml_dim)
    syz_y_old = zeros(settings.float, pml_dim)
    szz_z_old = zeros(settings.float, pml_dim)

    # damping (assemble a,b & K in lookup tables)
    x_evn = zeros(settings.float, 3,domain.nx)
    x_evn[1,:] .= a_x_evn
    x_evn[2,:] .= b_x_evn
    x_evn[3,:] .= K_x_evn

    x_odd = zeros(settings.float,3,domain.nx)
    x_odd[1,:] .= a_x_odd
    x_odd[2,:] .= b_x_odd
    x_odd[3,:] .= K_x_odd

    y_evn = zeros(settings.float,3,domain.ny)
    y_evn[1, :] .= a_y_evn
    y_evn[2, :] .= b_y_evn
    y_evn[3, :] .= K_y_evn

    y_odd = zeros(settings.float,3,domain.ny)
    y_odd[1, :] .= a_y_odd
    y_odd[2, :] .= b_y_odd
    y_odd[3, :] .= K_y_odd

    z_evn = zeros(settings.float,3,domain.nz)
    z_evn[1, :] .= a_z_evn
    z_evn[2, :] .= b_z_evn
    z_evn[3, :] .= K_z_evn

    z_odd = zeros(settings.float,3,domain.nz)
    z_odd[1, :] .= a_z_odd
    z_odd[2, :] .= b_z_odd
    z_odd[3, :] .= K_z_odd

    # types 
    T1 = typeof(sxx_x_old)
    T2 = typeof(x_evn)
    
    pml = Pml{T1, T2}(vx_x_old,
                    vy_y_old,
                    vz_z_old,
                    vy_x_old,
                    vx_y_old,
                    vz_x_old,
                    vx_z_old,
                    vz_y_old,
                    vy_z_old,
                    sxx_x_old,
                    sxy_y_old,
                    sxz_z_old,
                    sxy_x_old,
                    syy_y_old,
                    syz_z_old,
                    sxz_x_old,
                    syz_y_old,
                    szz_z_old,
                    x_evn, x_odd,
                    y_evn, y_odd, 
                    z_evn, z_odd)

    return pml
end  
