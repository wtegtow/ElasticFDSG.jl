mutable struct Pml{T, U}

    # fields 
    vx_x_old::T
    vy_y_old::T
    vz_z_old::T
    vy_x_old::T
    vx_y_old::T
    vz_x_old::T
    vx_z_old::T
    vz_y_old::T
    vy_z_old::T
    sxx_x_old::T
    sxy_y_old::T
    sxz_z_old::T
    sxy_x_old::T
    syy_y_old::T
    syz_z_old::T
    sxz_x_old::T
    syz_y_old::T
    szz_z_old::T

    # damping 
    K_x_evn::U 
    K_x_odd::U
    a_x_evn::U
    a_x_odd::U
    b_x_evn::U
    b_x_odd::U

    K_y_evn::U
    K_y_odd::U
    a_y_evn::U
    a_y_odd::U
    b_y_evn::U
    b_y_odd::U

    K_z_evn::U
    K_z_odd::U
    a_z_evn::U
    a_z_odd::U
    b_z_evn::U
    b_z_odd::U

end


function init_pml(settings::Settings,
                  domain::Domain,
                  elastic::Elastic,
                  time::Time,
                  source::Source)

    N = settings.N
    npoints_pml = settings.config["pml"]["nlayer"]
    rcoeff = settings.config["pml"]["reflection_coefficient"]
    
    # use pml? 
    xstart = settings.config["boundaries"]["xstart"] == "absorbing" ? true : false
    xend   = settings.config["boundaries"]["xend"]   == "absorbing" ? true : false
    ystart = settings.config["boundaries"]["ystart"] == "absorbing" ? true : false
    yend   = settings.config["boundaries"]["yend"]   == "absorbing" ? true : false
    zstart = settings.config["boundaries"]["zstart"] == "absorbing" ? true : false
    zend   = settings.config["boundaries"]["zend"]   == "absorbing" ? true : false

    K_x_evn, K_x_odd, 
    a_x_evn, a_x_odd, 
    b_x_evn, b_x_odd = cmpl(N, npoints_pml, xstart, xend, domain.xcoords, elastic.vpmax, source.fdom, time.dt, rcoeff, settings.float);

    K_y_evn, K_y_odd, 
    a_y_evn, a_y_odd, 
    b_y_evn, b_y_odd = cmpl(N, npoints_pml, ystart, yend, domain.ycoords, elastic.vpmax, source.fdom, time.dt, rcoeff, settings.float);

    K_z_evn, K_z_odd, 
    a_z_evn, a_z_odd, 
    b_z_evn, b_z_odd = cmpl(N, npoints_pml, zstart, zend, domain.zcoords, elastic.vpmax, source.fdom, time.dt, rcoeff, settings.float);

    pml_dim = length(domain.pml_points)
    ref_dim = pml_dim 
    
    vx_x_old = zeros(settings.float,ref_dim);
    vy_y_old = zeros(settings.float,ref_dim);
    vz_z_old = zeros(settings.float,ref_dim);
    vy_x_old = zeros(settings.float,ref_dim);
    vx_y_old = zeros(settings.float,ref_dim);
    vz_x_old = zeros(settings.float,ref_dim);
    vx_z_old = zeros(settings.float,ref_dim);
    vz_y_old = zeros(settings.float,ref_dim);
    vy_z_old = zeros(settings.float,ref_dim);
    sxx_x_old = zeros(settings.float,ref_dim);
    sxy_y_old = zeros(settings.float,ref_dim);
    sxz_z_old = zeros(settings.float,ref_dim);
    sxy_x_old = zeros(settings.float,ref_dim);
    syy_y_old = zeros(settings.float,ref_dim);
    syz_z_old = zeros(settings.float,ref_dim);
    sxz_x_old = zeros(settings.float,ref_dim);
    syz_y_old = zeros(settings.float,ref_dim);
    szz_z_old = zeros(settings.float,ref_dim);

    # types 
    T = typeof(sxx_x_old)
    U = typeof(K_x_evn)
    
    pml = Pml{T, U}(vx_x_old,
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
                    K_x_evn, K_x_odd, 
                    a_x_evn, a_x_odd, 
                    b_x_evn, b_x_odd,
                    K_y_evn, K_y_odd, 
                    a_y_evn, a_y_odd, 
                    b_y_evn, b_y_odd,
                    K_z_evn, K_z_odd, 
                    a_z_evn, a_z_odd, 
                    b_z_evn, b_z_odd)

    return pml
end  
