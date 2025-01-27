mutable struct Pml{T, U}

    # fields 
    vx_x_old::T
    vy_y_old::T
    vy_x_old::T
    vx_y_old::T
    sxx_x_old::T
    sxy_y_old::T
    sxy_x_old::T
    syy_y_old::T

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

end


function init_pml(settings::Settings,
                  domain::Domain,
                  elastic::Elastic,
                  time::Time,
                  source::Source)

    N = settings.N
    npoints_pml = settings.config["pml"]["nlayer"]
    rcoeff = settings.config["pml"]["reflection_coefficient"]
    RIGHT = settings.config["boundaries"]["right"] == "absorbing" ? true : false
    LEFT  = settings.config["boundaries"]["left"] == "absorbing" ? true : false
    TOP   = settings.config["boundaries"]["top"] == "absorbing" ? true : false
    BOTTOM = settings.config["boundaries"]["bottom"] == "absorbing" ? true : false

    K_x_evn, K_x_odd, 
    a_x_evn, a_x_odd, 
    b_x_evn, b_x_odd = cmpl(N, npoints_pml, RIGHT, LEFT, domain.xcoords, elastic.vpmax, source.fdom, time.dt, rcoeff, settings.float);

    K_y_evn, K_y_odd, 
    a_y_evn, a_y_odd, 
    b_y_evn, b_y_odd = cmpl(N, npoints_pml, TOP, BOTTOM, domain.ycoords, elastic.vpmax, source.fdom, time.dt, rcoeff, settings.float);


    pml_dim = length(domain.pml_points)
    ref_dim = pml_dim 
    

    vx_x_old  = zeros(settings.float,ref_dim)
    vy_y_old  = zeros(settings.float,ref_dim);
    vy_x_old  = zeros(settings.float,ref_dim);
    vx_y_old  = zeros(settings.float,ref_dim);
    sxx_x_old = zeros(settings.float,ref_dim);
    sxy_y_old = zeros(settings.float,ref_dim);
    sxy_x_old = zeros(settings.float,ref_dim);
    syy_y_old = zeros(settings.float,ref_dim);

    # types 
    T = typeof(sxx_x_old)
    U = typeof(K_x_evn)
    
    pml = Pml{T, U}(
              vx_x_old, vy_y_old, vy_x_old, vx_y_old, 
              sxx_x_old, sxy_y_old, sxy_x_old, syy_y_old,
              K_x_evn, K_x_odd, a_x_evn, a_x_odd, b_x_evn, b_x_odd,
              K_y_evn, K_y_odd, a_y_evn, a_y_odd, b_y_evn, b_y_odd)

    return pml
end  
