mutable struct Pml{T1, T2}
    # fields  
    vx_x_old::T1
    vz_z_old::T1
    vz_x_old::T1
    vx_z_old::T1
    sxx_x_old::T1
    sxz_z_old::T1
    sxz_x_old::T1
    szz_z_old::T1
    # damping 
    x_evn::T2
    x_odd::T2
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
    zstart = settings.config["boundaries"]["zstart"] == "absorbing" ? true : false
    zend   = settings.config["boundaries"]["zend"]   == "absorbing" ? true : false

    K_x_evn, K_x_odd, 
    a_x_evn, a_x_odd, 
    b_x_evn, b_x_odd = cmpl(N, npoints_pml, xstart, xend, domain.xcoords, elastic.vmax, source.fdom, time.dt, rcoeff, settings.float);

    K_z_evn, K_z_odd, 
    a_z_evn, a_z_odd, 
    b_z_evn, b_z_odd = cmpl(N, npoints_pml, zstart, zend, domain.zcoords, elastic.vmax, source.fdom, time.dt, rcoeff, settings.float);
 
    pml_dim = length(domain.pml_lookup)

    vx_x_old = zeros(settings.float, pml_dim)
    vz_z_old = zeros(settings.float, pml_dim)
    vz_x_old = zeros(settings.float, pml_dim)
    vx_z_old = zeros(settings.float, pml_dim)
    sxx_x_old = zeros(settings.float, pml_dim)
    sxz_z_old = zeros(settings.float, pml_dim)
    sxz_x_old = zeros(settings.float, pml_dim)
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
                vz_z_old,
                vz_x_old,
                vx_z_old,
                sxx_x_old,
                sxz_z_old,
                sxz_x_old,
                szz_z_old,
                x_evn, x_odd,
                z_evn, z_odd)
    return pml
end  
