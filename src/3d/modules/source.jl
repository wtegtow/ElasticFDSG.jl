struct Source{T1,T2,T3}
    x::T1
    y::T1
    z::T1
    sx::Int64
    sy::Int64
    sz::Int64

    fdom::T1
    rhosrc::T1
    stf::T2
    stf_d1::T2

    Mxx::T3
    Mxy::T3
    Mxz::T3
    Myy::T3
    Myz::T3
    Mzz::T3
end

function init_source(settings::Settings,
                    domain::Domain,
                    elastic::Elastic,
                    time::Time)

    x0 = domain.xcoords[begin]
    xend = domain.xcoords[end]
    y0 = domain.ycoords[begin]
    yend = domain.ycoords[end]
    z0 = domain.zcoords[begin]
    zend = domain.zcoords[end]

    fdom = settings.float(settings.config["source"]["dominant_frequency"])
    ts = settings.float(settings.config["source"]["wavelet_center"])

    source_x = settings.float(settings.config["source"]["location"]["x"])
    source_y = settings.float(settings.config["source"]["location"]["y"])
    source_z = settings.float(settings.config["source"]["location"]["z"])

    # is source within the domain?
    if source_x >= x0   && 
       source_x <= xend && 
       source_y >= y0   && 
       source_y <= yend &&
       source_z >= z0   && 
       source_z <= zend 
    ;
    else
        throw(DimensionMismatch("Source location is not within the defined domain."))
    end

    # is the wavelet truncated?
    if ts < (1.25 / fdom)
        @warn """Wavelet is not fully included in the source-time function (ts !< 1.25/fdom).
        This may cause oscillatory behavior in the medium's response."""  _module=nothing _file=nothing _line=nothing
    end

    # is the grid spacing sufficient fine ? 
    λ_dom = elastic.vmin / fdom
    dx_safe = λ_dom / 6
    dx_min = min(domain.dx, domain.dy, domain.dz)

    if !(dx_min ≤ dx_safe)
        @warn """Grid might be too coarse and the solution may be affected by numerical dispersion.
        At least 6 grid points per slowest S-wave wavelength are recommended (safe Δh ≤ $(round(dx_safe,digits=2)), current Δh = $dx_min).""" _module=nothing _file=nothing _line=nothing
    end

    # indices of source
    sx = argmin(abs.(domain.xcoords .- source_x))
    sy = argmin(abs.(domain.ycoords .- source_y))
    sz = argmin(abs.(domain.zcoords .- source_z))

    # density at source location
    rhosrc = elastic.c_tensors[elastic.c_lookup[sz,sy,sx],10]
    
    # source time function 
    wavelet_type = settings.config["source"]["wavelet_type"]
    if wavelet_type == "ricker"
        stf = ricker(time.t, ts, fdom)
    elseif wavelet_type == "gauss1d"
        stf = gauss1d(time.t, ts, fdom)
    end; 
    μ0  = settings.float(settings.config["source"]["seismic_moment"])
    stf = settings.float.(stf) .* μ0
    stf_d1 = zeros(eltype(stf), size(stf,1))
    stf_d1[begin:end-1] .= diff(stf) ./ time.dt

    # moment tensor
    Mxx = settings.float(settings.config["source"]["moment_tensor"]["Mxx"])
    Mxy = settings.float(settings.config["source"]["moment_tensor"]["Mxy"])
    Mxz = settings.float(settings.config["source"]["moment_tensor"]["Mxz"])
    Myy = settings.float(settings.config["source"]["moment_tensor"]["Myy"])
    Myz = settings.float(settings.config["source"]["moment_tensor"]["Myz"])
    Mzz = settings.float(settings.config["source"]["moment_tensor"]["Mzz"])

    # if anisotropic
    if settings.config["source"]["moment_tensor"]["anisotropic"]

        c11, c12, c13, c22, c23, c33, c44, c55, c66, rho = elastic.c_tensors[elastic.c_lookup[sz,sy,sx],:]

        C_source = [c11 c12 c13 0 0 0; 
                    c12 c22 c23 0 0 0;
                    c13 c23 c33 0 0 0;
                    0 0 0      c44 0 0;
                    0 0 0      0 c55 0;
                    0 0 0      0 0 c66]
        M = [Mxx Myy Mzz Myz Mxz Mxy] * C_source 
        Mxx = M[1]
        Myy = M[2]
        Mzz = M[3]
        Myz = M[4]
        Mxz = M[5]
        Mxy = M[6]
    end;

    # normalize 
    M = [Mxx Mxy Mxz;
         Mxy Myy Myz;
         Mxz Myz Mzz] 
    M = normalize(M)
    Mxx = M[1,1]
    Myy = M[2,2]
    Mzz = M[3,3]
    Myz = M[2,3]
    Mxz = M[1,3]
    Mxy = M[1,2]

    # types 
    T1 = typeof(source_x)
    T2 = typeof(stf)
    T3 = typeof(Mxx)
    source = Source{T1, T2, T3}(source_x, source_y, source_z, 
                            sx, sy, sz, 
                            fdom, rhosrc, stf,stf_d1,
                            Mxx, Mxy, Mxz, Myy, Myz, Mzz)

    return source 
end;
