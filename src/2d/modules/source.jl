struct Source{T1,T2,T3}
    x::T1
    z::T1
    sx::Int64
    sz::Int64

    fdom::T1
    rhosrc::T1
    stf::T2
    stf_d1::T2

    Mxx::T3
    Mxz::T3
    Mzz::T3
end

function init_source(settings::Settings,
                    domain::Domain,
                    elastic::Elastic,
                    time::Time)

    x0 = domain.xcoords[begin]
    xend = domain.xcoords[end]
    z0 = domain.zcoords[begin]
    zend = domain.zcoords[end]

    fdom = settings.float(settings.config["source"]["dominant_frequency"])
    ts = settings.float(settings.config["source"]["wavelet_center"])

    source_x = settings.float(settings.config["source"]["location"]["x"])
    source_z = settings.float(settings.config["source"]["location"]["z"])

    # is source within the domain?
    if source_x >= x0   && 
       source_x <= xend && 
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
    dx_min = min(domain.dx, domain.dz)

    if !(dx_min ≤ dx_safe)
        @warn """Grid might be too coarse and the solution may be affected by numerical dispersion.
        At least 6 grid points per slowest S-wave wavelength are recommended (safe Δh ≤ $(round(dx_safe,digits=2)), current Δh = $dx_min).""" _module=nothing _file=nothing _line=nothing
    end

    # indices of source
    sx = argmin(abs.(domain.xcoords .- source_x))
    sz = argmin(abs.(domain.zcoords .- source_z))

    # density at source location
    rhosrc = elastic.c_tensors[elastic.c_lookup[sz,sx],5]
    
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
    Mxz = settings.float(settings.config["source"]["moment_tensor"]["Mxz"])
    Mzz = settings.float(settings.config["source"]["moment_tensor"]["Mzz"])

    # if anisotropic
    if settings.config["source"]["moment_tensor"]["anisotropic"]

        c11, c13, c33, c55, rho = elastic.c_tensors[elastic.c_lookup[sz,sx],:]

        C_source = [c11 c13 0; 
                    c13 c33 0;
                     0   0 c55]
        M = [Mxx Mzz Mxz] * C_source 
        Mxx = M[1]
        Mzz = M[2]
        Mxz = M[3]
    end;

    # normalize 
    M = [Mxx Mxz;
         Mxz Mzz] 
    M = normalize(M)
    Mxx = M[1,1]
    Mzz = M[2,2]
    Mxz = M[1,2]

    # types 
    T1 = typeof(source_x)
    T2 = typeof(stf)
    T3 = typeof(Mxx)
    source = Source{T1, T2, T3}(source_x, source_z, 
                            sx, sz, 
                            fdom, rhosrc, stf, stf_d1,
                            Mxx, Mzz, Mxz)

    return source 

end;
