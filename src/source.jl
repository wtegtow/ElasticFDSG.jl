# wavelets 
function get_wavelet(wavelet_type::String, t, ts, fdom)
    if wavelet_type == "ricker"
        return ricker(t, ts, fdom)
    elseif wavelet_type == "gauss1d"
        return gauss1d(t, ts, fdom)
    # register new wavelet types here

    else
        error("Unsupported wavelet type: $wavelet_type")
    end
end

function ricker(t, ts, fdom)  
    wavelet = @. (1 - (2*π^2 * fdom^2) * (t - ts)^2) * (exp(  (-π^2 * fdom^2 ) * (t - ts)^2 ) ) 
    wavelet ./= maximum(abs.(wavelet))
    return wavelet 
end

function gauss1d(t, ts, fdom) 
    wavelet = @. -2 * π^2 * fdom^2 * (t .- ts) * exp.(-π^2 * fdom^2 * (t - ts)^2)
    wavelet ./= maximum(abs.(wavelet))
    return wavelet 
end

# source 
struct Source2D{T, A<:AbstractVector}
    x::T;     z::T
    sx::Int;  sz::Int
    fdom::T;  rhosrc::T
    stf::A;   stf_d1::A
    Mxx::T;   Mxz::T;   Mzz::T
end

struct Source3D{T, A<:AbstractVector}
    x::T;     y::T;     z::T
    sx::Int;  sy::Int;  sz::Int
    fdom::T;  rhosrc::T
    stf::A;   stf_d1::A
    Mxx::T;   Mxy::T;   Mxz::T
    Myy::T;   Myz::T;   Mzz::T
end

function _check_source_in_domain(locs, domain::Domain)
    for (loc, coords) in zip(locs, domain.coordinates)
        if loc < first(coords) || loc > last(coords)
            throw(DomainError(loc, "Source location is outside the defined domain."))
        end
    end
end

function _check_wavelet(ts, fdom)
    if ts < (1.25 / fdom)
        @logger :warn "Wavelet is not fully included in the source-time function (ts < 1.25/fdom). This may cause numerical artefacts."
    end
end

# NOTE: This is dominant wavelength, not minimum wavelength of fmax. So need to set pts_per_lambda a bit larger than recommended 
const MIN_PTS_PER_DOM_LAMBDA_2D = 8
const MIN_PTS_PER_DOM_LAMBDA_3D = 12
function _check_dispersion(domain::Domain, elastic::Elastic, fdom, pts_per_lambda) 
    λ_dom   = elastic.vmin / fdom
    dx_safe = λ_dom / pts_per_lambda
    dx_min  = minimum(map(c -> abs(step(c)), domain.coordinates))
    if dx_min > dx_safe
        @logger :warn "Grid may be too coarse for numerical accuracy. Safe Δh ≤ $(round(dx_safe, digits=2)), current Δh = $dx_min." 
    end
end

function _source_indices(locs, coordinates)
    return map((loc, c) -> argmin(abs.(collect(c) .- loc)), locs, coordinates)
end

function _build_stf(wavelet_type, t, ts, fdom, μ0::T, dt) where T
    stf    = T.(get_wavelet(wavelet_type, t, ts, fdom)) .* μ0
    n      = length(stf)
    stf_d1 = zeros(T, n)
    stf_d1[2:end-1] .= (stf[3:end] .- stf[1:end-2]) ./ (2*dt)
    stf_d1[1]   = (-3*stf[1] + 4*stf[2]   - stf[3])   / (2*dt)
    stf_d1[end] = ( 3*stf[end] - 4*stf[end-1] + stf[end-2]) / (2*dt)
    return stf, stf_d1
end

function init_source(config::Config, domain::Domain{2}, elastic::Elastic, time::SimTime)
    fp   = eval(Symbol(config.dict["settings"]["precision"]))
    scfg = config.dict["source"]

    fdom = fp(scfg["dominant_frequency"])
    ts   = fp(scfg["wavelet_center"])
    x    = fp(scfg["location"]["x"])
    z    = fp(scfg["location"]["z"])

    _check_source_in_domain((x, z), domain)
    _check_wavelet(ts, fdom)
    _check_dispersion(domain, elastic, fdom, MIN_PTS_PER_DOM_LAMBDA_2D)

    sx, sz = _source_indices((x, z), domain.coordinates)

    s      = elastic.c_tensors[elastic.c_lookup[sx, sz]].fields
    rhosrc = fp(s.rho)

    μ0           = fp(scfg["seismic_moment"])
    wavelet_type = scfg["wavelet_type"]
    stf, stf_d1  = _build_stf(wavelet_type, time.t, ts, fdom, μ0, time.dt)
    Mxx = fp(scfg["moment_tensor"]["Mxx"])
    Mxz = fp(scfg["moment_tensor"]["Mxz"])
    Mzz = fp(scfg["moment_tensor"]["Mzz"])

    return Source2D(x, z, sx, sz, fdom, rhosrc, stf, stf_d1, Mxx, Mxz, Mzz)
end

function init_source(config::Config, domain::Domain{3}, elastic::Elastic, time::SimTime)
    fp   = eval(Symbol(config.dict["settings"]["precision"]))
    scfg = config.dict["source"]

    fdom = fp(scfg["dominant_frequency"])
    ts   = fp(scfg["wavelet_center"])
    x    = fp(scfg["location"]["x"])
    y    = fp(scfg["location"]["y"])
    z    = fp(scfg["location"]["z"])

    _check_source_in_domain((x, y, z), domain)
    _check_wavelet(ts, fdom)
    _check_dispersion(domain, elastic, fdom, MIN_PTS_PER_DOM_LAMBDA_3D)

    sx, sy, sz = _source_indices((x, y, z), domain.coordinates)

    s      = elastic.c_tensors[elastic.c_lookup[sx, sy, sz]].fields
    rhosrc = fp(s.rho)

    μ0           = fp(scfg["seismic_moment"])
    wavelet_type = scfg["wavelet_type"]
    stf, stf_d1  = _build_stf(wavelet_type, time.t, ts, fdom, μ0, time.dt)

    Mxx = fp(scfg["moment_tensor"]["Mxx"])
    Mxy = fp(scfg["moment_tensor"]["Mxy"])
    Mxz = fp(scfg["moment_tensor"]["Mxz"])
    Myy = fp(scfg["moment_tensor"]["Myy"])
    Myz = fp(scfg["moment_tensor"]["Myz"])
    Mzz = fp(scfg["moment_tensor"]["Mzz"])

    return Source3D(
        x, y, z, 
        sx, sy, sz, 
        fdom, rhosrc, stf, stf_d1, 
        Mxx, Mxy, Mxz, Myy, Myz, Mzz)
end