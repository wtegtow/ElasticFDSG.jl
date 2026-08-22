struct Stiffness2D{T<:AbstractFloat}
    c11::T
    c13::T
    c33::T
    c44::T
    rho::T
end

struct Stiffness3D{T<:AbstractFloat}
    c11::T
    c12::T
    c13::T
    c22::T
    c23::T
    c33::T
    c44::T
    c55::T
    c66::T
    rho::T
end

struct Stiffness
    fields::Union{Stiffness2D, Stiffness3D}
    dim::Int
end

struct Elastic{N}
    c_tensors::Vector{Stiffness}    # unique stiffness tensors in the model
    c_lookup::AbstractArray{Int, N} # maps each grid point to an index in c_tensors
    vmax::AbstractFloat
    vmin::AbstractFloat
end

function _fill_outer_domain!(arr::AbstractArray{<:Any, N}, domain::Domain{N}) where N
    inner_mins = map(ids -> ids[begin], domain.inner_ids)
    inner_maxs = map(ids -> ids[end],   domain.inner_ids)
    for idx in CartesianIndices(domain.shape)
        nearest = CartesianIndex(ntuple(k -> clamp(idx[k], inner_mins[k], inner_maxs[k]), N))
        arr[idx] = arr[nearest]
    end
end

function _velocity_bounds(vp, vs, ::Type{T}) where {T}
    
    # Maximum for CFL stability
    vmax = T(maximum(vp) * T(1.1)) # add 10% since phase velocity is function of direction. I dont want to compute the max here explicitly.

    # Minimum relevant wave speed for dispersion
    tol = T(1e-4) # tolerance for none solids
    candidates = vcat(vp[vp .> tol], vs[vs .> tol])

    isempty(candidates) && throw(ArgumentError("No valid wave velocities found."))
    vmin = T(minimum(candidates) * T(0.9)) # remove 10%, same as above

    @logger :debug "vmax/vmin: $(vmax)/$(vmin)"
    return vmax, vmin
end

function _check_cpml_stability(unique_c::Vector{Stiffness2D{T}}) where {T<:AbstractFloat}

    # Bécache 2003 stability conditions for (X,Z)-plane
    @inbounds for s in unique_c
        c11_ = Float64(s.c11) # f64 required here to prevent "catastrophic cancellation".
        c13_ = Float64(s.c13)
        c33_ = Float64(s.c33)
        c44_ = Float64(s.c44)

        s1 = ((c13_ + c44_)^2 - c11_ * (c33_ - c44_)) *
             ((c13_ + c44_)^2 + c44_ * (c33_ - c44_))

        s2 = (c13_ + 2*c44_)^2 - c11_ * c33_
        s3 = (c13_ + c44_)^2 - c11_ * c33_ - c44_^2

        if s1 > 0 || s2 > 0 || s3 > 0
            @logger :warn "Stiffness tensor violates C-PML stability at some grid points."
            return 
        end
    end
    @logger :debug "Stiffness tensor does not violates C-PML stability"
    return 
end

function _init_stiffness(vm::VelocityModel2D, domain::Domain{2}, fp)
    shape = domain.shape
    inner = domain.inner_ids

    vp  = zeros(fp, shape);  vp[inner...]  .= vm.vp
    vs  = zeros(fp, shape);  vs[inner...]  .= vm.vs
    rho = zeros(fp, shape);  rho[inner...] .= vm.rho
    eps = zeros(fp, shape);  eps[inner...] .= vm.eps
    del = zeros(fp, shape);  del[inner...] .= vm.del

    for arr in (vp, vs, rho, eps, del)
        _fill_outer_domain!(arr, domain)
    end
    vmax, vmin = _velocity_bounds(vp, vs, fp)

    c33 = @. rho * vp^2
    c44 = @. rho * vs^2
    c11 = @. c33 * (2*eps + 1)
    c13 = @. (c33 - 2*c44) + del * c33
   
    # determine and check unique c-tensors
    cdata    = [Stiffness2D(c11[i,j], c13[i,j], c33[i,j], c44[i,j], rho[i,j]) for i in 1:shape[1], j in 1:shape[2]]
    unique_c = unique(vec(cdata)) 
    _check_cpml_stability(unique_c)

    # index mapping 
    c_lud    = Dict(t => idx for (idx, t) in enumerate(unique_c))
    c_lookup = [c_lud[cdata[i,j]] for i in 1:shape[1], j in 1:shape[2]]

    return Stiffness.(unique_c,2), c_lookup, vmax, vmin
end


function _check_cpml_stability(unique_c::Vector{Stiffness3D{T}}) where {T<:AbstractFloat}
    # Bécache 2003 stability conditions, all three propagation planes.
    @inbounds for s in unique_c
        c11_ = Float64(s.c11)  # f64 required here to prevent "catastrophic cancellation".
        c12_ = Float64(s.c12)
        c13_ = Float64(s.c13)
        c22_ = Float64(s.c22)
        c23_ = Float64(s.c23)
        c33_ = Float64(s.c33)
        c44_ = Float64(s.c44)
        c55_ = Float64(s.c55)
        c66_ = Float64(s.c66)

        # XZ plane
        xz1 = ((c13_ + c55_)^2 - c11_ * (c33_ - c55_)) *
              ((c13_ + c55_)^2 + c55_ * (c33_ - c55_))

        xz2 = (c13_ + 2*c55_)^2 - c11_ * c33_
        xz3 = (c13_ + c55_)^2 - c11_ * c33_ - c55_^2

        # XY plane
        xy1 = ((c12_ + c66_)^2 - c11_ * (c22_ - c66_)) *
              ((c12_ + c66_)^2 + c66_ * (c22_ - c66_))

        xy2 = (c12_ + 2*c66_)^2 - c11_ * c22_
        xy3 = (c12_ + c66_)^2 - c11_ * c22_ - c66_^2

        # YZ plane
        yz1 = ((c23_ + c44_)^2 - c22_ * (c33_ - c44_)) *
              ((c23_ + c44_)^2 + c44_ * (c33_ - c44_))

        yz2 = (c23_ + 2*c44_)^2 - c22_ * c33_
        yz3 = (c23_ + c44_)^2 - c22_ * c33_ - c44_^2

        if xz1 > 0 || xz2 > 0 || xz3 > 0 ||
           xy1 > 0 || xy2 > 0 || xy3 > 0 ||
           yz1 > 0 || yz2 > 0 || yz3 > 0

            @logger :warn "Stiffness tensor violates C-PML stability at some grid points."
            return 
        end
    end
    @logger :debug "Stiffness tensor does not violates C-PML stability"
    return 
end

function _init_stiffness(vm::VelocityModel3D, domain::Domain{3}, fp)
    shape = domain.shape
    inner = domain.inner_ids

    # NOTE: This function temporarly makes many copys of arrays... 
    # But this is fine here, since later the application anyway needs to allocate even more memeory (fields). 
    vp   = zeros(fp, shape);  vp[inner...]   .= vm.vp
    vs   = zeros(fp, shape);  vs[inner...]   .= vm.vs
    rho  = zeros(fp, shape);  rho[inner...]  .= vm.rho
    eps1 = zeros(fp, shape);  eps1[inner...] .= vm.eps1
    eps2 = zeros(fp, shape);  eps2[inner...] .= vm.eps2
    gam1 = zeros(fp, shape);  gam1[inner...] .= vm.gam1
    gam2 = zeros(fp, shape);  gam2[inner...] .= vm.gam2
    del1 = zeros(fp, shape);  del1[inner...] .= vm.del1
    del2 = zeros(fp, shape);  del2[inner...] .= vm.del2
    del3 = zeros(fp, shape);  del3[inner...] .= vm.del3

    for arr in (vp, vs, rho, eps1, eps2, gam1, gam2, del1, del2, del3)
        _fill_outer_domain!(arr, domain)
    end
    vmax, vmin = _velocity_bounds(vp, vs, fp)

    c33 = @. vp^2 * rho
    c55 = @. vs^2 * rho
    c11 = @. (2*eps2 + 1) * c33
    c22 = @. (2*eps1 + 1) * c33
    c66 = @. (2*gam1 + 1) * c55
    c44 = @. c66 / (1 + 2*gam2)
    c13 = @. sqrt(2*c33*(c33-c55)*del2 + (c33-c55)^2) - c55
    c23 = @. sqrt(2*c33*(c33-c44)*del1 + (c33-c44)^2) - c44
    c12 = @. sqrt(2*c11*(c11-c66)*del3 + (c11-c66)^2) - c66

    # determine and check unique c-tensors
    cdata = [Stiffness3D(c11[i,j,k], c12[i,j,k], c13[i,j,k],
                         c22[i,j,k], c23[i,j,k], c33[i,j,k],
                         c44[i,j,k], c55[i,j,k], c66[i,j,k], rho[i,j,k]) 
            for i in 1:shape[1], j in 1:shape[2], k in 1:shape[3]]
    unique_c = unique(vec(cdata))
    _check_cpml_stability(unique_c)
    
    # index mapping 
    c_lud    = Dict(t => idx for (idx, t) in enumerate(unique_c))
    c_lookup = [c_lud[cdata[i,j,k]] for i in 1:shape[1], j in 1:shape[2], k in 1:shape[3]]

    return Stiffness.(unique_c,3), c_lookup, vmax, vmin
end

function init_elastic(config::Config, domain::Domain, velmod::VelocityModel)
    fp = eval(Symbol(config.dict["settings"]["precision"]))
    c_tensors, c_lookup, vmax, vmin = _init_stiffness(velmod.fields, domain, fp)
    elastic = Elastic(c_tensors, c_lookup, vmax, vmin)
    @logger :info "Elastic initialized"
    return elastic
end;