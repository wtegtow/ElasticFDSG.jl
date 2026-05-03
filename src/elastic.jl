struct Stiffness2D
    c11::Real
    c13::Real
    c33::Real
    c44::Real
    rho::Real
end

struct Stiffness3D
    c11::Real
    c12::Real
    c13::Real
    c22::Real
    c23::Real
    c33::Real
    c44::Real
    c55::Real
    c66::Real
    rho::Real
end

struct Stiffness
    fields::Union{Stiffness2D, Stiffness3D}
    dim::Int
end

struct Elastic{N}
    c_tensors::Vector{Stiffness} # unique stiffness tensors in the model
    c_lookup::AbstractArray{Int32, N} # maps each grid point to an index in c_tensors
    vmax::Real
    vmin::Real
end


function _fill_outer_domain!(arr::AbstractArray{<:Any, N}, domain::Domain{N}) where N
    inner_mins = map(ids -> ids[begin], domain.inner_ids)
    inner_maxs = map(ids -> ids[end],   domain.inner_ids)
    for idx in CartesianIndices(domain.shape)
        nearest = CartesianIndex(ntuple(k -> clamp(idx[k], inner_mins[k], inner_maxs[k]), N))
        arr[idx] = arr[nearest]
    end
end

function _init_stiffness(vm::VelocityModel2D, domain::Domain{2}, fp)
    shape = domain.shape
    inner = domain.inner_ids

    vp  = zeros(fp, shape);  vs  = zeros(fp, shape)
    rho = zeros(fp, shape);  eps = zeros(fp, shape);  del = zeros(fp, shape)

    vp[inner...]  .= fp.(vm.vp)
    vs[inner...]  .= fp.(vm.vs)
    rho[inner...] .= fp.(vm.rho)
    eps[inner...] .= fp.(vm.eps)
    del[inner...] .= fp.(vm.del)

    for arr in (vp, vs, rho, eps, del)
        _fill_outer_domain!(arr, domain)
    end

    c33 = @. rho * vp^2
    c44 = @. rho * vs^2
    c11 = @. c33 * (2*eps + 1)
    c13 = @. (c33 - 2*c44) + del * c33

    # Bécache 2003 stability conditions for XZ-plane C-PML
    s1 = @. ((c13+c44)^2 - c11*(c33-c44)) * ((c13+c44)^2 + c44*(c33-c44))
    s2 = @. (c13 + 2*c44)^2 - c11*c33
    s3 = @. (c13+c44)^2 - c11*c33 - c44^2
    if any(s1 .> 0) || any(s2 .> 0) || any(s3 .> 0)
        @warn "Stiffness tensor violates C-PML stability at some grid points." _module=nothing _file=nothing _line=nothing
    end

    cdata    = [Stiffness(Stiffness2D(c11[i,j], c13[i,j], c33[i,j], c44[i,j], rho[i,j]), 2)
                for i in 1:shape[1], j in 1:shape[2]]
    unique_c = unique(vec(cdata))
    c_lud    = Dict(t => Int32(idx) for (idx, t) in enumerate(unique_c))
    c_lookup = [c_lud[cdata[i,j]] for i in 1:shape[1], j in 1:shape[2]]

    vmax = fp(maximum(vp) * 1.1) # simply 10% since phase velocity is function of direction, and i dont want to compute it explicitly
    vmin = fp(any(vs .> 0) ? minimum(vs[vs .> 0]) : minimum(vp) * 0.9)

    return unique_c, c_lookup, vmax, vmin
end


function _init_stiffness(vm::VelocityModel3D, domain::Domain{3}, fp)
    shape = domain.shape
    inner = domain.inner_ids

    vp   = zeros(fp, shape);  vs   = zeros(fp, shape);  rho  = zeros(fp, shape)
    eps1 = zeros(fp, shape);  eps2 = zeros(fp, shape)
    gam1 = zeros(fp, shape);  gam2 = zeros(fp, shape)
    del1 = zeros(fp, shape);  del2 = zeros(fp, shape);  del3 = zeros(fp, shape)

    vp[inner...]   .= fp.(vm.vp);    vs[inner...]   .= fp.(vm.vs)
    rho[inner...]  .= fp.(vm.rho)
    eps1[inner...] .= fp.(vm.eps1);  eps2[inner...] .= fp.(vm.eps2)
    gam1[inner...] .= fp.(vm.gam1);  gam2[inner...] .= fp.(vm.gam2)
    del1[inner...] .= fp.(vm.del1);  del2[inner...] .= fp.(vm.del2);  del3[inner...] .= fp.(vm.del3)

    for arr in (vp, vs, rho, eps1, eps2, gam1, gam2, del1, del2, del3)
        _fill_outer_domain!(arr, domain)
    end

    c33 = @. vp^2 * rho
    c55 = @. vs^2 * rho
    c11 = @. (2*eps2 + 1) * c33
    c22 = @. (2*eps1 + 1) * c33
    c66 = @. (2*gam1 + 1) * c55
    c44 = @. c66 / (1 + 2*gam2)
    c13 = @. sqrt(2*c33*(c33-c55)*del2 + (c33-c55)^2) - c55
    c23 = @. sqrt(2*c33*(c33-c44)*del1 + (c33-c44)^2) - c44
    c12 = @. sqrt(2*c11*(c11-c66)*del3 + (c11-c66)^2) - c66

    # Bécache 2003 stability conditions, all three propagation planes.

    # XZ: P=c11,c33 / cross=c13 / shear=c55
    xz1 = @. ((c13+c55)^2 - c11*(c33-c55)) * ((c13+c55)^2 + c55*(c33-c55))
    xz2 = @. (c13 + 2*c55)^2 - c11*c33
    xz3 = @. (c13+c55)^2 - c11*c33 - c55^2

    # XY: P=c11,c22 / cross=c12 / shear=c66
    xy1 = @. ((c12+c66)^2 - c11*(c22-c66)) * ((c12+c66)^2 + c66*(c22-c66))
    xy2 = @. (c12 + 2*c66)^2 - c11*c22
    xy3 = @. (c12+c66)^2 - c11*c22 - c66^2

    # YZ: P=c22,c33 / cross=c23 / shear=c44
    yz1 = @. ((c23+c44)^2 - c22*(c33-c44)) * ((c23+c44)^2 + c44*(c33-c44))
    yz2 = @. (c23 + 2*c44)^2 - c22*c33
    yz3 = @. (c23+c44)^2 - c22*c33 - c44^2

    if any(xz1.>0)||any(xz2.>0)||any(xz3.>0) ||
       any(xy1.>0)||any(xy2.>0)||any(xy3.>0) ||
       any(yz1.>0)||any(yz2.>0)||any(yz3.>0)
        @warn "Stiffness tensor violates C-PML stability at some grid points." _module=nothing _file=nothing _line=nothing
    end

    cdata = [Stiffness(Stiffness3D(c11[i,j,k], c12[i,j,k], c13[i,j,k],
                                   c22[i,j,k], c23[i,j,k], c33[i,j,k],
                                   c44[i,j,k], c55[i,j,k], c66[i,j,k], rho[i,j,k]), 3)
                for i in 1:shape[1], j in 1:shape[2], k in 1:shape[3]]
    unique_c = unique(vec(cdata))
    c_lud    = Dict(t => Int32(idx) for (idx, t) in enumerate(unique_c))
    c_lookup = [c_lud[cdata[i,j,k]] for i in 1:shape[1], j in 1:shape[2], k in 1:shape[3]]

    vmax = fp(maximum(vp) * 1.1) # simply 10% since phase velocity is function of direction, and i dont want to compute it explicitly
    vmin = fp(any(vs .> 0) ? minimum(vs[vs .> 0]) : minimum(vp) * 0.9)

    return unique_c, c_lookup, vmax, vmin
end

function init_elastic(config::Config, domain::Domain, velmod::VelocityModel)
    fp = eval(Symbol(config.dict["settings"]["precision"]))
    c_tensors, c_lookup, vmax, vmin = _init_stiffness(velmod.fields, domain, fp)
    return Elastic(c_tensors, c_lookup, vmax, vmin)
end;