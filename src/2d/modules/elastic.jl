struct Elastic{T1,T2,T3}
    c_tensors::T1     
    c_lookup::T2              
    vmax::T3  
    vmin::T3 
end

function fill_outer_domain!(domain::Domain, array)

    nz = domain.nz
    nx = domain.nx 
    inner_id = domain.inner_id

    z_min, z_max = inner_id[1][begin], inner_id[1][end]
    x_min, x_max = inner_id[2][begin], inner_id[2][end]

    for z in 1:nz 
        for x in 1:nx 
            if z_min <= z <= z_max && 
               x_min <= x <= x_max
                continue  
            else

            nearest_z = clamp(z, z_min, z_max)
            nearest_x = clamp(x, x_min, x_max)

            array[z,x] = array[nearest_z, nearest_x]

            end
        end 
    end
end;

# Christoffel matrix Γ(n) = n C n 
function Γn(c11, c13, c33, c55, n)
    Γ = @SMatrix[
         c11 * n[1]^2 + c55 * n[2]^2  (c13 + c55) * n[1] * n[2]    ;
        (c13 + c55) * n[1] * n[2]      c55 * n[1]^2 + c33 * n[2]^2
        ]
    return Γ
end

function solve_christoffel!(VpVs, UpUs, c11, c13, c33, c55, n) 
    Γ = Γn(c11, c13, c33, c55, n) 
    F = eigen(Γ)           
    V = F.values
    U = F.vectors
    # P -> 1, S1 -> 2
    VpVs[1] = sqrt(real(V[2]))  
    VpVs[2] = sqrt(real(V[1]))  
    UpUs[:,1] .= U[:,2]         
    UpUs[:,2] .= U[:,1]       
end

function group_velocity(VpVs, UpUs, c11, c13, c33, c55, n)

    # solve for slowness vector P
    solve_christoffel!(VpVs, UpUs, c11, c13, c33, c55, n) 
    Pp  = n ./ VpVs[1]
    Ps1 = n ./ VpVs[2]

    # solve for group velocity vector
    ΓUp  = Γn(c11, c13, c33, c55, UpUs[:,1])
    ΓUs1 = Γn(c11, c13, c33, c55, UpUs[:,2])
    @einsum  gp[i] := ΓUp[i,j] * Pp[j]
    @einsum gs1[i] := ΓUs1[i,j] * Ps1[j]
    @assert isapprox(dot(gp, Pp), 1; rtol=1e-3)
    @assert isapprox(dot(gs1, Ps1), 1; rtol=1e-3)
    return (gp, gs1)
end;


function init_elastic(settings::Settings, domain::Domain, velmod)

    dim = domain.dim
    inner_id = domain.inner_id

    vp = zeros(eltype(velmod),dim);
    vs = zeros(eltype(velmod),dim);
    rho = zeros(eltype(velmod),dim);
    eps = zeros(eltype(velmod),dim);
    del = zeros(eltype(velmod),dim);

    # fill inner domain with velmod 
    vp[inner_id...]   .= permutedims(velmod[3,:,:], (2,1));
    vs[inner_id...]   .= permutedims(velmod[4,:,:], (2,1));
    rho[inner_id...]  .= permutedims(velmod[5,:,:], (2,1));
    eps[inner_id...]  .= permutedims(velmod[6,:,:], (2,1));
    del[inner_id...]  .= permutedims(velmod[7,:,:], (2,1));

    # fill outer domain of nearest values of inner domain
    fill_outer_domain!(domain,vp);
    fill_outer_domain!(domain,vs);
    fill_outer_domain!(domain,rho);
    fill_outer_domain!(domain,eps);
    fill_outer_domain!(domain,del);

    # stiffness 
    lam = similar(vp);
    mu = similar(vp);
    lam[:,:] .= rho .* (vp.^2 - 2 .* vs.^2);
    mu[:,:]  .= rho .* vs.^2;

    c11 = (lam + 2 .* mu) .* (2*eps .+ 1)
    c13 = lam + del .* (lam + 2*mu)
    c33 = lam + 2*mu
    c55 = mu 

    # find unique c-tensors and create corresponding index lookup table
    cdata = [(c11[z,x], c13[z,x], c33[z,x], c55[z,x],rho[z,x])
            for z in 1:dim[1], x in 1:dim[2]]

    unique_c = unique(vec(cdata))
    c_lud = Dict(t => idx for (idx,t) in enumerate(unique_c))

    c_lookup = zeros(Int32, dim)  # nz,ny,nx
    for i in 1:dim[1], j in 1:dim[2]
        c_lookup[i,j] = c_lud[cdata[i,j]]
    end

    c_tensors = zeros(settings.float, length(unique_c), length(unique_c[1]))
    for (idx, tup) in enumerate(unique_c)
        c_tensors[idx, :] .= collect(tup)
    end

    # calculate fastest and slowest group velocity 
    if any(vp .== 0) == true  
        @warn "P-wave velocity somewhere zero. This may cause numerical issues."
    end

    angles_theta = deg2rad.(0:2:360)
    VpVs = MVector{2,Float64}(undef)
    UpUs = MMatrix{2,2,Float64}(undef)

    vmax = 0.0
    vmin = Inf
    for c in eachrow(deepcopy(c_tensors))
        c ./= c[end] # density normalized stiffness
        c11, c13, c33, c55, rho_ = c 

        if c55 == 0  # skip liquid layer
            continue 
        end

        for theta in angles_theta
            
            n = SVector(cos(theta), sin(theta))
            n = n/norm(n)
            gp, gs1 = group_velocity(VpVs, UpUs, c11, c13, c33, c55, n)
            comb = [norm(gp), norm(gs1)]
            vmax_ = maximum(comb)
            vmin_ = minimum(comb)
            if vmax_ > vmax 
                vmax = vmax_
            end
            if vmin_ < vmin 
                vmin = vmin_
            end
             
        end
    end;

    # fallback (if all nodes are liquids)
    if vmax == 0.0 
        vmax = maximum([maximum(vp), maximum(vs)])
    end;
    if vmin == Inf 
        vmin = minimum(vp)
    end;

    # types 
    T1 = typeof(c_tensors)
    T2 = typeof(c_lookup)
    T3 = typeof(vmin)
    elastic = Elastic{T1,T2,T3}(c_tensors, c_lookup, vmax, vmin)

    return elastic 
end
