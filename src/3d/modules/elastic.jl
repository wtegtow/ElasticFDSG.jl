struct Elastic{T1,T2,T3}
    c_tensors::T1     
    c_lookup::T2              
    vmax::T3  
    vmin::T3 
end

function fill_outer_domain!(domain::Domain, array)

    nz = domain.nz
    ny = domain.ny 
    nx = domain.nx 
    inner_id = domain.inner_id

    z_min, z_max = inner_id[1][begin], inner_id[1][end]
    y_min, y_max = inner_id[2][begin], inner_id[2][end]
    x_min, x_max = inner_id[3][begin], inner_id[3][end]

    for z in 1:nz 
        for y in 1:ny 
            for x in 1:nx 
                if z_min <= z <= z_max && 
                   y_min <= y <= y_max && 
                   x_min <= x <= x_max
                    continue  
                else

                nearest_z = clamp(z, z_min, z_max)
                nearest_y = clamp(y, y_min, y_max)
                nearest_x = clamp(x, x_min, x_max)

                array[z,y,x] = array[nearest_z,nearest_y, nearest_x]

                end
            end 
        end 
    end
end;

# Christoffel matrix Γ(n) = n C n 
function Γn(c11, c12, c13, c22, c23, c33, c44, c55, c66, n)
    n1, n2, n3 = n[1], n[2], n[3]
    Γ = @SMatrix[
         c11 * n1^2 + c66 * n2^2 + c55 * n3^2  (c12 + c66) * n1 * n2                   (c13 + c55) * n1 * n3;
        (c12 + c66) * n1 * n2                   c66 * n1^2 + c22 * n2^2 + c44 * n3^2   (c23 + c44) * n2 * n3;
        (c13 + c55) * n1 * n3                  (c23 + c44) * n2 * n3                    c55 * n1^2 + c44 * n2^2 + c33 * n3^2
    ]
    return Γ
end

function solve_christoffel!(VpVs, UpUs, c11, c12, c13, c22, c23, c33, c44, c55, c66, n) 
    Γ = Γn(c11, c12, c13, c22, c23, c33, c44, c55, c66, n) 
    F = eigen(Γ)           
    V = F.values
    U = F.vectors
    # P -> 1, S1 -> 2, S2 -> 3
    VpVs[1] = sqrt(real(V[3]))  
    VpVs[2] = sqrt(real(V[2]))  
    VpVs[3] = sqrt(real(V[1])) 
    UpUs[:,1] .= U[:,3]         
    UpUs[:,2] .= U[:,2]       
    UpUs[:,3] .= U[:,1]   
end

function group_velocity(VpVs, UpUs, c11, c12, c13, c22, c23, c33, c44, c55, c66, n)

    # solve for slowness vector P
    solve_christoffel!(VpVs, UpUs, c11, c12, c13, c22, c23, c33, c44, c55, c66, n) 
    Pp  = n ./ VpVs[1]
    Ps1 = n ./ VpVs[2]
    Ps2 = n ./ VpVs[3]

    # solve for group velocity vector
    ΓUp  = Γn(c11, c12, c13, c22, c23, c33, c44, c55, c66, UpUs[:,1])
    ΓUs1 = Γn(c11, c12, c13, c22, c23, c33, c44, c55, c66, UpUs[:,2])
    ΓUs2 = Γn(c11, c12, c13, c22, c23, c33, c44, c55, c66, UpUs[:,3])
    @einsum  gp[i] := ΓUp[i,j] * Pp[j]
    @einsum gs1[i] := ΓUs1[i,j] * Ps1[j]
    @einsum gs2[i] := ΓUs2[i,j] * Ps2[j]
    @assert isapprox(dot(gp, Pp), 1; rtol=1e-3)
    @assert isapprox(dot(gs1, Ps1), 1; rtol=1e-3)
    @assert isapprox(dot(gs2, Ps2), 1; rtol=1e-3)
    return (gp, gs1, gs2)
end;


function init_elastic(settings::Settings, domain::Domain, velmod)

    dim = domain.dim
    inner_id = domain.inner_id

    vp = zeros(eltype(velmod),dim);
    vs = zeros(eltype(velmod),dim);
    rho = zeros(eltype(velmod),dim);
    eps1 = zeros(eltype(velmod),dim);
    eps2 = zeros(eltype(velmod),dim);
    gam1 = zeros(eltype(velmod),dim);
    gam2 = zeros(eltype(velmod),dim);
    del1 = zeros(eltype(velmod),dim);
    del2 = zeros(eltype(velmod),dim);
    del3 = zeros(eltype(velmod),dim);

    # fill inner domain with velmod 
    vp[inner_id...]   .= permutedims(velmod[4, :, :, :], (3, 2, 1));
    vs[inner_id...]   .= permutedims(velmod[5, :, :, :], (3, 2, 1));
    rho[inner_id...]  .= permutedims(velmod[6, :, :, :], (3, 2, 1));
    eps1[inner_id...] .= permutedims(velmod[7, :, :, :], (3, 2, 1));
    eps2[inner_id...] .= permutedims(velmod[8, :, :, :], (3, 2, 1));
    gam1[inner_id...] .= permutedims(velmod[9, :, :, :], (3, 2, 1));
    gam2[inner_id...] .= permutedims(velmod[10, :, :, :], (3, 2, 1));
    del1[inner_id...] .= permutedims(velmod[11, :, :, :], (3, 2, 1));
    del2[inner_id...] .= permutedims(velmod[12, :, :, :], (3, 2, 1));
    del3[inner_id...] .= permutedims(velmod[13, :, :, :], (3, 2, 1));

    # fill outer domain of nearest values of inner domain
    fill_outer_domain!(domain,vp);
    fill_outer_domain!(domain,vs);
    fill_outer_domain!(domain,rho);
    fill_outer_domain!(domain,eps1);
    fill_outer_domain!(domain,eps2);
    fill_outer_domain!(domain,gam1);
    fill_outer_domain!(domain,gam2);
    fill_outer_domain!(domain,del1);
    fill_outer_domain!(domain,del2);
    fill_outer_domain!(domain,del3);

    # stiffness 
    c33 = vp.^2 .* rho
    c55 = vs.^2 .* rho
    c11 = (2*eps2 .+ 1) .* c33
    c22 = c33 .* (2*eps1 .+ 1)
    c66 = c55 .* (2*gam1 .+ 1)
    c44 = c66 ./ ( 1 .+ gam2)
    c13 = sqrt.( 2*c33 .* (c33-c55) .* del2 +(c33-c55).^2 ) - c55
    c23 = sqrt.( 2*c33 .* (c33-c44) .* del1 +(c33-c44).^2 ) - c44
    c12 = sqrt.( 2*c11 .* (c11-c66) .* del3 +(c11-c66).^2 ) - c66

    # find unique c-tensors and create corresponding index lookup table
    cdata = [(c11[z,y,x], c12[z,y,x], c13[z,y,x],
            c22[z,y,x], c23[z,y,x], c33[z,y,x],
            c44[z,y,x], c55[z,y,x], c66[z,y,x],
            rho[z,y,x])
            for z in 1:dim[1], y in 1:dim[2], x in 1:dim[3]]

    unique_c = unique(vec(cdata))
    c_lud = Dict(t => idx for (idx,t) in enumerate(unique_c))

    c_lookup = Array{Int32}(undef, dim...)  # nz,ny,nx
    for i in 1:dim[1], j in 1:dim[2], k in 1:dim[3]
        c_lookup[i,j,k] = c_lud[cdata[i,j,k]]
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
    angles_phi = deg2rad.(0:2:90)
    VpVs = MVector{3,Float64}(undef)
    UpUs = MMatrix{3,3,Float64}(undef)

    vmax = 0.0
    vmin = Inf
    for c in eachrow(deepcopy(c_tensors))
        c ./= c[end] # density normalized stiffness
        c11, c12, c13, c22, c23, c33, c44, c55, c66, rho_ = c 

        if c55 == 0 || c44 == 0 # skip liquid layer
            continue 
        end

        for theta in angles_theta
            for phi in angles_phi
                n = SVector(sin(phi)*cos(theta), sin(phi)*sin(theta), cos(phi))
                n = n/norm(n)
                gp, gs1, gs2 = group_velocity(VpVs, UpUs, c11, c12, c13, c22, c23, c33, c44, c55, c66, n)
                comb = [norm(gp), norm(gs1), norm(gs2)]
                vmax_ = maximum(comb)
                vmin_ = minimum(comb)
                if vmax_ > vmax 
                    vmax = vmax_
                end
                if vmin_ < vmin 
                    vmin = vmin_
                end
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
