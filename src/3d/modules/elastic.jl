struct Elastic{T, U}
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
    vpmax::U
end


function fill_outer_domain!(domain::Domain,array)

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


function init_elastic(settings::Settings, domain::Domain)

    dim = domain.dim
    FLOAT = settings.float
    inner_id = domain.inner_id

    vp = zeros(FLOAT,dim);
    vs = zeros(FLOAT,dim);
    rho = zeros(FLOAT,dim);
    eps1 = zeros(FLOAT,dim);
    eps2 = zeros(FLOAT,dim);
    gam1 = zeros(FLOAT,dim);
    gam2 = zeros(FLOAT,dim);
    del1 = zeros(FLOAT,dim);
    del2 = zeros(FLOAT,dim);
    del3 = zeros(FLOAT,dim);

    # fill inner domain with velmod 
    vp[inner_id...]   .= permutedims(domain.velmod[4, :, :, :], (3, 2, 1))
    vs[inner_id...]   .= permutedims(domain.velmod[5, :, :, :], (3, 2, 1))
    rho[inner_id...]  .= permutedims(domain.velmod[6, :, :, :], (3, 2, 1))
    eps1[inner_id...] .= permutedims(domain.velmod[7, :, :, :], (3, 2, 1))
    eps2[inner_id...] .= permutedims(domain.velmod[8, :, :, :], (3, 2, 1))
    gam1[inner_id...] .= permutedims(domain.velmod[9, :, :, :], (3, 2, 1))
    gam2[inner_id...] .= permutedims(domain.velmod[10, :, :, :], (3, 2, 1))
    del1[inner_id...] .= permutedims(domain.velmod[11, :, :, :], (3, 2, 1))
    del2[inner_id...] .= permutedims(domain.velmod[12, :, :, :], (3, 2, 1))
    del3[inner_id...] .= permutedims(domain.velmod[13, :, :, :], (3, 2, 1))

    # fill outer domain of nearest values of inner domain
    fill_outer_domain!(domain,vp)
    fill_outer_domain!(domain,vs)
    fill_outer_domain!(domain,rho)
    fill_outer_domain!(domain,eps1)
    fill_outer_domain!(domain,eps2)
    fill_outer_domain!(domain,gam1)
    fill_outer_domain!(domain,gam2)
    fill_outer_domain!(domain,del1)
    fill_outer_domain!(domain,del2)
    fill_outer_domain!(domain,del3)

    if any(vp .== 0) == true  println("vp somewhere zero") end
    vpmax = maximum(vp)

    # C-TENSOR
    c33 = vp.^2 .* rho
    c55 = vs.^2 .* rho
    c11 = (2*eps2 .+ 1) .* c33
    c22 = c33 .* (2*eps1 .+ 1)
    c66 = c55 .* (2*gam1 .+ 1)
    c44 = c66 ./ ( 1 .+ gam2)
    c13 = sqrt.( 2*c33 .* (c33-c55) .* del2 +(c33-c55).^2 ) - c55
    c23 = sqrt.( 2*c33 .* (c33-c44) .* del1 +(c33-c44).^2 ) - c44
    c12 = sqrt.( 2*c11 .* (c11-c66) .* del3 +(c11-c66).^2 ) - c66

    T = typeof(c11)
    U = typeof(vpmax)
    elastic = Elastic{T, U}(c11, c12, c13, c22, c23, c33, c44, c55, c66, rho, vpmax)

    return elastic 
end
