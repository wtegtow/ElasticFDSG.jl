struct Elastic{T, U}
    c11::T 
    c13::T 
    c33::T 
    c55::T 
    rho::T 
    vpmax::U

end


function fill_outer_domain!(domain::Domain, array)

    ny = domain.ny 
    nx = domain.nx 
    inner_id = domain.inner_id

    y_min, y_max = inner_id[1][begin], inner_id[1][end]
    x_min, x_max = inner_id[2][begin], inner_id[2][end]

    for y in 1:ny 
        for x in 1:nx 
            if  y_min <= y <= y_max && 
                x_min <= x <= x_max
                continue  
            else
            nearest_y = clamp(y, y_min, y_max)
            nearest_x = clamp(x, x_min, x_max)
            array[y,x] = array[nearest_y, nearest_x]
            end 
        end 
    end
end;


function init_elastic(settings::Settings, domain::Domain)
    
    vp = zeros(settings.float,domain.dim)
    vs = similar(vp);
    rho = similar(vp);
    eps0 = similar(vp);
    del0 = similar(vp);

    # fill inner domain of elastic 
    vp[domain.inner_id...]   .= domain.velmod[3,:,:]
    vs[domain.inner_id...]   .= domain.velmod[4,:,:]
    rho[domain.inner_id...]  .= domain.velmod[5,:,:]
    eps0[domain.inner_id...] .= domain.velmod[6,:,:]
    del0[domain.inner_id...] .= domain.velmod[7,:,:]

    # fill outer domain with nearest inner values
    fill_outer_domain!(domain,vp)
    fill_outer_domain!(domain,vs)
    fill_outer_domain!(domain,rho)
    fill_outer_domain!(domain,eps0)
    fill_outer_domain!(domain,del0)

    if any(vp .== 0) == true  
        error("Error in elastic.jl. vp somewhere zero") 
    end
    
    vpmax = maximum(vp)

    lam = similar(vp);
    mu = similar(vp);
    lam[:,:] .= rho .* (vp.^2 - 2 .* vs.^2);
    mu[:,:]  .= rho .* vs.^2;
  
    c11 = (lam + 2 .* mu) .* (2*eps0 .+ 1)
    c13 = lam + del0 .* (lam + 2*mu)
    c33 = lam + 2*mu
    c55 = mu 

    # types 
    T = typeof(c11)
    U = typeof(vpmax)

    elastic = Elastic{T, U}(c11, c13, c33, c55, rho, vpmax)
    return elastic

end; 