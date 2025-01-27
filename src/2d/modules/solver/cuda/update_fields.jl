# =============================================================================
# update field cuda kernel
# =============================================================================

function update_normal_stresses_cuda(
    N, c, dt,
    ny, nx, dy, dx, 
    c11,c13,c33,
    vx,vy,
    sxx,syy,
    a_x_odd, b_x_odd, K_x_odd, 
    a_y_evn, b_y_evn, K_y_evn, 
    vx_x_old, vy_y_old,
    pml_hashmap)

    y = threadIdx().x + (blockIdx().x - 1) * blockDim().x
    x = threadIdx().y + (blockIdx().y - 1) * blockDim().y

    if  2+(N-1) <= y <= ny-(N-1) &&
        1+(N-1) <= x <= nx-N

        @inbounds begin

        # EFFECTIVE ELASTIC PARAMETER
        c11eff = 0.5f0 * (c11[y,x+1] + c11[y,x])
        c13eff = 0.5f0 * (c13[y,x+1] + c13[y,x])
        c33eff = 0.5f0 * (c33[y,x+1] + c33[y,x])

        # SPACE DERIVATIVE
        vx_x = 0
        vy_y = 0
        for i in 1:N
            vx_x += c[i]/dx * (vx[y,x+i]-vx[y,x-(i-1)])
            vy_y += c[i]/dy * (vy[y+(i-1),x]-vy[y-i,x])
        end

        # UPDATE PML REGIONS 
        if pml_hashmap[y,x] != -1 

            pml_idx = pml_hashmap[y,x]

            vx_x_old[pml_idx] = b_x_odd[x] * vx_x_old[pml_idx] + a_x_odd[x] * vx_x
            vy_y_old[pml_idx] = b_y_evn[y] * vy_y_old[pml_idx] + a_y_evn[y] * vy_y

            # UPDATE SPACE DERIVATIVE
            vx_x = vx_x/K_x_odd[x] + vx_x_old[pml_idx]
            vy_y = vy_y/K_y_evn[y] + vy_y_old[pml_idx]

        end

        # UPDATE STRESSES
        sxx[y,x] += dt * (c11eff * vx_x +  c13eff * vy_y)
        syy[y,x] += dt * (c13eff * vx_x +  c33eff * vy_y)

        end
    end

    return nothing
end;

function update_sxy_cuda(
    N, c, dt,
    ny, nx, dy, dx, 
    c55,
    vx,vy,
    sxy,
    a_x_evn, b_x_evn, K_x_evn, 
    a_y_odd, b_y_odd, K_y_odd, 
    vy_x_old,vx_y_old,
    pml_hashmap)

    y = threadIdx().x + (blockIdx().x - 1) * blockDim().x
    x = threadIdx().y + (blockIdx().y - 1) * blockDim().y

    if  1+(N-1) <= y <= ny-N &&
        2+(N-1) <= x <= nx-(N-1)

        @inbounds begin

        # EFFECTIVE ELASTIC PARAMETER
        c55eff = 0.5f0 * (c55[y+1,x] + c55[y,x])

        # SPACE DERIVATIVE
        vy_x = 0
        vx_y = 0
        for i in 1:N
            vy_x +=  c[i]/dx * (vy[y,x+(i-1)]-vy[y,x-i])
            vx_y +=  c[i]/dy * (vx[y+i,x]-vx[y-(i-1),x])
        end

        # UPDATE PML REGION
        if pml_hashmap[y,x] != -1 

            pml_idx = pml_hashmap[y,x]

            vy_x_old[pml_idx] = b_x_evn[x] * vy_x_old[pml_idx] + a_x_evn[x] * vy_x
            vx_y_old[pml_idx] = b_y_odd[y] * vx_y_old[pml_idx] + a_y_odd[y] * vx_y

            # UPDATE SPACE DERIVATIVE
            vy_x = vy_x/K_x_evn[x] + vy_x_old[pml_idx]
            vx_y = vx_y/K_y_odd[y] + vx_y_old[pml_idx]

        end

        # UPDATE STRESSES
        sxy[y,x] += dt * c55eff * (vy_x + vx_y)

        end
    end

    return nothing
end




function update_vx_cuda(N, c, dt,
    ny, nx, dy, dx, 
    rho,
    sxx,sxy,
    vx,
    a_x_evn, b_x_evn, K_x_evn, 
    a_y_evn, b_y_evn, K_y_evn, 
    sxx_x_old,sxy_y_old,
    pml_hashmap)

    y = threadIdx().x + (blockIdx().x - 1) * blockDim().x
    x = threadIdx().y + (blockIdx().y - 1) * blockDim().y

    if  2+(N-1) <= y <= ny-(N-1) &&
        2+(N-1) <= x <= nx-(N-1)

        @inbounds begin

        # EFFECTIVE ELASTIC PARAMETER
        # rho is defined on vx-grid points

        # SPACE DERIVATIVE
        sxx_x = 0
        sxy_y = 0
        for i in 1:N 
            sxx_x += c[i]/dx * (sxx[y,x+(i-1)]-sxx[y,x-i])
            sxy_y += c[i]/dy * (sxy[y+(i-1),x]-sxy[y-i,x])
        end

        # UPDATE PML REGION
        if pml_hashmap[y,x] != -1 

            pml_idx = pml_hashmap[y,x]

            sxx_x_old[pml_idx] = b_x_evn[x] * sxx_x_old[pml_idx] + a_x_evn[x] * sxx_x
            sxy_y_old[pml_idx] = b_y_evn[y] * sxy_y_old[pml_idx] + a_y_evn[y] * sxy_y

            # UPDATE SPACE DERIVATIVE
            sxx_x = sxx_x/K_x_evn[x] + sxx_x_old[pml_idx]
            sxy_y = sxy_y/K_y_evn[y] + sxy_y_old[pml_idx]
        end

        # UPDATE VELOCITY
        vx[y,x] += dt/rho[y,x] * (sxx_x + sxy_y)

        end
    end

    return nothing
end



function update_vy_cuda(N, c, dt,
    ny, nx, dy, dx, 
    rho,
    sxy,syy,
    vy,
    a_x_odd, b_x_odd, K_x_odd, 
    a_y_odd, b_y_odd, K_y_odd, 
    sxy_x_old,syy_y_old,
    pml_hashmap)

    y = threadIdx().x + (blockIdx().x - 1) * blockDim().x
    x = threadIdx().y + (blockIdx().y - 1) * blockDim().y

    if  1+(N-1) <= y <= ny-N &&
        1+(N-1) <= x <= nx-N

        @inbounds begin

        # EFFECTIVE ELASTIC PARAMETER
        rhoeff = 0.25f0 * (rho[y+1,x] + rho[y,x] + rho[y,x+1] + rho[y+1,x+1])
        
        # SPACE DERIVATIVE
        sxy_x = 0
        syy_y = 0
        for i in 1:N 
            sxy_x += c[i]/dx * (sxy[y,x+i] - sxy[y,x-(i-1)])
            syy_y += c[i]/dy * (syy[y+i,x] - syy[y-(i-1),x])
        end

        # UPDATE PML REGION
        if pml_hashmap[y,x] != -1 

            pml_idx = pml_hashmap[y,x]

            sxy_x_old[pml_idx] = b_x_odd[x] * sxy_x_old[pml_idx] + a_x_odd[x] * sxy_x
            syy_y_old[pml_idx] = b_y_odd[y] * syy_y_old[pml_idx] + a_y_odd[y] * syy_y

            # UPDATE SPACE DERIVATIVE
            sxy_x = sxy_x/K_x_odd[x] + sxy_x_old[pml_idx]
            syy_y = syy_y/K_y_odd[y] + syy_y_old[pml_idx]

        end
        
        # UPDATE VELOCITY
        vy[y,x] += dt/rhoeff * (sxy_x + syy_y)
        
        end
    end

    return nothing
end;

