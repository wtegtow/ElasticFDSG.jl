function kernel_args(f::FDSG3D)
    # cpu: points to memory
    # gpu: allocates new memory on gpu 
    return (

    N = f.settings.N,
    nz = f.domain.nz,
    ny = f.domain.ny,
    nx = f.domain.nx,
    dx = f.domain.dx,
    dy = f.domain.dy,
    dz = f.domain.dz,
    dt = f.time.dt,
    dim = f.domain.dim,
    c = f.settings.array(f.settings.c),

    x_evn = f.settings.array(f.pml.x_evn),
    x_odd = f.settings.array(f.pml.x_odd),
    y_evn = f.settings.array(f.pml.y_evn),
    y_odd = f.settings.array(f.pml.y_odd),
    z_evn = f.settings.array(f.pml.z_evn),
    z_odd = f.settings.array(f.pml.z_odd),

    c_tensors  = f.settings.array(f.elastic.c_tensors), 
    c_lookup   = f.settings.array(f.elastic.c_lookup),
    pml_lookup = f.settings.array{Int32}(f.domain.pml_lookup),

    )
end

function field_args(f::FDSG3D)
    return (
        vx = f.settings.array(f.fields.vx),
        vy = f.settings.array(f.fields.vy),
        vz = f.settings.array(f.fields.vz),

        sxx = f.settings.array(f.fields.sxx),
        sxy = f.settings.array(f.fields.sxy),
        sxz = f.settings.array(f.fields.sxz),
        syy = f.settings.array(f.fields.syy),
        syz = f.settings.array(f.fields.syz),
        szz = f.settings.array(f.fields.szz),

        vx_x_old = f.settings.array(f.pml.vx_x_old),
        vy_y_old = f.settings.array(f.pml.vy_y_old),
        vz_z_old = f.settings.array(f.pml.vz_z_old),
        vy_x_old = f.settings.array(f.pml.vy_x_old),
        vx_y_old = f.settings.array(f.pml.vx_y_old),
        vz_x_old = f.settings.array(f.pml.vz_x_old),
        vx_z_old = f.settings.array(f.pml.vx_z_old),
        vz_y_old = f.settings.array(f.pml.vz_y_old),
        vy_z_old = f.settings.array(f.pml.vy_z_old),

        sxx_x_old = f.settings.array(f.pml.sxx_x_old),
        sxy_y_old = f.settings.array(f.pml.sxy_y_old),
        sxz_z_old = f.settings.array(f.pml.sxz_z_old),
        sxy_x_old = f.settings.array(f.pml.sxy_x_old),
        syy_y_old = f.settings.array(f.pml.syy_y_old),
        syz_z_old = f.settings.array(f.pml.syz_z_old),
        sxz_x_old = f.settings.array(f.pml.sxz_x_old),
        syz_y_old = f.settings.array(f.pml.syz_y_old),
        szz_z_old = f.settings.array(f.pml.szz_z_old),
        )
end

# forward kernel: sxx,syy,szz 
@kernel inbounds=true function forward_normal_stresses_kernel!(kernel_params, field_params)

    @unpack N, c, dx, dy, dz, dt, nx, ny, nz, 
            c_tensors, c_lookup, pml_lookup, x_odd, y_evn, z_evn = kernel_params

    @unpack vx, vy, vz, sxx, syy, szz,
            vx_x_old, vy_y_old, vz_z_old = field_params

    z,y,x = @index(Global, NTuple) 
    
    if  2+(N-1) <= z <= nz-N &&
        2+(N-1) <= y <= ny-N &&
        1+(N-1) <= x <= nx-N

        # EFFECTIVE MEDIA PARAMETER
        cidx = c_lookup[z,y,x]
        cidx_xp = c_lookup[z,y,x+1]  
        c11eff = 0.5f0 * (c_tensors[cidx, 1] + c_tensors[cidx_xp, 1])
        c12eff = 0.5f0 * (c_tensors[cidx, 2] + c_tensors[cidx_xp, 2])
        c13eff = 0.5f0 * (c_tensors[cidx, 3] + c_tensors[cidx_xp, 3])
        c22eff = 0.5f0 * (c_tensors[cidx, 4] + c_tensors[cidx_xp, 4])
        c23eff = 0.5f0 * (c_tensors[cidx, 5] + c_tensors[cidx_xp, 5])
        c33eff = 0.5f0 * (c_tensors[cidx, 6] + c_tensors[cidx_xp, 6])

        # SPACE DERIVATIVE
        vx_x = 0
        vy_y = 0
        vz_z = 0
        for i in 1:N
            vx_x += c[i]/dx * (vx[z,y,x+i]-vx[z,y,x-(i-1)])
            vy_y += c[i]/dy * (vy[z,y+(i-1),x]-vy[z,y-i,x])
            vz_z += c[i]/dz * (vz[z+(i-1),y,x]-vz[z-i,y,x])
        end

        # UPDATE PML REGION
        if pml_lookup[z,y,x] != -1 
            pml_idx = pml_lookup[z,y,x]

            # UPDATE SPACE DERIVATIVE
            vx_x_old[pml_idx] = x_odd[2,x] * vx_x_old[pml_idx] + x_odd[1,x] * vx_x
            vy_y_old[pml_idx] = y_evn[2,y] * vy_y_old[pml_idx] + y_evn[1,y] * vy_y
            vz_z_old[pml_idx] = z_evn[2,z] * vz_z_old[pml_idx] + z_evn[1,z] * vz_z
            
            vx_x = vx_x/x_odd[3,x] + vx_x_old[pml_idx]
            vy_y = vy_y/y_evn[3,y] + vy_y_old[pml_idx]
            vz_z = vz_z/z_evn[3,z] + vz_z_old[pml_idx]

        end

        # UPDATE STRESSES
        sxx[z,y,x] += dt * (c11eff * vx_x +  c12eff * vy_y + c13eff * vz_z) 
        syy[z,y,x] += dt * (c12eff * vx_x +  c22eff * vy_y + c23eff * vz_z)
        szz[z,y,x] += dt * (c13eff * vx_x +  c23eff * vy_y + c33eff * vz_z)

    end

end;
function forward_normal_stresses!(kernel_params, field_params; backend=CPU(), block_size=(7,7,7))
    kernel = forward_normal_stresses_kernel!(backend, block_size)
    kernel(kernel_params, field_params, ndrange = kernel_params.dim)
end;

# forward kernel: sxy 
@kernel inbounds=true function forward_sxy_kernel!(kernel_params, field_params)

    @unpack N, c, dx, dy, dz, dt, nx, ny, nz, 
            c_tensors, c_lookup, pml_lookup, x_evn, y_odd = kernel_params

    @unpack vx, vy, sxy, vy_x_old, vx_y_old = field_params

    z,y,x = @index(Global, NTuple) 

    if         1  <= z <= nz      &&
          1+(N-1) <= y <= ny-N    &&
          2+(N-1) <= x <= nx-(N-1)

        # EFFECTIVE MEDIA PARAMETER
        cidx = c_lookup[z,y,x]
        cidx_yp = c_lookup[z,y+1,x]  
        c66eff = 0.5f0 * (c_tensors[cidx, 9] + c_tensors[cidx_yp, 9])
        
        vy_x = 0
        vx_y = 0
        for i in 1:N
            vy_x +=  c[i]/dx * (vy[z,y,x+(i-1)]-vy[z,y,x-i])
            vx_y +=  c[i]/dy * (vx[z,y+i,x]-vx[z,y-(i-1),x])
        end
        
        # UPDATE PML VALUES
        if pml_lookup[z,y,x] != -1 
            pml_idx = pml_lookup[z,y,x]

            # UPDATE SPACE DERIVATIVE
            vy_x_old[pml_idx] = x_evn[2,x] * vy_x_old[pml_idx] + x_evn[1,x] * vy_x
            vx_y_old[pml_idx] = y_odd[2,y] * vx_y_old[pml_idx] + y_odd[1,y] * vx_y
            
            vy_x = vy_x/x_evn[3,x] + vy_x_old[pml_idx]
            vx_y = vx_y/y_odd[3,y] + vx_y_old[pml_idx]

        end
        
        # UPDATE STRESSES
        sxy[z,y,x] += dt * c66eff * (vy_x + vx_y)
       
    end
end;
function forward_sxy!(kernel_params, field_params; backend=CPU(), block_size=(7,7,7))
    kernel = forward_sxy_kernel!(backend, block_size)
    kernel(kernel_params, field_params, ndrange = kernel_params.dim)
end;

# forward kernel: sxz 
@kernel inbounds=true function forward_sxz_kernel!(kernel_params, field_params)

    @unpack N, c, dx, dy, dz, dt, nx, ny, nz, 
            c_tensors, c_lookup, pml_lookup, x_evn, z_odd = kernel_params

    @unpack vx, vz, sxz, vz_x_old, vx_z_old = field_params

    z,y,x = @index(Global, NTuple) 

    if  1+(N-1) <= z <= nz-N     &&
        1       <= y <= ny       && 
        2+(N-1) <= x <= nx-(N-1)
    
        # EFFECTIVE MEDIA PARAMETER 
        cidx = c_lookup[z,y,x]
        cidx_zp = c_lookup[z+1,y,x]  
        c55eff = 0.5f0 * (c_tensors[cidx, 8] + c_tensors[cidx_zp, 8])
      
        # SPACE DERIVATIVE
        vz_x = 0
        vx_z = 0
        for i in 1:N
            vz_x += c[i]/dx * (vz[z,y,x+(i-1)]-vz[z,y,x-i])
            vx_z += c[i]/dz * (vx[z+i,y,x]-vx[z-(i-1),y,x])
        end
        
        # UPDATE PML VALUES
        if pml_lookup[z,y,x] != -1 
            pml_idx = pml_lookup[z,y,x]

            # UPDATE SPACE DERIVATIVE
            vz_x_old[pml_idx] = x_evn[2,x] * vz_x_old[pml_idx] + x_evn[1,x] * vz_x
            vx_z_old[pml_idx] = z_odd[2,z] * vx_z_old[pml_idx] + z_odd[1,z] * vx_z
        
            vz_x = vz_x/x_evn[3,x] + vz_x_old[pml_idx]
            vx_z = vx_z/z_odd[3,z] + vx_z_old[pml_idx] 

        end
        
        # UPDATE STRESSES
        sxz[z,y,x] += dt * c55eff * (vz_x + vx_z)
    
    end
end;
function forward_sxz!(kernel_params, field_params; backend=CPU(), block_size=(7,7,7))
    kernel = forward_sxz_kernel!(backend, block_size)
    kernel(kernel_params, field_params, ndrange = kernel_params.dim)
end;

# forward kernel: syz
@kernel inbounds=true function forward_syz_kernel!(kernel_params, field_params)

    @unpack N, c, dx, dy, dz, dt, nx, ny, nz, 
            c_tensors, c_lookup, pml_lookup, y_odd, z_odd = kernel_params

    @unpack vy, vz, syz, vz_y_old, vy_z_old = field_params

    z,y,x = @index(Global, NTuple) 

    if  1+(N-1) <= z <= nz-N  &&
        1+(N-1) <= y <= ny-N  && 
        1       <= x <= nx
    
        # EFFECTIVE MEDIA PARAMETER
        cidx      = c_lookup[z,y,x]
        cidx_zp   = c_lookup[z+1,y,x]  
        cidx_yp   = c_lookup[z,y+1,x]
        cidx_zpyp = c_lookup[z+1,y+1,x]
        c44eff = 0.25f0 * (c_tensors[cidx, 7]    + c_tensors[cidx_zp, 7] +
                           c_tensors[cidx_yp, 7] + c_tensors[cidx_zpyp, 7])
        
        # SPACE DERIVATIVE
        vz_y = 0
        vy_z = 0
        for i in 1:N 
            vz_y += c[i]/dy * (vz[z,y+i,x]-vz[z,y-(i-1),x])
            vy_z += c[i]/dz * (vy[z+i,y,x]-vy[z-(i-1),y,x])
        end
        
        # UPDATE PML VALUES
        if pml_lookup[z,y,x] != -1 
            pml_idx = pml_lookup[z,y,x] 

            # UPDATE SPACE DERIVATIVE
            vz_y_old[pml_idx] = y_odd[2,y] * vz_y_old[pml_idx] + y_odd[1,y] * vz_y
            vy_z_old[pml_idx] = z_odd[2,z] * vy_z_old[pml_idx] + z_odd[1,z] * vy_z
            
            vz_y = vz_y/y_odd[3,y] + vz_y_old[pml_idx]
            vy_z = vy_z/z_odd[3,z] + vy_z_old[pml_idx]

        end
        
        # UPDATE STRESS
        syz[z,y,x] += dt * c44eff * (vz_y + vy_z)

    end
end;
function forward_syz!(kernel_params, field_params; backend=CPU(), block_size=(7,7,7))
    kernel = forward_syz_kernel!(backend, block_size)
    kernel(kernel_params, field_params, ndrange = kernel_params.dim)
end;

# forward kernel: vx
@kernel inbounds=true function forward_vx_kernel!(kernel_params, field_params)

    @unpack N, c, dx, dy, dz, dt, nx, ny, nz, 
            c_tensors, c_lookup, pml_lookup, x_evn, y_evn, z_evn = kernel_params

    @unpack vx, sxx, sxy, sxz, 
            sxx_x_old, sxy_y_old, sxz_z_old = field_params


    z,y,x = @index(Global, NTuple) 
    
    if  2+(N-1) <= z <= nz-(N-1) && 
        2+(N-1) <= y <= ny-(N-1) && 
        2+(N-1) <= x <= nx-(N-1)
    
        # EFFECTIVE MEDIA PARAMETER
        cidx = c_lookup[z,y,x]
        rho  = c_tensors[cidx, 10] 

        # SPACE DERIVATIVE
        sxx_x = 0
        sxy_y = 0
        sxz_z = 0
        for i in 1:N 
            sxx_x += c[i]/dx * (sxx[z,y,x+(i-1)]-sxx[z,y,x-i])
            sxy_y += c[i]/dy * (sxy[z,y+(i-1),x]-sxy[z,y-i,x])
            sxz_z += c[i]/dz * (sxz[z+(i-1),y,x]-sxz[z-i,y,x])
        end
        
        # UPDATE PML VALUES 
        if pml_lookup[z,y,x] != -1 
            pml_idx = pml_lookup[z,y,x] 

            # UPDATE SPACE DERIVATIVE
            sxx_x_old[pml_idx] = x_evn[2,x] * sxx_x_old[pml_idx] + x_evn[1,x] * sxx_x
            sxy_y_old[pml_idx] = y_evn[2,y] * sxy_y_old[pml_idx] + y_evn[1,y] * sxy_y
            sxz_z_old[pml_idx] = z_evn[2,z] * sxz_z_old[pml_idx] + z_evn[1,z] * sxz_z
            
            sxx_x = sxx_x/x_evn[3,x] + sxx_x_old[pml_idx]
            sxy_y = sxy_y/y_evn[3,y] + sxy_y_old[pml_idx]
            sxz_z = sxz_z/z_evn[3,z] + sxz_z_old[pml_idx]
        end
        
        # UPDATE VELOCITY
        vx[z,y,x] += dt/rho * (sxx_x + sxy_y + sxz_z)
    
    end
end;
function forward_vx!(kernel_params, field_params; backend=CPU(), block_size=(7,7,7))
    kernel = forward_vx_kernel!(backend, block_size)
    kernel(kernel_params, field_params, ndrange = kernel_params.dim)
end;

# forward kernel: vy
@kernel inbounds=true function forward_vy_kernel!(kernel_params, field_params)

    @unpack N, c, dx, dy, dz, dt, nx, ny, nz, 
            c_tensors, c_lookup, pml_lookup, x_odd, y_odd, z_evn  = kernel_params

    @unpack vy, sxy, syy, syz, 
            sxy_x_old, syy_y_old, syz_z_old = field_params
    
    z,y,x = @index(Global, NTuple) 
    
    if  2+(N-1) <= z <= nz-(N-1)  &&
        1+(N-1) <= y <= ny-N      &&
        1+(N-1) <= x <= nx-N
    
        # EFFECTIVE MEDIA PARAMETER
        cidx      = c_lookup[z,y,x]
        cidx_yp   = c_lookup[z,y+1,x]  
        cidx_xp   = c_lookup[z,y,x+1]
        cidx_ypxp = c_lookup[z,y+1,x+1]
        rhoeff = 0.25f0 * (c_tensors[cidx, 10]    + c_tensors[cidx_yp, 10] +
                           c_tensors[cidx_xp, 10] + c_tensors[cidx_ypxp, 10])
        
        # SPACE DERIVATIVE
        sxy_x = 0
        syy_y = 0
        syz_z  = 0
        for i in 1:N 
            sxy_x += c[i]/dx * (sxy[z,y,x+i] - sxy[z,y,x-(i-1)])
            syy_y += c[i]/dy * (syy[z,y+i,x]-syy[z,y-(i-1),x])
            syz_z += c[i]/dz * (syz[z+(i-1),y,x]-syz[z-i,y,x])
        end
        
        # UPDATE PML VALUES 
        if pml_lookup[z,y,x] != -1 
            pml_idx = pml_lookup[z,y,x] 

            # UPDATE SPACE DERIVATIVE
            sxy_x_old[pml_idx] = x_odd[2,x] * sxy_x_old[pml_idx] + x_odd[1,x] * sxy_x
            syy_y_old[pml_idx] = y_odd[2,y] * syy_y_old[pml_idx] + y_odd[1,y] * syy_y
            syz_z_old[pml_idx] = z_evn[2,z] * syz_z_old[pml_idx] + z_evn[1,z] * syz_z
        
            sxy_x = sxy_x/x_odd[3,x] + sxy_x_old[pml_idx]
            syy_y = syy_y/y_odd[3,y] + syy_y_old[pml_idx]
            syz_z = syz_z/z_evn[3,z] + syz_z_old[pml_idx]
        end
        
        # UPDATE VELOCITY
        vy[z,y,x] += dt/rhoeff * (sxy_x + syy_y + syz_z)
    
    end
end;
function forward_vy!(kernel_params, field_params; backend=CPU(), block_size=(7,7,7))
    kernel = forward_vy_kernel!(backend, block_size)
    kernel(kernel_params, field_params, ndrange = kernel_params.dim)
end;

# forward kernel: vz
@kernel inbounds=true function forward_vz_kernel!(kernel_params, field_params)

    @unpack N, c, dx, dy, dz, dt, nx, ny, nz, 
            c_tensors, c_lookup, pml_lookup, x_odd, y_evn, z_odd  = kernel_params

    @unpack vz, sxz, syz, szz, 
            sxz_x_old, syz_y_old, szz_z_old = field_params

    z,y,x = @index(Global, NTuple) 
    
    if  1+(N-1) <= z <= nz-N      && 
        2+(N-1) <= y <= ny-(N-1)  &&
        1+(N-1) <= x <= nx-N
    
        # EFFECTIVE MEDIA PARAMETER
        cidx      = c_lookup[z,y,x]
        cidx_zp   = c_lookup[z+1,y,x]  
        cidx_xp   = c_lookup[z,y,x+1]
        cidx_zpxp = c_lookup[z+1,y,x+1]
        rhoeff = 0.25f0 * (c_tensors[cidx, 10]    + c_tensors[cidx_zp, 10] +
                           c_tensors[cidx_xp, 10] + c_tensors[cidx_zpxp, 10])
        
        
        # SPACE DERIVATIVE
        sxz_x = 0
        syz_y = 0
        szz_z = 0
        for i in 1:N 
            sxz_x += c[i]/dx * (sxz[z,y,x+i]-sxz[z,y,x-(i-1)])
            syz_y += c[i]/dy * (syz[z,y+(i-1),x]-syz[z,y-i,x])
            szz_z += c[i]/dz * (szz[z+i,y,x]-szz[z-(i-1),y,x]) 
        end
        
        # UPDATE PML VALUES 
        if pml_lookup[z,y,x] != -1 
            pml_idx = pml_lookup[z,y,x] 

            # UPDATE SPACE DERIVATIVE
            sxz_x_old[pml_idx] = x_odd[2,x] * sxz_x_old[pml_idx] + x_odd[1,x] * sxz_x
            syz_y_old[pml_idx] = y_evn[2,y] * syz_y_old[pml_idx] + y_evn[1,y] * syz_y
            szz_z_old[pml_idx] = z_odd[2,z] * szz_z_old[pml_idx] + z_odd[1,z] * szz_z
            
            sxz_x = sxz_x/x_odd[3,x] + sxz_x_old[pml_idx]
            syz_y = syz_y/y_evn[3,y] + syz_y_old[pml_idx]
            szz_z = szz_z/z_odd[3,z] + szz_z_old[pml_idx]
        end
        
        # UPDATE VELOCITY
        vz[z,y,x] += dt/rhoeff * (sxz_x + syz_y + szz_z)
    
    end
end
function forward_vz!(kernel_params, field_params; backend=CPU(), block_size=(7,7,7))
    kernel = forward_vz_kernel!(backend, block_size)
    kernel(kernel_params, field_params, ndrange = kernel_params.dim)
end;