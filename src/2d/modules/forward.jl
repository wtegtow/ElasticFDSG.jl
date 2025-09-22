function kernel_args(f::FDSG2D)
    # cpu: points to memory
    # gpu: allocates new memory on gpu 
    return (

    N = f.settings.N,
    nz = f.domain.nz,
    nx = f.domain.nx,
    dx = f.domain.dx,
    dz = f.domain.dz,
    dt = f.time.dt,
    dim = f.domain.dim,
    c = f.settings.array(f.settings.c),

    x_evn = f.settings.array(f.pml.x_evn),
    x_odd = f.settings.array(f.pml.x_odd),
    z_evn = f.settings.array(f.pml.z_evn),
    z_odd = f.settings.array(f.pml.z_odd),

    c_tensors  = f.settings.array(f.elastic.c_tensors), 
    c_lookup   = f.settings.array(f.elastic.c_lookup),
    pml_lookup = f.settings.array{Int32}(f.domain.pml_lookup),

    )
end

function field_args(f::FDSG2D)
    return (
        vx = f.settings.array(f.fields.vx),
        vz = f.settings.array(f.fields.vz),

        sxx = f.settings.array(f.fields.sxx),
        sxz = f.settings.array(f.fields.sxz),
        szz = f.settings.array(f.fields.szz),

        vx_x_old = f.settings.array(f.pml.vx_x_old),
        vz_z_old = f.settings.array(f.pml.vz_z_old),
        vz_x_old = f.settings.array(f.pml.vz_x_old),
        vx_z_old = f.settings.array(f.pml.vx_z_old),
      
        sxx_x_old = f.settings.array(f.pml.sxx_x_old),
        sxz_z_old = f.settings.array(f.pml.sxz_z_old),
        sxz_x_old = f.settings.array(f.pml.sxz_x_old),
        szz_z_old = f.settings.array(f.pml.szz_z_old),

        )
end

# forward kernel: sxx,syy,szz  y
@kernel inbounds=true function forward_normal_stresses_kernel!(kernel_params, field_params)

    @unpack N, c, dx, dz, dt, nx, nz, 
            c_tensors, c_lookup, pml_lookup, x_odd, z_evn = kernel_params

    @unpack vx, vz, sxx, szz,
            vx_x_old, vz_z_old = field_params

    z,x = @index(Global, NTuple) 
    
    if  2+(N-1) <= z <= nz-N &&
        1+(N-1) <= x <= nx-N

        # EFFECTIVE MEDIA PARAMETER
        cidx = c_lookup[z,x]
        cidx_xp = c_lookup[z,x+1]  
        c11eff = 0.5f0 * (c_tensors[cidx, 1] + c_tensors[cidx_xp, 1])
        c13eff = 0.5f0 * (c_tensors[cidx, 2] + c_tensors[cidx_xp, 2])
        c33eff = 0.5f0 * (c_tensors[cidx, 3] + c_tensors[cidx_xp, 3])

        # SPACE DERIVATIVE
        vx_x = 0
        vz_z = 0
        for i in 1:N
            vx_x += c[i]/dx * (vx[z,x+i]-vx[z,x-(i-1)])
            vz_z += c[i]/dz * (vz[z+(i-1),x]-vz[z-i,x])
        end

        # UPDATE PML REGION
        if pml_lookup[z,x] != -1 
            pml_idx = pml_lookup[z,x]

            # UPDATE SPACE DERIVATIVE
            vx_x_old[pml_idx] = x_odd[2,x] * vx_x_old[pml_idx] + x_odd[1,x] * vx_x
            vz_z_old[pml_idx] = z_evn[2,z] * vz_z_old[pml_idx] + z_evn[1,z] * vz_z
            
            vx_x = vx_x/x_odd[3,x] + vx_x_old[pml_idx]
            vz_z = vz_z/z_evn[3,z] + vz_z_old[pml_idx]

        end

        # UPDATE STRESSES
        sxx[z,x] += dt * (c11eff * vx_x + c13eff * vz_z) 
        szz[z,x] += dt * (c13eff * vx_x + c33eff * vz_z)

    end

end;
function forward_normal_stresses!(kernel_params, field_params; backend=CPU(), block_size=(32,32))
    kernel = forward_normal_stresses_kernel!(backend, block_size)
    kernel(kernel_params, field_params, ndrange = kernel_params.dim)
end;

# forward kernel: sxz 
@kernel inbounds=true function forward_sxz_kernel!(kernel_params, field_params)

    @unpack N, c, dx, dz, dt, nx, nz, 
            c_tensors, c_lookup, pml_lookup, x_evn, z_odd = kernel_params

    @unpack vx, vz, sxz, vz_x_old, vx_z_old = field_params

    z,x = @index(Global, NTuple) 

    if  1+(N-1) <= z <= nz-N     &&
        2+(N-1) <= x <= nx-(N-1)
    
        # EFFECTIVE MEDIA PARAMETER 
        cidx = c_lookup[z,x]
        cidx_zp = c_lookup[z+1,x]  
        c55eff = 0.5f0 * (c_tensors[cidx, 4] + c_tensors[cidx_zp, 4])
      
        # SPACE DERIVATIVE
        vz_x = 0
        vx_z = 0
        for i in 1:N
            vz_x += c[i]/dx * (vz[z,x+(i-1)]-vz[z,x-i])
            vx_z += c[i]/dz * (vx[z+i,x]-vx[z-(i-1),x])
        end
        
        # UPDATE PML VALUES
        if pml_lookup[z,x] != -1 
            pml_idx = pml_lookup[z,x]

            # UPDATE SPACE DERIVATIVE
            vz_x_old[pml_idx] = x_evn[2,x] * vz_x_old[pml_idx] + x_evn[1,x] * vz_x
            vx_z_old[pml_idx] = z_odd[2,z] * vx_z_old[pml_idx] + z_odd[1,z] * vx_z
        
            vz_x = vz_x/x_evn[3,x] + vz_x_old[pml_idx]
            vx_z = vx_z/z_odd[3,z] + vx_z_old[pml_idx] 

        end
        
        # UPDATE STRESSES
        sxz[z,x] += dt * c55eff * (vz_x + vx_z)
    
    end
end;
function forward_sxz!(kernel_params, field_params; backend=CPU(), block_size=(32,32))
    kernel = forward_sxz_kernel!(backend, block_size)
    kernel(kernel_params, field_params, ndrange = kernel_params.dim)
end;

# forward kernel: vx
@kernel inbounds=true function forward_vx_kernel!(kernel_params, field_params)

    @unpack N, c, dx, dz, dt, nx, nz, 
            c_tensors, c_lookup, pml_lookup, x_evn, z_evn = kernel_params

    @unpack vx, sxx, sxz, 
            sxx_x_old, sxz_z_old = field_params


    z,x = @index(Global, NTuple) 
    
    if  2+(N-1) <= z <= nz-(N-1) && 
        2+(N-1) <= x <= nx-(N-1)
    
        # EFFECTIVE MEDIA PARAMETER
        cidx = c_lookup[z,x]
        rho  = c_tensors[cidx, 5] 

        # SPACE DERIVATIVE
        sxx_x = 0
        sxz_z = 0
        for i in 1:N 
            sxx_x += c[i]/dx * (sxx[z,x+(i-1)]-sxx[z,x-i])
            sxz_z += c[i]/dz * (sxz[z+(i-1),x]-sxz[z-i,x])
        end
        
        # UPDATE PML VALUES 
        if pml_lookup[z,x] != -1 
            pml_idx = pml_lookup[z,x] 

            # UPDATE SPACE DERIVATIVE
            sxx_x_old[pml_idx] = x_evn[2,x] * sxx_x_old[pml_idx] + x_evn[1,x] * sxx_x
            sxz_z_old[pml_idx] = z_evn[2,z] * sxz_z_old[pml_idx] + z_evn[1,z] * sxz_z
            
            sxx_x = sxx_x/x_evn[3,x] + sxx_x_old[pml_idx]
            sxz_z = sxz_z/z_evn[3,z] + sxz_z_old[pml_idx]
        end
        
        # UPDATE VELOCITY
        vx[z,x] += dt/rho * (sxx_x + sxz_z)
    
    end
end;
function forward_vx!(kernel_params, field_params; backend=CPU(), block_size=(32,32))
    kernel = forward_vx_kernel!(backend, block_size)
    kernel(kernel_params, field_params, ndrange = kernel_params.dim)
end;

# forward kernel: vz
@kernel inbounds=true function forward_vz_kernel!(kernel_params, field_params)

    @unpack N, c, dx, dz, dt, nx, nz, 
            c_tensors, c_lookup, pml_lookup, x_odd, z_odd  = kernel_params

    @unpack vz, sxz, szz, 
            sxz_x_old, szz_z_old = field_params

    z,x = @index(Global, NTuple) 
    
    if  1+(N-1) <= z <= nz-N      && 
        1+(N-1) <= x <= nx-N
    
        # EFFECTIVE MEDIA PARAMETER
        cidx      = c_lookup[z,x]
        cidx_zp   = c_lookup[z+1,x]  
        cidx_xp   = c_lookup[z,x+1]
        cidx_zpxp = c_lookup[z+1,x+1]
        rhoeff = 0.25f0 * (c_tensors[cidx, 5]    + c_tensors[cidx_zp, 5] +
                           c_tensors[cidx_xp, 5] + c_tensors[cidx_zpxp, 5])
        
        # SPACE DERIVATIVE
        sxz_x = 0
        szz_z = 0
        for i in 1:N 
            sxz_x += c[i]/dx * (sxz[z,x+i]-sxz[z,x-(i-1)])
            szz_z += c[i]/dz * (szz[z+i,x]-szz[z-(i-1),x]) 
        end
        
        # UPDATE PML VALUES 
        if pml_lookup[z,x] != -1 
            pml_idx = pml_lookup[z,x] 

            # UPDATE SPACE DERIVATIVE
            sxz_x_old[pml_idx] = x_odd[2,x] * sxz_x_old[pml_idx] + x_odd[1,x] * sxz_x
            szz_z_old[pml_idx] = z_odd[2,z] * szz_z_old[pml_idx] + z_odd[1,z] * szz_z
            
            sxz_x = sxz_x/x_odd[3,x] + sxz_x_old[pml_idx]
            szz_z = szz_z/z_odd[3,z] + szz_z_old[pml_idx]
        end
        
        # UPDATE VELOCITY
        vz[z,x] += dt/rhoeff * (sxz_x + szz_z)
    
    end
end
function forward_vz!(kernel_params, field_params; backend=CPU(), block_size=(32,32))
    kernel = forward_vz_kernel!(backend, block_size)
    kernel(kernel_params, field_params, ndrange = kernel_params.dim)
end;