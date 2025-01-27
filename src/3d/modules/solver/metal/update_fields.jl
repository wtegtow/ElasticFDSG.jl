# =============================================================================
#  METAL UPDATE STRESS KERNEL 
# 1: sxx,syy,szz
# 2: sxy
# 3: sxz
# 4: syz
# =============================================================================

# 1 ========================= sxx syy szz =====================================

function update_normal_stresses_mkernel(
    vx, vy, vz,
    sxx, syy, szz,
    c11, c12, c13, c22, c23, c33,
    vx_x_old, vy_y_old, vz_z_old,
    x_odd,
    y_evn, 
    z_evn, 
    nz, ny, nx,
    dx, dy, dz,
    dt, c, N,
    pml_points_hash)

    z,y,x  = Tuple(thread_position_in_grid_3d())
    stride = Tuple(threads_per_grid_3d())
    
    while 2+(N-1) <= z <= nz-(N-1) &&
          2+(N-1) <= y <= ny-(N-1) &&
          1+(N-1) <= x <= nx-N
        
        # EFFECTIVE MEDIA PARAMETER
        c11eff = 0.5f0 * ( c11[z,y,x] + c11[z,y,x+1])
        c12eff = 0.5f0 * ( c12[z,y,x] + c12[z,y,x+1])
        c13eff = 0.5f0 * ( c13[z,y,x] + c13[z,y,x+1])
        c22eff = 0.5f0 * ( c22[z,y,x] + c22[z,y,x+1])
        c23eff = 0.5f0 * ( c23[z,y,x] + c23[z,y,x+1])
        c33eff = 0.5f0 * ( c33[z,y,x] + c33[z,y,x+1])
            
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
        if pml_points_hash[z,y,x] != -1 

            pml_idx = pml_points_hash[z,y,x]
            vx_x_old[pml_idx] = x_odd[2,x] * vx_x_old[pml_idx] + x_odd[1,x] * vx_x
            vy_y_old[pml_idx] = y_evn[2,y] * vy_y_old[pml_idx] + y_evn[1,y] * vy_y
            vz_z_old[pml_idx] = z_evn[2,z] * vz_z_old[pml_idx] + z_evn[1,z] * vz_z
            
            # UPDATE SPACE DERIVATIVE
            vx_x = vx_x/x_odd[3,x] + vx_x_old[pml_idx]
            vy_y = vy_y/y_evn[3,y] + vy_y_old[pml_idx]
            vz_z = vz_z/z_evn[3,z] + vz_z_old[pml_idx]

        end
        
        # UPDATE STRESSES
        sxx[z,y,x] += dt * (c11eff * vx_x +  c12eff * vy_y + c13eff * vz_z)
        syy[z,y,x] += dt * (c12eff * vx_x +  c22eff * vy_y + c23eff * vz_z)
        szz[z,y,x] += dt * (c13eff * vx_x +  c23eff * vy_y + c33eff * vz_z)
    
        z += stride[1]
        y += stride[2]
        x += stride[3]
    end
end;



    
    
# 2 =========================== sxy ===============================================
function update_sxy_mkernel(
    vx, vy,
    sxy,
    c66,
    vy_x_old, vx_y_old,
    x_evn, 
    y_odd,
    nz, ny, nx,
    dx, dy,
    dt, c, N,
    pml_points_hash)

    z,y,x  = Tuple(thread_position_in_grid_3d())
    stride = Tuple(threads_per_grid_3d())
    
    while      1  <= z <= nz      &&
          1+(N-1) <= y <= ny-N    &&
          2+(N-1) <= x <= nx-(N-1)
    
        # EFFECTIVE MEDIA PARAMETER
        c66eff = 0.5f0 * ( c66[z,y,x] + c66[z,y+1,x])
        
        vy_x = 0
        vx_y = 0
        for i in 1:N
            vy_x +=  c[i]/dx * (vy[z,y,x+(i-1)]-vy[z,y,x-i])
            vx_y +=  c[i]/dy * (vx[z,y+i,x]-vx[z,y-(i-1),x])
        end
        
        # UPDATE PML VALUES
        if pml_points_hash[z,y,x] != -1 

            pml_idx = pml_points_hash[z,y,x]
            vy_x_old[pml_idx] = x_evn[2,x] * vy_x_old[pml_idx] + x_evn[1,x] * vy_x
            vx_y_old[pml_idx] = y_odd[2,y] * vx_y_old[pml_idx] + y_odd[1,y] * vx_y
            
            # UPDATE SPACE DERIVATIVE
            vy_x = vy_x/x_evn[3,x] + vy_x_old[pml_idx]
            vx_y = vx_y/y_odd[3,y] + vx_y_old[pml_idx]

        end
        
        # UPDATE STRESSES
        sxy[z,y,x] += dt * c66eff * (vy_x + vx_y)
        
        z += stride[1]
        y += stride[2]
        x += stride[3]
    
    end
end;
    
    
# 3 =========================== sxz ===============================================
    
function update_sxz_mkernel(
    vx, vz,
    sxz,
    c55,
    vz_x_old, vx_z_old,
    x_evn, 
    z_odd,
    nz, ny, nx,
    dx, dz,
    dt, c, N,
    pml_points_hash)
    
    z,y,x  = Tuple(thread_position_in_grid_3d())
    stride = Tuple(threads_per_grid_3d())

    while 1+(N-1) <= z <= nz-N     &&
          1       <= y <= ny       && 
          2+(N-1) <= x <= nx-(N-1)
    
        # EFFECTIVE MEDIA PARAMETER
        c55eff = 0.5f0 * ( c55[z,y,x] + c55[z+1,y,x])
        
        # SPACE DERIVATIVE
        vz_x = 0
        vx_z = 0
        for i in 1:N
            vz_x += c[i]/dx * (vz[z,y,x+(i-1)]-vz[z,y,x-i])
            vx_z += c[i]/dz * (vx[z+i,y,x]-vx[z-(i-1),y,x])
        end
        
        # UPDATE PML VALUES
        if pml_points_hash[z,y,x] != -1 

            pml_idx = pml_points_hash[z,y,x]
            vz_x_old[pml_idx] = x_evn[2,x] * vz_x_old[pml_idx] + x_evn[1,x] * vz_x
            vx_z_old[pml_idx] = z_odd[2,z] * vx_z_old[pml_idx] + z_odd[1,z] * vx_z
            
            # UPDATE SPACE DERIVATIVE
            vz_x = vz_x/x_evn[3,x] + vz_x_old[pml_idx]
            vx_z = vx_z/z_odd[3,z] + vx_z_old[pml_idx] 

        end
        
        # UPDATE STRESSES
        sxz[z,y,x] += dt * c55eff * (vz_x + vx_z)
    
        z += stride[1]
        y += stride[2]
        x += stride[3]
    end
end;
    
    
# 4 =========================== syz ===============================================
    
function update_syz_mkernel(
    vy, vz,
    syz,
    c44,
    vz_y_old, vy_z_old,
    y_odd,
    z_odd,
    nz, ny, nx,
    dy, dz,
    dt, c, N,
    pml_points_hash)

    z,y,x  = Tuple(thread_position_in_grid_3d())
    stride = Tuple(threads_per_grid_3d())
    
    while 1+(N-1) <= z <= nz-N  &&
          1+(N-1) <= y <= ny-N  && 
          1       <= x <= nx
    
        # EFFECTIVE MEDIA PARAMETER
        c44eff = 0.25f0 *(c44[z,y,x]     + c44[z+1,y,x]+
                          c44[z,y+1,x]   + c44[z+1,y+1,x])
        
        # SPACE DERIVATIVE
        vz_y = 0
        vy_z = 0
        for i in 1:N 
            vz_y += c[i]/dy * (vz[z,y+i,x]-vz[z,y-(i-1),x])
            vy_z += c[i]/dz * (vy[z+i,y,x]-vy[z-(i-1),y,x])
        end
        
        # UPDATE PML VALUES
        if pml_points_hash[z,y,x] != -1 

            pml_idx = pml_points_hash[z,y,x] 
            vz_y_old[pml_idx] = y_odd[2,y] * vz_y_old[pml_idx] + y_odd[1,y] * vz_y
            vy_z_old[pml_idx] = z_odd[2,z] * vy_z_old[pml_idx] + z_odd[1,z] * vy_z
            
            # UPDATE SPACE DERIVATIVE
            vz_y = vz_y/y_odd[3,y] + vz_y_old[pml_idx]
            vy_z = vy_z/z_odd[3,z] + vy_z_old[pml_idx]

        end
        
        # UPDATE STRESS
        syz[z,y,x] += dt * c44eff * (vz_y + vy_z)
        
        z += stride[1]
        y += stride[2]
        x += stride[3]

    end
end;
    
    
# =============================================================================
#  CPU UPDATE VELOCITY FUNCTIONS
# 1: vx
# 2: vy
# 3: vz
# =============================================================================

# 1 =========================== vx ===============================================
function update_vx_mkernel(
    vx,
    sxx, sxy, sxz,
    rho,
    sxx_x_old, sxy_y_old, sxz_z_old,
    x_evn, 
    y_evn,
    z_evn,
    nz, ny, nx,
    dx, dy, dz,
    dt, c, N,
    pml_points_hash)

    z,y,x  = Tuple(thread_position_in_grid_3d())
    stride = Tuple(threads_per_grid_3d())
    
    while 2+(N-1) <= z <= nz-(N-1) && 
          2+(N-1) <= y <= ny-(N-1) && 
          2+(N-1) <= x <= nx-(N-1)
    
        # EFFECTIVE MEDIA PARAMETER
        # rho defined on vx grid points :)
        
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
        if pml_points_hash[z,y,x] != -1 

            pml_idx = pml_points_hash[z,y,x] 
            sxx_x_old[pml_idx] = x_evn[2,x] * sxx_x_old[pml_idx] + x_evn[1,x] * sxx_x
            sxy_y_old[pml_idx] = y_evn[2,y] * sxy_y_old[pml_idx] + y_evn[1,y] * sxy_y
            sxz_z_old[pml_idx] = z_evn[2,z] * sxz_z_old[pml_idx] + z_evn[1,z] * sxz_z
            
            # UPDATE SPACE DERIVATIVE
            sxx_x = sxx_x/x_evn[3,x] + sxx_x_old[pml_idx]
            sxy_y = sxy_y/y_evn[3,y] + sxy_y_old[pml_idx]
            sxz_z = sxz_z/z_evn[3,z] + sxz_z_old[pml_idx]
        end
        
        # UPDATE VELOCITY
        vx[z,y,x] += dt/rho[z,y,x] * (sxx_x + sxy_y + sxz_z)
    
        z += stride[1]
        y += stride[2]
        x += stride[3]

    end
end;
    
# 2 =========================== vy ===============================================
function update_vy_mkernel(
    vy,
    sxy, syy, syz,
    rho,
    sxy_x_old, syy_y_old, syz_z_old,
    x_odd,
    y_odd, 
    z_evn,
    nz, ny, nx,
    dx, dy, dz,
    dt, c, N,
    pml_points_hash)

    z,y,x  = Tuple(thread_position_in_grid_3d())
    stride = Tuple(threads_per_grid_3d())
    
    while 2+(N-1) <= z <= nz-(N-1)  &&
          1+(N-1) <= y <= ny-N      &&
          1+(N-1) <= x <= nx-N
    
        # EFFECTIVE MEDIA PARAMETER
        rhoeff = 0.25f0 * (rho[z,y,x] + rho[z,y+1,x] + rho[z,y,x+1] + rho[z,y+1,x+1])
        
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
        if pml_points_hash[z,y,x] != -1 

            pml_idx = pml_points_hash[z,y,x] 
            sxy_x_old[pml_idx] = x_odd[2,x] * sxy_x_old[pml_idx] + x_odd[1,x] * sxy_x
            syy_y_old[pml_idx] = y_odd[2,y] * syy_y_old[pml_idx] + y_odd[1,y] * syy_y
            syz_z_old[pml_idx] = z_evn[2,z] * syz_z_old[pml_idx] + z_evn[1,z] * syz_z
            
            # UPDATE SPACE DERIVATIVE
            sxy_x = sxy_x/x_odd[3,x] + sxy_x_old[pml_idx]
            syy_y = syy_y/y_odd[3,y] + syy_y_old[pml_idx]
            syz_z = syz_z/z_evn[3,z] + syz_z_old[pml_idx]
        end
        
        # UPDATE VELOCITY
        vy[z,y,x] += dt/rhoeff * (sxy_x + syy_y + syz_z)
    
        z += stride[1]
        y += stride[2]
        x += stride[3]
    end
end;
    
    
# 3 =========================== vz ===============================================
    
function update_vz_mkernel(
    vz,
    sxz, syz, szz,
    rho,
    sxz_x_old, syz_y_old, szz_z_old,
    x_odd,
    y_evn,
    z_odd,
    nz, ny, nx,
    dx, dy, dz,
    dt, c, N,
    pml_points_hash)

    z,y,x  = Tuple(thread_position_in_grid_3d())
    stride = Tuple(threads_per_grid_3d())
    
    while 1+(N-1) <= z <= nz-N      && 
          2+(N-1) <= y <= ny-(N-1)  &&
          1+(N-1) <= x <= nx-N
    
        # EFFECTIVE MEDIA PARAMETER
        rhoeff = 0.25f0 * (rho[z,y,x] + rho[z+1,y,x] + rho[z,y,x+1] + rho[z+1,y,x+1])
        
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
        if pml_points_hash[z,y,x] != -1 

            pml_idx = pml_points_hash[z,y,x] 
            sxz_x_old[pml_idx] = x_odd[2,x] * sxz_x_old[pml_idx] + x_odd[1,x] * sxz_x
            syz_y_old[pml_idx] = y_evn[2,y] * syz_y_old[pml_idx] + y_evn[1,y] * syz_y
            szz_z_old[pml_idx] = z_odd[2,z] * szz_z_old[pml_idx] + z_odd[1,z] * szz_z
            
            # UPDATE SPACE DERIVATIVE
            sxz_x = sxz_x/x_odd[3,x] + sxz_x_old[pml_idx]
            syz_y = syz_y/y_evn[3,y] + syz_y_old[pml_idx]
            szz_z = szz_z/z_odd[3,z] + szz_z_old[pml_idx]
        end
        
        # UPDATE VELOCITY
        vz[z,y,x] += dt/rhoeff * (sxz_x + syz_y + szz_z)
    
        z += stride[1]
        y += stride[2]
        x += stride[3]
    end
end
    