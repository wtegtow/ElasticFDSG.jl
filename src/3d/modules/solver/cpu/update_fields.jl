# =============================================================================
#  CPU UPDATE STRESS FUNCTIONS 
# 1: sxx,syy,szz
# 2: sxy
# 3: sxz
# 4: syz
# =============================================================================

# 1 ========================= sxx syy szz =====================================

function update_normal_stresses!(
    vx, vy, vz,
    sxx, syy, szz,
    c11, c12, c13, c22, c23, c33,
    vx_x_old, vy_y_old, vz_z_old,
    b_x_odd, a_x_odd, K_x_odd,
    b_y_evn, a_y_evn, K_y_evn,
    b_z_evn, a_z_evn, K_z_evn,
    nz, ny, nx,
    dx, dy, dz,
    dt, c, N,
    pml_points_hash)
    
    Threads.@threads for z in 2+(N-1):nz-(N-1)
    for y in 2+(N-1):ny-(N-1)
    for x in 1+(N-1):nx-N
        
    
        # EFFECTIVE MEDIA PARAMETER
        c11eff = 0.5 * ( c11[z,y,x] + c11[z,y,x+1])
        c12eff = 0.5 * ( c12[z,y,x] + c12[z,y,x+1])
        c13eff = 0.5 * ( c13[z,y,x] + c13[z,y,x+1])
        c22eff = 0.5 * ( c22[z,y,x] + c22[z,y,x+1])
        c23eff = 0.5 * ( c23[z,y,x] + c23[z,y,x+1])
        c33eff = 0.5 * ( c33[z,y,x] + c33[z,y,x+1])
            
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
            vx_x_old[pml_idx] = b_x_odd[x] * vx_x_old[pml_idx] + a_x_odd[x] * vx_x
            vy_y_old[pml_idx] = b_y_evn[y] * vy_y_old[pml_idx] + a_y_evn[y] * vy_y
            vz_z_old[pml_idx] = b_z_evn[z] * vz_z_old[pml_idx] + a_z_evn[z] * vz_z
            
            # UPDATE SPACE DERIVATIVE
            vx_x = vx_x/K_x_odd[x] + vx_x_old[pml_idx]
            vy_y = vy_y/K_y_evn[y] + vy_y_old[pml_idx]
            vz_z = vz_z/K_z_evn[z] + vz_z_old[pml_idx]

        end
        
        # UPDATE STRESSES
        sxx[z,y,x] += dt * (c11eff * vx_x +  c12eff * vy_y + c13eff * vz_z)
        syy[z,y,x] += dt * (c12eff * vx_x +  c22eff * vy_y + c23eff * vz_z)
        szz[z,y,x] += dt * (c13eff * vx_x +  c23eff * vy_y + c33eff * vz_z)
        
    end
    end
    end
end;
    
    
# 2 =========================== sxy ===============================================
function update_sxy!(
    vx, vy,
    sxy,
    c66,
    vy_x_old, vx_y_old,
    b_x_evn, a_x_evn, K_x_evn,
    b_y_odd, a_y_odd, K_y_odd,
    nz, ny, nx,
    dx, dy,
    dt, c, N,
    pml_points_hash)
    
    Threads.@threads for z in 1:nz
    for y in 1+(N-1):ny-N
    for x in 2+(N-1):nx-(N-1)
    
        # EFFECTIVE MEDIA PARAMETER
        c66eff = 0.5 * ( c66[z,y,x] + c66[z,y+1,x])
        
        vy_x = 0
        vx_y = 0
        for i in 1:N
            vy_x +=  c[i]/dx * (vy[z,y,x+(i-1)]-vy[z,y,x-i])
            vx_y +=  c[i]/dy * (vx[z,y+i,x]-vx[z,y-(i-1),x])
        end
        
        
        # UPDATE PML VALUES
        if pml_points_hash[z,y,x] != -1 

            pml_idx = pml_points_hash[z,y,x]
            vy_x_old[pml_idx] = b_x_evn[x] * vy_x_old[pml_idx] + a_x_evn[x] * vy_x
            vx_y_old[pml_idx] = b_y_odd[y] * vx_y_old[pml_idx] + a_y_odd[y] * vx_y
            
            # UPDATE SPACE DERIVATIVE
            vy_x = vy_x/K_x_evn[x] + vy_x_old[pml_idx]
            vx_y = vx_y/K_y_odd[y] + vx_y_old[pml_idx]

        end
        
        # UPDATE STRESSES
        sxy[z,y,x] += dt * c66eff * (vy_x + vx_y)
        
    end
    end
    end
end;
    
    
# 3 =========================== sxz ===============================================
    
function update_sxz!(
    vx, vz,
    sxz,
    c55,
    vz_x_old, vx_z_old,
    b_x_evn, a_x_evn, K_x_evn,
    b_z_odd, a_z_odd, K_z_odd,
    nz, ny, nx,
    dx, dz,
    dt, c, N,
    pml_points_hash)
    
    Threads.@threads for z in 1+(N-1):nz-N
    for y in 1:ny
    for x in 2+(N-1):nx-(N-1)
    
        # EFFECTIVE MEDIA PARAMETER
        c55eff = 0.5 * ( c55[z,y,x] + c55[z+1,y,x])
        
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
            vz_x_old[pml_idx] = b_x_evn[x] * vz_x_old[pml_idx] + a_x_evn[x] * vz_x
            vx_z_old[pml_idx] = b_z_odd[z] * vx_z_old[pml_idx] + a_z_odd[z] * vx_z
            
            # UPDATE SPACE DERIVATIVE
            vz_x = vz_x/K_x_evn[x] + vz_x_old[pml_idx]
            vx_z = vx_z/K_z_odd[z] + vx_z_old[pml_idx] 

        end
        
        # UPDATE STRESSES
        sxz[z,y,x] += dt * c55eff * (vz_x + vx_z)
    
    end
    end
    end
end;
    
    
# 4 =========================== syz ===============================================
    
function update_syz!(
    vy, vz,
    syz,
    c44,
    vz_y_old, vy_z_old,
    b_y_odd, a_y_odd, K_y_odd,
    b_z_odd, a_z_odd, K_z_odd,
    nz, ny, nx,
    dy, dz,
    dt, c, N,
    pml_points_hash)
    
    Threads.@threads for z in 1+(N-1):nz-N
    for y in 1+(N-1):ny-N
    for x in 1:nx
    
        # EFFECTIVE MEDIA PARAMETER
        c44eff = 0.25 *(c44[z,y,x]     + c44[z+1,y,x]+
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
            vz_y_old[pml_idx] = b_y_odd[y] * vz_y_old[pml_idx] + a_y_odd[y] * vz_y
            vy_z_old[pml_idx] = b_z_odd[z] * vy_z_old[pml_idx] + a_z_odd[z] * vy_z
            
            # UPDATE SPACE DERIVATIVE
            vz_y = vz_y/K_y_odd[y] + vz_y_old[pml_idx]
            vy_z = vy_z/K_z_odd[z] + vy_z_old[pml_idx]

        end
        
        # UPDATE STRESS
        syz[z,y,x] += dt * c44eff * (vz_y + vy_z)
        
    end
    end
    end
end;
    
    
# =============================================================================
#  CPU UPDATE VELOCITY FUNCTIONS
# 1: vx
# 2: vy
# 3: vz
# =============================================================================

# 1 =========================== vx ===============================================
function update_vx!(
    vx,
    sxx, sxy, sxz,
    rho,
    sxx_x_old, sxy_y_old, sxz_z_old,
    b_x_evn, a_x_evn, K_x_evn,
    b_y_evn, a_y_evn, K_y_evn,
    b_z_evn, a_z_evn, K_z_evn,
    nz, ny, nx,
    dx, dy, dz,
    dt, c, N,
    pml_points_hash)
    
    Threads.@threads for z in 2+(N-1):nz-(N-1)
    for y in 2+(N-1):ny-(N-1)
    for x in 2+(N-1):nx-(N-1)
    
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
            sxx_x_old[pml_idx] = b_x_evn[x] * sxx_x_old[pml_idx] + a_x_evn[x] * sxx_x
            sxy_y_old[pml_idx] = b_y_evn[y] * sxy_y_old[pml_idx] + a_y_evn[y] * sxy_y
            sxz_z_old[pml_idx] = b_z_evn[z] * sxz_z_old[pml_idx] + a_z_evn[z] * sxz_z
            
            # UPDATE SPACE DERIVATIVE
            sxx_x = sxx_x/K_x_evn[x] + sxx_x_old[pml_idx]
            sxy_y = sxy_y/K_y_evn[y] + sxy_y_old[pml_idx]
            sxz_z = sxz_z/K_z_evn[z] + sxz_z_old[pml_idx]
        end
        
        
        # UPDATE VELOCITY
        vx[z,y,x] += dt/rho[z,y,x] * (sxx_x + sxy_y + sxz_z)
    
    end
    end
    end
end;
    
# 2 =========================== vy ===============================================
function update_vy!(
    vy,
    sxy, syy, syz,
    rho,
    sxy_x_old, syy_y_old, syz_z_old,
    b_x_odd, a_x_odd, K_x_odd,
    b_y_odd, a_y_odd, K_y_odd,
    b_z_evn, a_z_evn, K_z_evn,
    nz, ny, nx,
    dx, dy, dz,
    dt, c, N,
    pml_points_hash)
    
    Threads.@threads for z in 2+(N-1):nz-(N-1)
    for y in 1+(N-1):ny-N
    for x in 1+(N-1):nx-N
    
        # EFFECTIVE MEDIA PARAMETER
        rhoeff = 0.25 * (rho[z,y,x] + rho[z,y+1,x] + rho[z,y,x+1] + rho[z,y+1,x+1])
        
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
            sxy_x_old[pml_idx] = b_x_odd[x] * sxy_x_old[pml_idx] + a_x_odd[x] * sxy_x
            syy_y_old[pml_idx] = b_y_odd[y] * syy_y_old[pml_idx] + a_y_odd[y] * syy_y
            syz_z_old[pml_idx] = b_z_evn[z] * syz_z_old[pml_idx] + a_z_evn[z] * syz_z
            
            # UPDATE SPACE DERIVATIVE
            sxy_x = sxy_x/K_x_odd[x] + sxy_x_old[pml_idx]
            syy_y = syy_y/K_y_odd[y] + syy_y_old[pml_idx]
            syz_z = syz_z/K_z_evn[z] + syz_z_old[pml_idx]
        end
        
        # UPDATE VELOCITY
        vy[z,y,x] += dt/rhoeff * (sxy_x + syy_y + syz_z)
    
    end
    end
    end
end;
    
    
# 3 =========================== vz ===============================================
    
function update_vz!(
    vz,
    sxz, syz, szz,
    rho,
    sxz_x_old, syz_y_old, szz_z_old,
    b_x_odd, a_x_odd, K_x_odd,
    b_y_evn, a_y_evn, K_y_evn,
    b_z_odd, a_z_odd, K_z_odd,
    nz, ny, nx,
    dx, dy, dz,
    dt, c, N,
    pml_points_hash)
    
    Threads.@threads for z in 1+(N-1):nz-N
    for y in 2+(N-1):ny-(N-1)
    for x in 1+(N-1):nx-N
    
        # EFFECTIVE MEDIA PARAMETER
        rhoeff = 0.25 * (rho[z,y,x] + rho[z+1,y,x] + rho[z,y,x+1] + rho[z+1,y,x+1])
        
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
            sxz_x_old[pml_idx] = b_x_odd[x] * sxz_x_old[pml_idx] + a_x_odd[x] * sxz_x
            syz_y_old[pml_idx] = b_y_evn[y] * syz_y_old[pml_idx] + a_y_evn[y] * syz_y
            szz_z_old[pml_idx] = b_z_odd[z] * szz_z_old[pml_idx] + a_z_odd[z] * szz_z
            
            # UPDATE SPACE DERIVATIVE
            sxz_x = sxz_x/K_x_odd[x] + sxz_x_old[pml_idx]
            syz_y = syz_y/K_y_evn[y] + syz_y_old[pml_idx]
            szz_z = szz_z/K_z_odd[z] + szz_z_old[pml_idx]
        end
        
        # UPDATE VELOCITY
        vz[z,y,x] += dt/rhoeff * (sxz_x + syz_y + szz_z)
    
    end
    end
    end
end
    
    