function moment_tensor(fdsg3d::FDSG3D)

    """
    δ = Dip, λ = Rake, Φ = Strike
    Coordinate System: x -> north,
                       y -> east,
                       z -> positive downward
  
        x
      /
    /
    ---- y
    |
    |
    z
  
    """
    anisotropic_mtensor = fdsg3d.settings.config["source"]["double_couple"]["anisotropic_moment_tensor"]
    δ = fdsg3d.settings.config["source"]["double_couple"]["dip"]
    λ = fdsg3d.settings.config["source"]["double_couple"]["rake"]
    Φ = fdsg3d.settings.config["source"]["double_couple"]["strike"]

    δ = deg2rad(δ)
    λ = deg2rad(λ)
    Φ = deg2rad(Φ)
  
    Mxx = -(sin(δ)*cos(λ)*sin(2*Φ) + sin(2*δ)*sin(λ)*sin(Φ)^2)
    Mxy = sin(δ)*cos(λ)*cos(2*Φ) + 1/2*sin(2*δ)*sin(λ)*sin(2*Φ)
    Mxz = -(cos(δ)*cos(λ)*cos(Φ) + cos(2*δ)*sin(λ)*sin(Φ))
    Myy = sin(δ)*cos(λ)*sin(2*Φ) - sin(2*δ)*sin(λ)*cos(Φ)^2
    Myz = -(cos(δ)*cos(λ)*sin(Φ) - cos(2*δ)*sin(λ)*cos(Φ))
    Mzz = sin(2*δ)*sin(λ)
  
    if anisotropic_mtensor

        src = (fdsg3d.source.zid, fdsg3d.source.yid, fdsg3d.source.xid)

        C_source = [fdsg3d.elastic.c11[src...] fdsg3d.elastic.c12[src...] fdsg3d.elastic.c13[src...]    0                   0                  0;
                    fdsg3d.elastic.c12[src...] fdsg3d.elastic.c22[src...] fdsg3d.elastic.c23[src...]    0                   0                  0;
                    fdsg3d.elastic.c13[src...] fdsg3d.elastic.c23[src...] fdsg3d.elastic.c33[src...]    0                   0                  0;
                                0                       0                           0            fdsg3d.elastic.c44[src...] 0                  0;
                                0                       0                           0                   0           fdsg3d.elastic.c55[src...] 0;
                                0                       0                           0                   0                   0           fdsg3d.elastic.c66[src...]]
        
        M = [Mxx Myy Mzz Myz Mxz Mxy] * C_source 
        M = M / norm(M)
        Mxx = M[1]
        Myy = M[2]
        Mzz = M[3]
        Myz = M[4]
        Mxz = M[5]
        Mxy = M[6]

    end

    return fdsg3d.settings.float.([Mxx, Mxy, Mxz, Myy, Myz, Mzz])
end;


function init_forces(fdsg3d::FDSG3D) 
    
    use_point = fdsg3d.settings.config["source"]["point_source"]["use"]
    use_double = fdsg3d.settings.config["source"]["double_couple"]["use"]
    @assert use_point != use_double "Only one source type can be used at a time"

    if use_point
        phi = deg2rad(fdsg3d.settings.config["source"]["point_source"]["phi"])
        theta = deg2rad(fdsg3d.settings.config["source"]["point_source"]["theta"])

        force_x = sin(theta) .* cos(phi) .* fdsg3d.source.STF .* fdsg3d.time.dt/fdsg3d.source.rhosrc
        force_y = sin(theta) .* sin(phi) .* fdsg3d.source.STF .* fdsg3d.time.dt/fdsg3d.source.rhosrc
        force_z = cos(theta) .* fdsg3d.source.STF .* fdsg3d.time.dt/fdsg3d.source.rhosrc

        forces = [fdsg3d.settings.float.(force_x), 
                  fdsg3d.settings.float.(force_y), 
                  fdsg3d.settings.float.(force_z)]

    elseif use_double
        forces = moment_tensor(fdsg3d)
    end

    return forces
end;


function apply_forces!(forces,
    vx,vy,vz,
    sx,sy,sz,
    dx,dy,dz,
    STF,rhosrc,dt,ti)

    if length(forces) == 3 # point source
        force_x, force_y, force_z = forces
        vx[sz,sy,sx] += force_x[ti]
        vy[sz,sy,sx] += force_y[ti]
        vz[sz,sy,sx] += force_z[ti]

    elseif length(forces) == 6 # double couple source
        Mxx, Mxy, Mxz, Myy, Myz, Mzz = forces
        STF_ti = STF[ti]

         # Force X
        # Mxx
        vx[sz,sy,sx+1] += dt/rhosrc * (Mxx * STF_ti)/(dx^2 * dy * dz)
        vx[sz,sy,sx]   -= dt/rhosrc * (Mxx * STF_ti)/(dx^2 * dy * dz)

        # Mxy 
        vx[sz,sy+1,sx+1] += dt/rhosrc * (Mxy * STF_ti)/ (4*dx * dy^2 * dz)
        vx[sz,sy+1,sx]   += dt/rhosrc * (Mxy * STF_ti)/ (4*dx * dy^2 * dz)
        vx[sz,sy-1,sx+1] -= dt/rhosrc * (Mxy * STF_ti)/ (4*dx * dy^2 * dz)
        vx[sz,sy-1,sx]   -= dt/rhosrc * (Mxy * STF_ti)/ (4*dx * dy^2 * dz)

        # Mxz
        vx[sz+1,sy,sx+1] += dt/rhosrc * (Mxz * STF_ti)/ (4*dx * dy * dz^2)
        vx[sz+1,sy,sx]   += dt/rhosrc * (Mxz * STF_ti)/ (4*dx * dy * dz^2)
        vx[sz-1,sy,sx+1] -= dt/rhosrc * (Mxz * STF_ti)/ (4*dx * dy * dz^2)
        vx[sz-1,sy,sx]   -= dt/rhosrc * (Mxz * STF_ti)/ (4*dx * dy * dz^2)

        # Force Y
        # Myy
        vy[sz,sy+1,sx] += dt/rhosrc * (Myy * STF_ti)/(dx * dy^2 * dz)
        vy[sz,sy,sx]   -= dt/rhosrc * (Myy * STF_ti)/(dx * dy^2 * dz)

        # Mxy 
        vy[sz,sy+1,sx+1] += dt/rhosrc * (Mxy * STF_ti)/ (dx^2 * 4*dy * dz)
        vy[sz,sy,sx+1]   += dt/rhosrc * (Mxy * STF_ti)/ (dx^2 * 4*dy * dz)
        vy[sz,sy+1,sx-1] -= dt/rhosrc * (Mxy * STF_ti)/ (dx^2 * 4*dy * dz)
        vy[sz,sy,sx-1]   -= dt/rhosrc * (Mxy * STF_ti)/ (dx^2 * 4*dy * dz)

        # Myz
        vy[sz+1,sy+1,sx] += dt/rhosrc * (Myz * STF_ti)/ (dx * 4*dy * dz^2)
        vy[sz+1,sy,sx]   += dt/rhosrc * (Myz * STF_ti)/ (dx * 4*dy * dz^2)
        vy[sz-1,sy+1,sx] -= dt/rhosrc * (Myz * STF_ti)/ (dx * 4*dy * dz^2)
        vy[sz-1,sy,sx]   -= dt/rhosrc * (Myz * STF_ti)/ (dx * 4*dy * dz^2)

        # Force Z
        # Mzz
        vz[sz+1,sy,sx] += dt/rhosrc * (Mzz * STF_ti)/(dx * dy * dz^2)
        vz[sz,sy,sx]   -= dt/rhosrc * (Mzz * STF_ti)/(dx * dy * dz^2)

        # Mxz
        vz[sz+1,sy,sx+1] += dt/rhosrc * (Mxz * STF_ti)/ (dx^2 * dy * 4*dz)
        vz[sz,sy,sx+1]   += dt/rhosrc * (Mxz * STF_ti)/ (dx^2 * dy * 4*dz)
        vz[sz+1,sy,sx-1] -= dt/rhosrc * (Mxz * STF_ti)/ (dx^2 * dy * 4*dz)
        vz[sz,sy,sx-1]   -= dt/rhosrc * (Mxz * STF_ti)/ (dx^2 * dy * 4*dz)

        # Myz
        vz[sz+1,sy+1,sx] += dt/rhosrc * (Myz * STF_ti)/ (dx * dy^2 * 4*dz)
        vz[sz,sy+1,sx]   += dt/rhosrc * (Myz * STF_ti)/ (dx * dy^2 * 4*dz)
        vz[sz+1,sy-1,sx] -= dt/rhosrc * (Myz * STF_ti)/ (dx * dy^2 * 4*dz)
        vz[sz,sy-1,sx]   -= dt/rhosrc * (Myz * STF_ti)/ (dx * dy^2 * 4*dz)

    end

end
