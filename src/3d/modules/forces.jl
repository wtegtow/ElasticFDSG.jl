function apply_forces!(fdsg::FDSG3D, field_params, ti)
    # incremental stress 
    
    @unpack Mxx, Mxy, Mxz, Myy, Myz, Mzz, sx, sy, sz, stf_d1 = fdsg.source
    @unpack dx, dy, dz = fdsg.domain
    @unpack dt = fdsg.time
    @unpack sxx, sxy, sxz, syy, syz, szz = field_params

    STF_d1_ti = stf_d1[ti]
    V = dx * dy * dz

    sxx[sz,sy,sx] -= dt / V * Mxx * STF_d1_ti
    syy[sz,sy,sx] -= dt / V * Myy * STF_d1_ti
    szz[sz,sy,sx] -= dt / V * Mzz * STF_d1_ti
   
    sxy[sz, sy,   sx  ] -= dt / (4V) * Mxy * STF_d1_ti   
    sxy[sz, sy+1, sx  ] -= dt / (4V) * Mxy * STF_d1_ti  
    sxy[sz, sy,   sx+1] -= dt / (4V) * Mxy * STF_d1_ti   
    sxy[sz, sy+1, sx+1] -= dt / (4V) * Mxy * STF_d1_ti   

    sxz[sz,   sy, sx  ] -= dt / (4V) * Mxz * STF_d1_ti 
    sxz[sz+1, sy, sx  ] -= dt / (4V) * Mxz * STF_d1_ti   
    sxz[sz,   sy, sx+1] -= dt / (4V) * Mxz * STF_d1_ti   
    sxz[sz+1, sy, sx+1] -= dt / (4V) * Mxz * STF_d1_ti   

    syz[sz,   sy,   sx] -= dt / (4V) * Myz * STF_d1_ti 
    syz[sz+1, sy,   sx] -= dt / (4V) * Myz * STF_d1_ti  
    syz[sz,   sy+1, sx] -= dt / (4V) * Myz * STF_d1_ti  
    syz[sz+1, sy+1, sx] -= dt / (4V) * Myz * STF_d1_ti 



end