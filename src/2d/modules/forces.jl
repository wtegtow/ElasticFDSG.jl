function apply_forces!(fdsg::FDSG2D, field_params, ti)
    # incremental stress 
    @unpack Mxx, Mxz, Mzz, sx, sz, stf_d1 = fdsg.source
    @unpack dx, dz = fdsg.domain
    @unpack dt = fdsg.time
    @unpack sxx, sxz, szz = field_params

    STF_d1_ti = stf_d1[ti]
    V = dx * dz

    sxx[sz,sx] -= dt / V * Mxx * STF_d1_ti
    szz[sz,sx] -= dt / V * Mzz * STF_d1_ti

    sxz[sz, sx  ] -= dt / (4V) * Mxz * STF_d1_ti 
    sxz[sz+1, sx  ] -= dt / (4V) * Mxz * STF_d1_ti   
    sxz[sz, sx+1] -= dt / (4V) * Mxz * STF_d1_ti   
    sxz[sz+1, sx+1] -= dt / (4V) * Mxz * STF_d1_ti   

end