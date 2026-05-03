function _stress_glut_source!(fields::Fields2D, source::Source2D,
                         domain::Domain{2}, time::SimTime, ti::Int)

    sx = source.sx;  sz = source.sz
    dx = T_spacing(domain.coordinates[1])
    dz = T_spacing(domain.coordinates[2])
    dt = time.dt
    V  = dx * dz
    S  = source.stf_d1[ti]

    fields.sxx[sx,   sz  ] -= dt / V       * source.Mxx * S
    fields.szz[sx,   sz  ] -= dt / V       * source.Mzz * S

    fields.sxz[sx,   sz  ] -= dt / (4*V)   * source.Mxz * S
    fields.sxz[sx+1, sz  ] -= dt / (4*V)   * source.Mxz * S
    fields.sxz[sx,   sz+1] -= dt / (4*V)   * source.Mxz * S
    fields.sxz[sx+1, sz+1] -= dt / (4*V)   * source.Mxz * S
end

function _stress_glut_source!(fields::Fields3D, source::Source3D,
                         domain::Domain{3}, time::SimTime, ti::Int)

    sx = source.sx;  sy = source.sy;  sz = source.sz
    dx = T_spacing(domain.coordinates[1])
    dy = T_spacing(domain.coordinates[2])
    dz = T_spacing(domain.coordinates[3])
    dt = time.dt
    V  = dx * dy * dz
    S  = source.stf_d1[ti]

    fields.sxx[sx,   sy,   sz  ] -= dt / V       * source.Mxx * S
    fields.syy[sx,   sy,   sz  ] -= dt / V       * source.Myy * S
    fields.szz[sx,   sy,   sz  ] -= dt / V       * source.Mzz * S

    fields.sxy[sx,   sy,   sz  ] -= dt / (4*V)   * source.Mxy * S
    fields.sxy[sx+1, sy,   sz  ] -= dt / (4*V)   * source.Mxy * S
    fields.sxy[sx,   sy+1, sz  ] -= dt / (4*V)   * source.Mxy * S
    fields.sxy[sx+1, sy+1, sz  ] -= dt / (4*V)   * source.Mxy * S

    fields.sxz[sx,   sy,   sz  ] -= dt / (4*V)   * source.Mxz * S
    fields.sxz[sx+1, sy,   sz  ] -= dt / (4*V)   * source.Mxz * S
    fields.sxz[sx,   sy,   sz+1] -= dt / (4*V)   * source.Mxz * S
    fields.sxz[sx+1, sy,   sz+1] -= dt / (4*V)   * source.Mxz * S

    fields.syz[sx,   sy,   sz  ] -= dt / (4*V)   * source.Myz * S
    fields.syz[sx,   sy,   sz+1] -= dt / (4*V)   * source.Myz * S
    fields.syz[sx,   sy+1, sz  ] -= dt / (4*V)   * source.Myz * S
    fields.syz[sx,   sy+1, sz+1] -= dt / (4*V)   * source.Myz * S
end
