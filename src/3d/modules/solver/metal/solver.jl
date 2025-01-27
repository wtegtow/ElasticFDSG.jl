# =============================================================================
# metal solver
# =============================================================================

include(joinpath(@__DIR__,"update_fields.jl"))

function solver_metal(fdsg3d::FDSG3D)

    # switch everything needed to gpu
    ARRAY = MtlArray

    vx = ARRAY(fdsg3d.fields.vx)
    vy = ARRAY(fdsg3d.fields.vy)
    vz = ARRAY(fdsg3d.fields.vz)

    sxx = ARRAY(fdsg3d.fields.sxx)
    sxy = ARRAY(fdsg3d.fields.sxy)
    sxz = ARRAY(fdsg3d.fields.sxz)
    syy = ARRAY(fdsg3d.fields.syy)
    syz = ARRAY(fdsg3d.fields.syz)
    szz = ARRAY(fdsg3d.fields.szz)

    c11 = ARRAY(fdsg3d.elastic.c11)
    c12 = ARRAY(fdsg3d.elastic.c12)
    c13 = ARRAY(fdsg3d.elastic.c13)
    c22 = ARRAY(fdsg3d.elastic.c22)
    c23 = ARRAY(fdsg3d.elastic.c23)
    c33 = ARRAY(fdsg3d.elastic.c33)
    c44 = ARRAY(fdsg3d.elastic.c44)
    c55 = ARRAY(fdsg3d.elastic.c55)
    c66 = ARRAY(fdsg3d.elastic.c66)
    rho = ARRAY(fdsg3d.elastic.rho)

    vx_x_old = ARRAY(fdsg3d.pml.vx_x_old)
    vy_y_old = ARRAY(fdsg3d.pml.vy_y_old)
    vz_z_old = ARRAY(fdsg3d.pml.vz_z_old)
    vy_x_old = ARRAY(fdsg3d.pml.vy_x_old)
    vx_y_old = ARRAY(fdsg3d.pml.vx_y_old)
    vz_x_old = ARRAY(fdsg3d.pml.vz_x_old)
    vx_z_old = ARRAY(fdsg3d.pml.vx_z_old)
    vz_y_old = ARRAY(fdsg3d.pml.vz_y_old)
    vy_z_old = ARRAY(fdsg3d.pml.vy_z_old)

    sxx_x_old = ARRAY(fdsg3d.pml.sxx_x_old)
    sxy_y_old = ARRAY(fdsg3d.pml.sxy_y_old)
    sxz_z_old = ARRAY(fdsg3d.pml.sxz_z_old)
    sxy_x_old = ARRAY(fdsg3d.pml.sxy_x_old)
    syy_y_old = ARRAY(fdsg3d.pml.syy_y_old)
    syz_z_old = ARRAY(fdsg3d.pml.syz_z_old)
    sxz_x_old = ARRAY(fdsg3d.pml.sxz_x_old)
    syz_y_old = ARRAY(fdsg3d.pml.syz_y_old)
    szz_z_old = ARRAY(fdsg3d.pml.szz_z_old)

    nz = fdsg3d.domain.nz
    ny = fdsg3d.domain.ny
    nx = fdsg3d.domain.nx
    dx = fdsg3d.domain.dx
    dy = fdsg3d.domain.dy
    dz = fdsg3d.domain.dz
    dt = fdsg3d.time.dt
    c = ARRAY(fdsg3d.settings.c)
    N = fdsg3d.settings.N
    dim = fdsg3d.domain.dim

    # because metal kernel only allow 31 arguments, we combine the damping profiles into single arrays
    x_evn = Metal.zeros(3,nx)
    x_evn[1,:] .= ARRAY(fdsg3d.pml.a_x_evn)
    x_evn[2,:] .= ARRAY(fdsg3d.pml.b_x_evn)
    x_evn[3,:] .= ARRAY(fdsg3d.pml.K_x_evn)

    x_odd = Metal.zeros(3,nx)
    x_odd[1,:] .= ARRAY(fdsg3d.pml.a_x_odd)
    x_odd[2,:] .= ARRAY(fdsg3d.pml.b_x_odd)
    x_odd[3,:] .= ARRAY(fdsg3d.pml.K_x_odd)

    y_evn = Metal.zeros(3, ny)
    y_evn[1, :] .= ARRAY(fdsg3d.pml.a_y_evn)
    y_evn[2, :] .= ARRAY(fdsg3d.pml.b_y_evn)
    y_evn[3, :] .= ARRAY(fdsg3d.pml.K_y_evn)

    y_odd = Metal.zeros(3, ny)
    y_odd[1, :] .= ARRAY(fdsg3d.pml.a_y_odd)
    y_odd[2, :] .= ARRAY(fdsg3d.pml.b_y_odd)
    y_odd[3, :] .= ARRAY(fdsg3d.pml.K_y_odd)

    z_evn = Metal.zeros(3, nz)
    z_evn[1, :] .= ARRAY(fdsg3d.pml.a_z_evn)
    z_evn[2, :] .= ARRAY(fdsg3d.pml.b_z_evn)
    z_evn[3, :] .= ARRAY(fdsg3d.pml.K_z_evn)

    z_odd = Metal.zeros(3, nz)
    z_odd[1, :] .= ARRAY(fdsg3d.pml.a_z_odd)
    z_odd[2, :] .= ARRAY(fdsg3d.pml.b_z_odd)
    z_odd[3, :] .= ARRAY(fdsg3d.pml.K_z_odd)

    pml_points_hash = ARRAY(fdsg3d.domain.pml_points_hash)

    # compile kernel
    # stresses
    update_normal_stresses_metal = @metal launch=false update_normal_stresses_mkernel(
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

    update_sxy_metal = @metal launch=false update_sxy_mkernel(
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

    update_sxz_metal = @metal launch=false update_sxz_mkernel(
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

    update_syz_metal = @metal launch=false update_syz_mkernel(
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

    # velocites 
    update_vx_metal = @metal launch=false update_vx_mkernel(
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

    update_vy_metal = @metal launch=false update_vy_mkernel(
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

    update_vz_metal = @metal launch=false update_vz_mkernel(
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


    dim_arg_sort = sort(collect(dim),rev=true)
    w = min(size(dim_arg_sort, 1), update_vz_metal.pipeline.threadExecutionWidth)
    h = min(size(dim_arg_sort, 2), update_vz_metal.pipeline.threadExecutionWidth,
                                   update_vz_metal.pipeline.maxTotalThreadsPerThreadgroup ÷ w)
    d = min(size(dim_arg_sort, 3), update_vz_metal.pipeline.maxTotalThreadsPerThreadgroup ÷ (w*h))

    threads = (w, h, d)
    groups = cld.(dim, threads)

    # forces 
    forces = init_forces(fdsg3d)

    # progressbar 
    if fdsg3d.settings.config["settings"]["show_progress_in_console"] 
        p = Progress(fdsg3d.time.nt; showspeed=true)
    end 


    # time loop 
    for ti in 1:fdsg3d.time.nt 

        # save receiver
        @Metal.allowscalar save_geophones!(fdsg3d.geophones, vx, vy, vz, ti)
        save_das!(fdsg3d.das, sxx, syy, szz, ti, fdsg3d.settings.float)
        save_snapshots!(fdsg3d.snapshots, vx, vy, vz, sxx, sxy, sxz, syy, syz, szz, ti, fdsg3d.settings.float)

        # forward
        update_normal_stresses_metal(
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
            pml_points_hash, threads=threads, groups=groups)
        
        update_sxy_metal(
            vx, vy,
            sxy,
            c66,
            vy_x_old, vx_y_old,
            x_evn, 
            y_odd,
            nz, ny, nx,
            dx, dy,
            dt, c, N,
            pml_points_hash, threads=threads, groups=groups)
        
        update_sxz_metal(
                vx, vz,
                sxz,
                c55,
                vz_x_old, vx_z_old,
                x_evn, 
                z_odd,
                nz, ny, nx,
                dx, dz,
                dt, c, N,
                pml_points_hash, threads=threads, groups=groups)
        
        update_syz_metal(
                    vy, vz,
                    syz,
                    c44,
                    vz_y_old, vy_z_old,
                    y_odd,
                    z_odd,
                    nz, ny, nx,
                    dy, dz,
                    dt, c, N,
                    pml_points_hash, threads=threads, groups=groups)
        
        # velocites 
        update_vx_metal(
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
            pml_points_hash, threads=threads, groups=groups)
        
        update_vy_metal(
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
                pml_points_hash, threads=threads, groups=groups)
        
        update_vz_metal(
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
                    pml_points_hash, threads=threads, groups=groups)

        @Metal.allowscalar apply_forces!(forces, 
                            vx, vy, vz,
                            fdsg3d.source.xid, fdsg3d.source.yid, fdsg3d.source.zid,
                            dx, dy, dz,
                            fdsg3d.source.STF, fdsg3d.source.rhosrc, dt, 
                            ti)

        # progress bar
        if fdsg3d.settings.config["settings"]["show_progress_in_console"] next!(p) end

        # garbage collection every 100 iterations
        if ti % 100 == 0
            GC.gc()
        end
    end
end;
