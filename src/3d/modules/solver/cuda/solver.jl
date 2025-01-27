# =============================================================================
# cuda solver
# =============================================================================

include(joinpath(@__DIR__,"update_fields.jl"))

function solver_cuda(fdsg3d::FDSG3D)

    # switch everything needed to gpu
    ARRAY = CuArray

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
    # for convinience, i use the same logic for cuda kernels
    x_evn = CUDA.zeros(3,nx)
    x_evn[1,:] .= ARRAY(fdsg3d.pml.a_x_evn)
    x_evn[2,:] .= ARRAY(fdsg3d.pml.b_x_evn)
    x_evn[3,:] .= ARRAY(fdsg3d.pml.K_x_evn)

    x_odd = CUDA.zeros(3,nx)
    x_odd[1,:] .= ARRAY(fdsg3d.pml.a_x_odd)
    x_odd[2,:] .= ARRAY(fdsg3d.pml.b_x_odd)
    x_odd[3,:] .= ARRAY(fdsg3d.pml.K_x_odd)

    y_evn = CUDA.zeros(3, ny)
    y_evn[1, :] .= ARRAY(fdsg3d.pml.a_y_evn)
    y_evn[2, :] .= ARRAY(fdsg3d.pml.b_y_evn)
    y_evn[3, :] .= ARRAY(fdsg3d.pml.K_y_evn)

    y_odd = CUDA.zeros(3, ny)
    y_odd[1, :] .= ARRAY(fdsg3d.pml.a_y_odd)
    y_odd[2, :] .= ARRAY(fdsg3d.pml.b_y_odd)
    y_odd[3, :] .= ARRAY(fdsg3d.pml.K_y_odd)

    z_evn = CUDA.zeros(3, nz)
    z_evn[1, :] .= ARRAY(fdsg3d.pml.a_z_evn)
    z_evn[2, :] .= ARRAY(fdsg3d.pml.b_z_evn)
    z_evn[3, :] .= ARRAY(fdsg3d.pml.K_z_evn)

    z_odd = CUDA.zeros(3, nz)
    z_odd[1, :] .= ARRAY(fdsg3d.pml.a_z_odd)
    z_odd[2, :] .= ARRAY(fdsg3d.pml.b_z_odd)
    z_odd[3, :] .= ARRAY(fdsg3d.pml.K_z_odd)

    pml_points_hash = ARRAY(fdsg3d.domain.pml_points_hash)

    THREADS_PER_BLOCK = (7, 7, 7)
    BLOCKS_PER_GRID   = (ceil(Int, nz/7), ceil(Int, ny/7), ceil(Int, nx/7))

    # forces 
    forces = init_forces(fdsg3d)

    # progressbar 
    if fdsg3d.settings.config["settings"]["show_progress_in_console"] 
        p = Progress(fdsg3d.time.nt; showspeed=true)
    end 


    # time loop 
    for ti in 1:fdsg3d.time.nt 

        # save receiver
        @CUDA.allowscalar save_geophones!(fdsg3d.geophones, vx, vy, vz, ti)
        save_das!(fdsg3d.das, sxx, syy, szz, ti, fdsg3d.settings.float)
        save_snapshots!(fdsg3d.snapshots, vx, vy, vz, sxx, sxy, sxz, syy, syz, szz, ti, fdsg3d.settings.float)

        # forward
        CUDA.@cuda threads=THREADS_PER_BLOCK blocks=BLOCKS_PER_GRID update_normal_stresses_cuda(
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

        CUDA.@cuda threads=THREADS_PER_BLOCK blocks=BLOCKS_PER_GRID update_sxy_cuda(
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

        CUDA.@cuda threads=THREADS_PER_BLOCK blocks=BLOCKS_PER_GRID update_sxz_cuda(
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

        CUDA.@cuda threads=THREADS_PER_BLOCK blocks=BLOCKS_PER_GRID update_syz_cuda(
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

        CUDA.@cuda threads=THREADS_PER_BLOCK blocks=BLOCKS_PER_GRID update_vx_cuda(
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

        CUDA.@cuda threads=THREADS_PER_BLOCK blocks=BLOCKS_PER_GRID update_vy_cuda(
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

        CUDA.@cuda threads=THREADS_PER_BLOCK blocks=BLOCKS_PER_GRID update_vz_cuda(
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

        @CUDA.allowscalar apply_forces!(forces, 
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