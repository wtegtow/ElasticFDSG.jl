# =============================================================================
# cuda solver
# =============================================================================

include(joinpath(@__DIR__,"update_fields.jl"))

function solver_cuda(fdsg2d::FDSG2D)
    
    # switch to GPU 
    N = fdsg2d.settings.N
    FLOAT = fdsg2d.settings.float
    c = CuArray(fdsg2d.settings.c)

    yid_src = fdsg2d.source.yid 
    xid_src = fdsg2d.source.xid 
    STF = fdsg2d.source.STF

    nx = fdsg2d.domain.nx
    ny = fdsg2d.domain.ny
    dx = Float32(fdsg2d.domain.dx)
    dy = Float32(fdsg2d.domain.dy)
    dt = Float32(fdsg2d.time.dt)

    c11 = CuArray(fdsg2d.elastic.c11)
    c13 = CuArray(fdsg2d.elastic.c13)
    c33 = CuArray(fdsg2d.elastic.c33)
    c55 = CuArray(fdsg2d.elastic.c55)
    rho = CuArray(fdsg2d.elastic.rho)

    vx = CuArray(fdsg2d.fields.vx)
    vy = CuArray(fdsg2d.fields.vy)
    sxx = CuArray(fdsg2d.fields.sxx) 
    sxy = CuArray(fdsg2d.fields.sxy)
    syy = CuArray(fdsg2d.fields.syy)

    vx_x_old = CuArray(fdsg2d.pml.vx_x_old)
    vy_y_old = CuArray(fdsg2d.pml.vy_y_old)
    vy_x_old = CuArray(fdsg2d.pml.vy_x_old)
    vx_y_old = CuArray(fdsg2d.pml.vx_y_old)
    sxx_x_old = CuArray(fdsg2d.pml.sxx_x_old)
    sxy_y_old = CuArray(fdsg2d.pml.sxy_y_old)
    sxy_x_old = CuArray(fdsg2d.pml.sxy_x_old)
    syy_y_old = CuArray(fdsg2d.pml.syy_y_old)

    K_x_evn = CuArray(fdsg2d.pml.K_x_evn)
    K_x_odd = CuArray(fdsg2d.pml.K_x_odd)
    a_x_evn = CuArray(fdsg2d.pml.a_x_evn)
    a_x_odd = CuArray(fdsg2d.pml.a_x_odd)
    b_x_evn = CuArray(fdsg2d.pml.b_x_evn)
    b_x_odd = CuArray(fdsg2d.pml.b_x_odd)

    K_y_evn = CuArray(fdsg2d.pml.K_y_evn)
    K_y_odd = CuArray(fdsg2d.pml.K_y_odd)
    a_y_evn = CuArray(fdsg2d.pml.a_y_evn)
    a_y_odd = CuArray(fdsg2d.pml.a_y_odd)
    b_y_evn = CuArray(fdsg2d.pml.b_y_evn)
    b_y_odd = CuArray(fdsg2d.pml.b_y_odd)

    pml_hashmap =  CuArray(Int32.(fdsg2d.domain.pml_points_hash))

    THREADS_PER_BLOCK = (16, 16)
    BLOCKS_PER_GRID = (ceil(Int, ny / THREADS_PER_BLOCK[1]),
                       ceil(Int, nx / THREADS_PER_BLOCK[2]))

    # forces 
    forces = get_forces(fdsg2d)

    # progress bar
    if fdsg2d.settings.config["settings"]["show_progress_in_console"]
        p = Progress(fdsg2d.time.nt; showspeed=true)
    end 

    # time loop
    for ti in 1:fdsg2d.time.nt 

        # save receiver
        @CUDA.allowscalar save_geophones!(fdsg2d.geophones, vx, vy, ti)
        save_snapshots!(fdsg2d.snapshots, vx, vy, sxx, sxy, syy, ti, FLOAT)
        save_das!(fdsg2d.das, sxx, syy, ti, FLOAT)

        # forward
        CUDA.@cuda threads=THREADS_PER_BLOCK blocks=BLOCKS_PER_GRID update_normal_stresses_cuda(
            N, c, dt,
            ny, nx, dy, dx, 
            c11,c13,c33,
            vx,vy,
            sxx,syy,
            a_x_odd, b_x_odd, K_x_odd, 
            a_y_evn, b_y_evn, K_y_evn, 
            vx_x_old, vy_y_old,
            pml_hashmap)

        CUDA.@cuda threads=THREADS_PER_BLOCK blocks=BLOCKS_PER_GRID update_sxy_cuda(
            N, c, dt,
            ny, nx, dy, dx, 
            c55,
            vx,vy,
            sxy,
            a_x_evn, b_x_evn, K_x_evn, 
            a_y_odd, b_y_odd, K_y_odd, 
            vy_x_old,vx_y_old,
            pml_hashmap)

        CUDA.@cuda threads=THREADS_PER_BLOCK blocks=BLOCKS_PER_GRID update_vx_cuda(N, c, dt,
            ny, nx, dy, dx, 
            rho,
            sxx,sxy,
            vx,
            a_x_evn, b_x_evn, K_x_evn, 
            a_y_evn, b_y_evn, K_y_evn, 
            sxx_x_old,sxy_y_old,
            pml_hashmap)

        CUDA.@cuda threads=THREADS_PER_BLOCK blocks=BLOCKS_PER_GRID update_vy_cuda(N, c, dt,
            ny, nx, dy, dx, 
            rho,
            sxy,syy,
            vy,
            a_x_odd, b_x_odd, K_x_odd, 
            a_y_odd, b_y_odd, K_y_odd, 
            sxy_x_old,syy_y_old,
            pml_hashmap)
            
        @CUDA.allowscalar apply_forces!(forces, vx, vy, yid_src, xid_src, dx, dy, STF, ti)
        
        # update progress bar
        if fdsg2d.settings.config["settings"]["show_progress_in_console"] next!(p) end

        # perform garbage collection every 100 iterations
        if ti % 100 == 0
            GC.gc()
        end   
    end 
 end

    