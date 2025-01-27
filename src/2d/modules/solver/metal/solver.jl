# =============================================================================
# metal solver
# =============================================================================

include(joinpath(@__DIR__,"update_fields.jl"))

function solver_metal(fdsg2d::FDSG2D)
    
    # switch to GPU 
    N = fdsg2d.settings.N
    float = fdsg2d.settings.float
    c = MtlArray(fdsg2d.settings.c)

    yid_src = fdsg2d.source.yid 
    xid_src = fdsg2d.source.xid 
    STF = fdsg2d.source.STF

    nx = fdsg2d.domain.nx
    ny = fdsg2d.domain.ny
    dx = Float32(fdsg2d.domain.dx)
    dy = Float32(fdsg2d.domain.dy)
    dt = Float32(fdsg2d.time.dt)

    c11 = MtlArray(fdsg2d.elastic.c11)
    c13 = MtlArray(fdsg2d.elastic.c13)
    c33 = MtlArray(fdsg2d.elastic.c33)
    c55 = MtlArray(fdsg2d.elastic.c55)
    rho = MtlArray(fdsg2d.elastic.rho)

    vx = MtlArray(fdsg2d.fields.vx)
    vy = MtlArray(fdsg2d.fields.vy)
    sxx = MtlArray(fdsg2d.fields.sxx) 
    sxy = MtlArray(fdsg2d.fields.sxy)
    syy = MtlArray(fdsg2d.fields.syy)

    vx_x_old = MtlArray(fdsg2d.pml.vx_x_old)
    vy_y_old = MtlArray(fdsg2d.pml.vy_y_old)
    vy_x_old = MtlArray(fdsg2d.pml.vy_x_old)
    vx_y_old = MtlArray(fdsg2d.pml.vx_y_old)
    sxx_x_old = MtlArray(fdsg2d.pml.sxx_x_old)
    sxy_y_old = MtlArray(fdsg2d.pml.sxy_y_old)
    sxy_x_old = MtlArray(fdsg2d.pml.sxy_x_old)
    syy_y_old = MtlArray(fdsg2d.pml.syy_y_old)

    K_x_evn = MtlArray(fdsg2d.pml.K_x_evn)
    K_x_odd = MtlArray(fdsg2d.pml.K_x_odd)
    a_x_evn = MtlArray(fdsg2d.pml.a_x_evn)
    a_x_odd = MtlArray(fdsg2d.pml.a_x_odd)
    b_x_evn = MtlArray(fdsg2d.pml.b_x_evn)
    b_x_odd = MtlArray(fdsg2d.pml.b_x_odd)

    K_y_evn = MtlArray(fdsg2d.pml.K_y_evn)
    K_y_odd = MtlArray(fdsg2d.pml.K_y_odd)
    a_y_evn = MtlArray(fdsg2d.pml.a_y_evn)
    a_y_odd = MtlArray(fdsg2d.pml.a_y_odd)
    b_y_evn = MtlArray(fdsg2d.pml.b_y_evn)
    b_y_odd = MtlArray(fdsg2d.pml.b_y_odd)

    pml_hashmap =  MtlArray(Int32.(fdsg2d.domain.pml_points_hash))

    # compile kernel
    update_normal_stresses = @metal launch=false update_normal_kernel!(N,c,dt,ny,nx,dy,dx,c11,c13,c33,vx,vy,sxx,syy,a_x_odd,b_x_odd,K_x_odd,a_y_evn,b_y_evn,K_y_evn,vx_x_old,vy_y_old,pml_hashmap)
    update_sxy = @metal launch=false update_sxy_kernel!(N,c,dt,ny,nx,dy,dx,c55,vx,vy,sxy,a_x_evn,b_x_evn,K_x_evn,a_y_odd,b_y_odd,K_y_odd,vy_x_old,vx_y_old,pml_hashmap)
    update_vx = @metal launch=false update_vx_kernel!(N,c,dt,ny,nx,dy,dx,rho,sxx,sxy,vx,a_x_evn,b_x_evn,K_x_evn,a_y_evn,b_y_evn,K_y_evn,sxx_x_old,sxy_y_old,pml_hashmap)
    update_vy = @metal launch=false update_vy_kernel!(N,c,dt,ny,nx,dy,dx,rho,sxy,syy,vy,a_x_odd,b_x_odd,K_x_odd,a_y_odd,b_y_odd,K_y_odd,sxy_x_old,syy_y_old,pml_hashmap)

    dim_arg_sort = sort(collect(size(vx)),rev=true)
    w = min(size(dim_arg_sort, 1), update_normal_stresses.pipeline.threadExecutionWidth)
    h = min(size(dim_arg_sort, 2), update_normal_stresses.pipeline.threadExecutionWidth,
                                   update_normal_stresses.pipeline.maxTotalThreadsPerThreadgroup ÷ w)

    threads = (w, h)
    groups  = cld.(size(vx), threads)

    # forces 
    forces = get_forces(fdsg2d)

    # progress bar
    if fdsg2d.settings.config["settings"]["show_progress_in_console"]
        p = Progress(fdsg2d.time.nt; showspeed=true)
    end 

    # time loop
    for ti in 1:fdsg2d.time.nt 

        # save receiver
        @Metal.allowscalar save_geophones!(fdsg2d.geophones, vx, vy, ti)
        save_snapshots!(fdsg2d.snapshots, vx, vy, sxx, sxy, syy, ti, float)
        save_das!(fdsg2d.das, sxx, syy, ti, float)

        # forward 
        update_normal_stresses(N,c,dt,ny,nx,dy,dx,c11,c13,c33,vx,vy,sxx,syy,a_x_odd,b_x_odd,K_x_odd,a_y_evn,b_y_evn,K_y_evn,vx_x_old,vy_y_old,pml_hashmap, threads=threads, groups=groups)
        update_sxy(N,c,dt,ny,nx,dy,dx,c55,vx,vy,sxy,a_x_evn,b_x_evn,K_x_evn,a_y_odd,b_y_odd,K_y_odd,vy_x_old,vx_y_old,pml_hashmap, threads=threads, groups=groups)
        update_vx(N,c,dt,ny,nx,dy,dx,rho,sxx,sxy,vx,a_x_evn,b_x_evn,K_x_evn,a_y_evn,b_y_evn,K_y_evn,sxx_x_old,sxy_y_old,pml_hashmap, threads=threads, groups=groups)
        update_vy(N,c,dt,ny,nx,dy,dx,rho,sxy,syy,vy,a_x_odd,b_x_odd,K_x_odd,a_y_odd,b_y_odd,K_y_odd,sxy_x_old,syy_y_old,pml_hashmap, threads=threads, groups=groups)
        @Metal.allowscalar apply_forces!(forces, vx, vy, yid_src, xid_src, dx, dy, STF, ti)
        
        # update progress bar
        if fdsg2d.settings.config["settings"]["show_progress_in_console"] next!(p) end

         # perform garbage collection every 100 iterations
        if ti % 100 == 0
            GC.gc()
        end   
    end 
end

