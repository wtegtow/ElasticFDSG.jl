# =============================================================================
# cpu solver
# =============================================================================

include(joinpath(@__DIR__,"update_fields.jl"))

function solver_cpu!(fdsg2d::FDSG2D)

    # update fields 
    function update_fields!(fdsg2d::FDSG2D)
        # update stresses
        update_normal_stresses!(fdsg2d.settings,fdsg2d.domain,fdsg2d.elastic,fdsg2d.fields,fdsg2d.pml,fdsg2d.time)
        update_sxy!(fdsg2d.settings,fdsg2d.domain,fdsg2d.elastic,fdsg2d.fields,fdsg2d.pml,fdsg2d.time)

        # update velocities
        update_vx!(fdsg2d.settings,fdsg2d.domain,fdsg2d.elastic,fdsg2d.fields,fdsg2d.pml,fdsg2d.time)   
        update_vy!(fdsg2d.settings,fdsg2d.domain,fdsg2d.elastic,fdsg2d.fields,fdsg2d.pml,fdsg2d.time)

    end;

    # forces 
    forces = get_forces(fdsg2d)

    # time loop
    if fdsg2d.settings.config["settings"]["show_progress_in_console"] 
        p = Progress(fdsg2d.time.nt; showspeed=true)
    end 
    
    for ti in 1:fdsg2d.time.nt 

        # save receiver 
        save_geophones!(fdsg2d.geophones, fdsg2d.fields.vx, fdsg2d.fields.vy, ti)
        save_das!(fdsg2d.das, fdsg2d.fields.sxx, fdsg2d.fields.syy, ti, fdsg2d.settings.float)
        save_snapshots!(fdsg2d.snapshots, fdsg2d.fields.vx, fdsg2d.fields.vy, fdsg2d.fields.sxx, fdsg2d.fields.sxy, fdsg2d.fields.syy, ti, fdsg2d.settings.float)
    
        # forward 
        update_fields!(fdsg2d)
        apply_forces!(forces, 
                      fdsg2d.fields.vx, fdsg2d.fields.vy, 
                      fdsg2d.source.yid, fdsg2d.source.xid, 
                      fdsg2d.domain.dx, fdsg2d.domain.dy, 
                      fdsg2d.source.STF, 
                      ti) 

        # update progress bar
        if fdsg2d.settings.config["settings"]["show_progress_in_console"] next!(p) end

        # perform garbage collection every 100 iterations
        if ti % 100 == 0
            GC.gc()
        end
    end
end