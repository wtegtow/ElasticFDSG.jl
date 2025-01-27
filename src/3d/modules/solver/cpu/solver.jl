# =============================================================================
# cpu solver
# =============================================================================

include(joinpath(@__DIR__,"update_fields.jl"))

function solver_cpu!(fdsg3d::FDSG3D)

    # update fields 
    function update_fields!(fdsg3d::FDSG3D)

        # update stresses
        update_normal_stresses!(
            fdsg3d.fields.vx, fdsg3d.fields.vy, fdsg3d.fields.vz,
            fdsg3d.fields.sxx, fdsg3d.fields.syy, fdsg3d.fields.szz,
            fdsg3d.elastic.c11, fdsg3d.elastic.c12, fdsg3d.elastic.c13,
            fdsg3d.elastic.c22, fdsg3d.elastic.c23, fdsg3d.elastic.c33,
            fdsg3d.pml.vx_x_old, fdsg3d.pml.vy_y_old, fdsg3d.pml.vz_z_old,
            fdsg3d.pml.b_x_odd, fdsg3d.pml.a_x_odd, fdsg3d.pml.K_x_odd,
            fdsg3d.pml.b_y_evn, fdsg3d.pml.a_y_evn, fdsg3d.pml.K_y_evn,
            fdsg3d.pml.b_z_evn, fdsg3d.pml.a_z_evn, fdsg3d.pml.K_z_evn,
            fdsg3d.domain.nz, fdsg3d.domain.ny, fdsg3d.domain.nx,
            fdsg3d.domain.dx, fdsg3d.domain.dy, fdsg3d.domain.dz,
            fdsg3d.time.dt, fdsg3d.settings.c, fdsg3d.settings.N,
            fdsg3d.domain.pml_points_hash
        )

        update_sxy!(
            fdsg3d.fields.vx, fdsg3d.fields.vy,
            fdsg3d.fields.sxy,
            fdsg3d.elastic.c66,
            fdsg3d.pml.vy_x_old, fdsg3d.pml.vx_y_old,
            fdsg3d.pml.b_x_evn, fdsg3d.pml.a_x_evn, fdsg3d.pml.K_x_evn,
            fdsg3d.pml.b_y_odd, fdsg3d.pml.a_y_odd, fdsg3d.pml.K_y_odd,
            fdsg3d.domain.nz, fdsg3d.domain.ny, fdsg3d.domain.nx,
            fdsg3d.domain.dx, fdsg3d.domain.dy,
            fdsg3d.time.dt, fdsg3d.settings.c, fdsg3d.settings.N,
            fdsg3d.domain.pml_points_hash
        )

        update_sxz!(
            fdsg3d.fields.vx, fdsg3d.fields.vz,
            fdsg3d.fields.sxz,
            fdsg3d.elastic.c55,
            fdsg3d.pml.vz_x_old, fdsg3d.pml.vx_z_old,
            fdsg3d.pml.b_x_evn, fdsg3d.pml.a_x_evn, fdsg3d.pml.K_x_evn,
            fdsg3d.pml.b_z_odd, fdsg3d.pml.a_z_odd, fdsg3d.pml.K_z_odd,
            fdsg3d.domain.nz, fdsg3d.domain.ny, fdsg3d.domain.nx,
            fdsg3d.domain.dx, fdsg3d.domain.dz,
            fdsg3d.time.dt, fdsg3d.settings.c, fdsg3d.settings.N,
            fdsg3d.domain.pml_points_hash
        )

        update_syz!(
            fdsg3d.fields.vy, fdsg3d.fields.vz,
            fdsg3d.fields.syz,
            fdsg3d.elastic.c44,
            fdsg3d.pml.vz_y_old, fdsg3d.pml.vy_z_old,
            fdsg3d.pml.b_y_odd, fdsg3d.pml.a_y_odd, fdsg3d.pml.K_y_odd,
            fdsg3d.pml.b_z_odd, fdsg3d.pml.a_z_odd, fdsg3d.pml.K_z_odd,
            fdsg3d.domain.nz, fdsg3d.domain.ny, fdsg3d.domain.nx,
            fdsg3d.domain.dy, fdsg3d.domain.dz,
            fdsg3d.time.dt, fdsg3d.settings.c, fdsg3d.settings.N,
            fdsg3d.domain.pml_points_hash
        )

        # update velocities
        update_vx!(
            fdsg3d.fields.vx,
            fdsg3d.fields.sxx, fdsg3d.fields.sxy, fdsg3d.fields.sxz,
            fdsg3d.elastic.rho,
            fdsg3d.pml.sxx_x_old, fdsg3d.pml.sxy_y_old, fdsg3d.pml.sxz_z_old,
            fdsg3d.pml.b_x_evn, fdsg3d.pml.a_x_evn, fdsg3d.pml.K_x_evn,
            fdsg3d.pml.b_y_evn, fdsg3d.pml.a_y_evn, fdsg3d.pml.K_y_evn,
            fdsg3d.pml.b_z_evn, fdsg3d.pml.a_z_evn, fdsg3d.pml.K_z_evn,
            fdsg3d.domain.nz, fdsg3d.domain.ny, fdsg3d.domain.nx,
            fdsg3d.domain.dx, fdsg3d.domain.dy, fdsg3d.domain.dz,
            fdsg3d.time.dt, fdsg3d.settings.c, fdsg3d.settings.N,
            fdsg3d.domain.pml_points_hash
        )

        update_vy!(
            fdsg3d.fields.vy,
            fdsg3d.fields.sxy, fdsg3d.fields.syy, fdsg3d.fields.syz,
            fdsg3d.elastic.rho,
            fdsg3d.pml.sxy_x_old, fdsg3d.pml.syy_y_old, fdsg3d.pml.syz_z_old,
            fdsg3d.pml.b_x_odd, fdsg3d.pml.a_x_odd, fdsg3d.pml.K_x_odd,
            fdsg3d.pml.b_y_odd, fdsg3d.pml.a_y_odd, fdsg3d.pml.K_y_odd,
            fdsg3d.pml.b_z_evn, fdsg3d.pml.a_z_evn, fdsg3d.pml.K_z_evn,
            fdsg3d.domain.nz, fdsg3d.domain.ny, fdsg3d.domain.nx,
            fdsg3d.domain.dx, fdsg3d.domain.dy, fdsg3d.domain.dz,
            fdsg3d.time.dt, fdsg3d.settings.c, fdsg3d.settings.N,
            fdsg3d.domain.pml_points_hash
        )

        update_vz!(
            fdsg3d.fields.vz,
            fdsg3d.fields.sxz, fdsg3d.fields.syz, fdsg3d.fields.szz,
            fdsg3d.elastic.rho,
            fdsg3d.pml.sxz_x_old, fdsg3d.pml.syz_y_old, fdsg3d.pml.szz_z_old,
            fdsg3d.pml.b_x_odd, fdsg3d.pml.a_x_odd, fdsg3d.pml.K_x_odd,
            fdsg3d.pml.b_y_evn, fdsg3d.pml.a_y_evn, fdsg3d.pml.K_y_evn,
            fdsg3d.pml.b_z_odd, fdsg3d.pml.a_z_odd, fdsg3d.pml.K_z_odd,
            fdsg3d.domain.nz, fdsg3d.domain.ny, fdsg3d.domain.nx,
            fdsg3d.domain.dx, fdsg3d.domain.dy, fdsg3d.domain.dz,
            fdsg3d.time.dt, fdsg3d.settings.c, fdsg3d.settings.N,
            fdsg3d.domain.pml_points_hash
        )

    end;

    # forces
    forces = init_forces(fdsg3d)

    # progressbar 
    if fdsg3d.settings.config["settings"]["show_progress_in_console"] 
        p = Progress(fdsg3d.time.nt; showspeed=true)
    end 

    # time loop 
    for ti in 1:fdsg3d.time.nt 

        # save receiver
        save_geophones!(fdsg3d.geophones, 
                        fdsg3d.fields.vx, fdsg3d.fields.vy, fdsg3d.fields.vz, 
                        ti)

        save_das!(fdsg3d.das, 
                fdsg3d.fields.sxx, fdsg3d.fields.syy, fdsg3d.fields.szz, 
                ti, fdsg3d.settings.float)

        save_snapshots!(fdsg3d.snapshots, 
                        fdsg3d.fields.vx, fdsg3d.fields.vy, fdsg3d.fields.vz, 
                        fdsg3d.fields.sxx, fdsg3d.fields.sxy, fdsg3d.fields.sxz, fdsg3d.fields.syy, fdsg3d.fields.syz, fdsg3d.fields.szz, 
                        ti, fdsg3d.settings.float)

        # forward
        update_fields!(fdsg3d)
        apply_forces!(forces, 
                    fdsg3d.fields.vx, fdsg3d.fields.vy, fdsg3d.fields.vz,
                    fdsg3d.source.xid, fdsg3d.source.yid, fdsg3d.source.zid,
                    fdsg3d.domain.dx, fdsg3d.domain.dy, fdsg3d.domain.dz,
                    fdsg3d.source.STF, fdsg3d.source.rhosrc, fdsg3d.time.dt, 
                    ti)

        # update progress bar
        if fdsg3d.settings.config["settings"]["show_progress_in_console"] next!(p) end

        # garbage collection every 100 iterations
        if ti % 100 == 0
            GC.gc()
        end
    end
end;

