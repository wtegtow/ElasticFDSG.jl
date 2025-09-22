"""
Leapfrog scheme:
- Time and field staggered grid. 
- Second order accurate in time. 
"""
function solve!(fdsg3d::FDSG3D; block_size=(7,7,7))

    kernel_params = kernel_args(fdsg3d);
    field_params = field_args(fdsg3d);
    backend = KernelAbstractions.get_backend(field_params.sxx);

    # init progressbar 
    if fdsg3d.settings.showinfo
        p = Progress(fdsg3d.time.nt; showspeed=true)
    end;

    # time loop 
    for ti in 1:fdsg3d.time.nt 
        # save receiver
        @unpack vx, vy, vz, sxx, sxy, sxz, syy, syz, szz = field_params
        GPUArrays.@allowscalar save_geophones!(fdsg3d.geophones,vx, vy, vz, ti)
        save_das!(fdsg3d.das, sxx, syy, szz, ti)
        save_snapshots!(fdsg3d.snapshots, vx, vy, vz, sxx, sxy, sxz, syy, syz, szz, ti)
        # forward stresses
        forward_normal_stresses!(kernel_params,field_params; backend=backend, block_size=block_size)
        forward_sxy!(kernel_params,field_params; backend=backend, block_size=block_size)
        forward_sxz!(kernel_params,field_params; backend=backend, block_size=block_size)
        forward_syz!(kernel_params,field_params; backend=backend, block_size=block_size)
        KernelAbstractions.synchronize(backend)
        # apply forces 
        GPUArrays.@allowscalar apply_forces!(fdsg3d, field_params, ti)
        # forward velocities
        forward_vx!(kernel_params,field_params; backend=backend, block_size=block_size)
        forward_vy!(kernel_params,field_params; backend=backend, block_size=block_size)
        forward_vz!(kernel_params,field_params; backend=backend, block_size=block_size)
        KernelAbstractions.synchronize(backend)
        # progress bar
        if fdsg3d.settings.showinfo next!(p) end;
    end
end