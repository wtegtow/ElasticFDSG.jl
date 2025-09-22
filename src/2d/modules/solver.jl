"""
Leapfrog scheme:
- Time and field staggered grid. 
- Second order accurate in time. 
"""
function solve!(fdsg2d::FDSG2D; block_size=(32,32))

    kernel_params = kernel_args(fdsg2d);
    field_params = field_args(fdsg2d);
    backend = KernelAbstractions.get_backend(field_params.sxx);

    # init progressbar 
    if fdsg2d.settings.showinfo
        p = Progress(fdsg2d.time.nt; showspeed=true)
    end;

    # time loop 
    for ti in 1:fdsg2d.time.nt 
        # save receiver
        @unpack vx, vz, sxx, sxz, szz = field_params
        GPUArrays.@allowscalar save_geophones!(fdsg2d.geophones, vx, vz, ti)
        save_das!(fdsg2d.das, sxx, szz, ti)
        save_snapshots!(fdsg2d.snapshots, vx, vz, sxx, sxz, szz, ti)
        # forward stresses
        forward_normal_stresses!(kernel_params,field_params; backend=backend, block_size=block_size)
        forward_sxz!(kernel_params,field_params; backend=backend, block_size=block_size)
        KernelAbstractions.synchronize(backend)
        # apply forces 
        GPUArrays.@allowscalar apply_forces!(fdsg2d, field_params, ti)
        # forward velocities
        forward_vx!(kernel_params,field_params; backend=backend, block_size=block_size)
        forward_vz!(kernel_params,field_params; backend=backend, block_size=block_size)
        KernelAbstractions.synchronize(backend)
        # progress bar
        if fdsg2d.settings.showinfo next!(p) end;
    end
end