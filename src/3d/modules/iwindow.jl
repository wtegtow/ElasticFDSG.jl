function iwindow(fdsg3d::FDSG3D)

    if fdsg3d.settings.showinfo
        
    # a rough memory approximation 
    grid_dim = prod(fdsg3d.domain.dim)
    num_fields = 9 
    num_lookup = 2
    num_pml_memory_variables = 18
    num_unique_ctensors = maximum(fdsg3d.elastic.c_lookup)
    pml_dim = length(fdsg3d.domain.pml_lookup)

    field_memory = grid_dim * num_fields * sizeof(fdsg3d.settings.float)
    lookup_memory = grid_dim * num_lookup * sizeof(fdsg3d.settings.float)
    pml_memory =  pml_dim * num_pml_memory_variables * sizeof(fdsg3d.settings.float)
    c_tensor_memory = (num_unique_ctensors * 10) * sizeof(fdsg3d.settings.float)

    memory_appr = (field_memory + lookup_memory + pml_memory + c_tensor_memory) / 1024^3
    memory_appr = round(memory_appr, digits=2)

    n_snap = 0
    if fdsg3d.snapshots.n > 0 
        n_snap = fdsg3d.snapshots.n * length(fdsg3d.snapshots.t) * length(fdsg3d.snapshots.fieldnames)
    end;

    # receiver 
    @printf("╔══════════════════════════════════════════════\n")
    @printf("║ FDSG 3D - Summary                   \n")
    @printf("╠══════════════════════════════════════════════\n")
    @printf("║ System: \n") 
    @printf("║   Device: %-1s                        \n", fdsg3d.settings.device_name)
    @printf("║   Approx. Memory: %-1.2f GB       \n", memory_appr)
    @printf("║ Grid: \n") 
    @printf("║   x: %-1s: %-1s: %-1s            \n", fdsg3d.domain.x0,round(fdsg3d.domain.dx,digits=2),round(fdsg3d.domain.xend,digits=2))
    @printf("║   y: %-1s: %-1s: %-1s            \n", fdsg3d.domain.y0,round(fdsg3d.domain.dy,digits=2),round(fdsg3d.domain.yend,digits=2))
    @printf("║   z: %-1s: %-1s: %-1s            \n", fdsg3d.domain.z0,round(fdsg3d.domain.dz,digits=2),round(fdsg3d.domain.zend,digits=2))
    @printf("║   Number of Nodes: %-1d               \n", grid_dim)
    @printf("║ Time: \n") 
    @printf("║   Interval: %-1s: %-1s: %-1s              \n", fdsg3d.time.t0,round(fdsg3d.time.dt,digits=5),round(fdsg3d.time.tend, digits=5))
    @printf("║   Number of Timesteps: %-1d           \n", fdsg3d.time.nt)
    @printf("║ Source: \n") 
    @printf("║   x,y,z: %-1s, %-1s, %-1s            \n", fdsg3d.source.x, fdsg3d.source.y, fdsg3d.source.z)
    @printf("║   fdom, t0: %-1s, %-1s           \n", fdsg3d.settings.config["source"]["dominant_frequency"], fdsg3d.settings.config["source"]["wavelet_center"])
    @printf("║ Receiver: \n")
    @printf("║   Number of Geophones: %-1d            \n", fdsg3d.geophones.n)
    @printf("║   Number of Snapshots: %-1s       \n", n_snap)
    @printf("║   Number of Fibers: %-1d, %-1d, %-1d            \n", fdsg3d.das.x_aligned.n, fdsg3d.das.y_aligned.n, fdsg3d.das.z_aligned.n)
    @printf("╚══════════════════════════════════════════════\n")
    end
end