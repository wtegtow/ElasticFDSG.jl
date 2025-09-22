function iwindow(fdsg2d::FDSG2D)

    if fdsg2d.settings.showinfo
        
    # a rough memory approximation 
    grid_dim = prod(fdsg2d.domain.dim)
    num_fields = 5 
    num_lookup = 2
    num_pml_memory_variables = 8
    num_unique_ctensors = maximum(fdsg2d.elastic.c_lookup)
    pml_dim = length(fdsg2d.domain.pml_lookup)

    field_memory = grid_dim * num_fields * sizeof(fdsg2d.settings.float)
    lookup_memory = grid_dim * num_lookup * sizeof(fdsg2d.settings.float)
    pml_memory =  pml_dim * num_pml_memory_variables * sizeof(fdsg2d.settings.float)
    c_tensor_memory = (num_unique_ctensors * 10) * sizeof(fdsg2d.settings.float)

    memory_appr = (field_memory + lookup_memory + pml_memory + c_tensor_memory) / 1024^3
    memory_appr = round(memory_appr, digits=2)

    n_snap = 0
    if fdsg2d.snapshots.n > 0 
        n_snap = fdsg2d.snapshots.n * length(fdsg2d.snapshots.t) * length(fdsg2d.snapshots.fieldnames)
    end;

    # receiver 
    @printf("╔══════════════════════════════════════════════\n")
    @printf("║ FDSG 2D - Summary                   \n")
    @printf("╠══════════════════════════════════════════════\n")
    @printf("║ System: \n") 
    @printf("║   Device: %-1s                        \n", fdsg2d.settings.device_name)
    @printf("║   Approx. Memory: %-1.2f GB       \n", memory_appr)
    @printf("║ Grid: \n") 
    @printf("║   x: %-1s: %-1s: %-1s            \n", fdsg2d.domain.x0,round(fdsg2d.domain.dx,digits=2),round(fdsg2d.domain.xend,digits=2))
    @printf("║   z: %-1s: %-1s: %-1s            \n", fdsg2d.domain.z0,round(fdsg2d.domain.dz,digits=2),round(fdsg2d.domain.zend,digits=2))
    @printf("║   Number of Nodes: %-1d               \n", grid_dim)
    @printf("║ Time: \n") 
    @printf("║   Interval: %-1s: %-1s: %-1s              \n", fdsg2d.time.t0,round(fdsg2d.time.dt,digits=5),round(fdsg2d.time.tend, digits=5))
    @printf("║   Number of Timesteps: %-1d           \n", fdsg2d.time.nt)
    @printf("║ Source: \n") 
    @printf("║   x,z: %-1s, %-1s            \n", fdsg2d.source.x, fdsg2d.source.z)
    @printf("║   fdom, t0: %-1s, %-1s           \n", fdsg2d.settings.config["source"]["dominant_frequency"], fdsg2d.settings.config["source"]["wavelet_center"])
    @printf("║ Receiver: \n")
    @printf("║   Number of Geophones: %-1d            \n", fdsg2d.geophones.n)
    @printf("║   Number of Snapshots: %-1s       \n", fdsg2d.snapshots.n)
    @printf("║   Number of Fibers: %-1d, %-1d            \n", fdsg2d.das.x_aligned.n, fdsg2d.das.z_aligned.n)
    @printf("╚══════════════════════════════════════════════\n")
    end
end