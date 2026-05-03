function _print_summary(fdsg::FDSG)
        # unpack
        c  = fdsg.config
        d  = fdsg.domain
        e  = fdsg.elastic
        t  = fdsg.time
        s  = fdsg.source
        N  = c.dim
        fp = eval(Symbol(c.dict["settings"]["precision"]))

        #  memory estimate 
        n_grid    = prod(d.shape)
        n_pml     = count(d.pml_lookup .> 0)
        n_unique  = maximum(e.c_lookup)
        n_stiff   = N == 2 ? 5 : 10
        n_fld     = N == 2 ? 5 : 9
        n_pmemvar = N == 2 ? 8 : 18
        nb        = sizeof(fp)
        mem_gb    = (n_grid * (n_fld + 2) + n_pml * n_pmemvar + n_unique * n_stiff) * nb / 1024^3

        # CFL / Nyquist
        spacings  = map(c -> step(c), d.coordinates)
        dx_min    = minimum(spacings)
        courant   = t.dt * e.vmax * sqrt(N / dx_min^2)
        λ_min     = e.vmin / s.fdom
        ppw       = λ_min / dx_min

        # boundaries
        bnd       = c.dict["boundaries"]
        axis_names = N == 2 ? ["x", "z"] : ["x", "y", "z"]
        pml_sides = [k * side for k in axis_names for side in ["start","end"]
                     if get(bnd, k*side, "none") == "absorbing"]

        # grid ranges
        coords = d.coordinates

        w    = 80
        div  = "├" * "─"^w
        top  = "┌" * "─"^w
        bot  = "└" * "─"^w
        head = k -> "├── $k"
        row  = (k, v) -> @sprintf("│    %-26s %s\n", k, string(v))

        println(top)
        println("│  ElasticFDSG — $(N)D Simulation Summary")
        println(div)
        println(head("System"))
        print(row("Device",    fdsg.device.name))
        print(row("Precision", fp))
        print(row("Est. Memory", @sprintf("%.2f GB", mem_gb)))
        println(div)
        println(head("Grid"))
        for (i, ax) in enumerate(axis_names)
            # inner (model) range
            i0 = round(d.coordinates[i][d.inner_ids[i][begin]], digits=1)
            i1 = round(d.coordinates[i][d.inner_ids[i][end]],   digits=1)
            dx = round(step(coords[i]), digits=2)
            ni = length(d.inner_ids[i])
            # extended (with PML/ghost padding)
            e0 = round(coords[i][begin], digits=1)
            e1 = round(coords[i][end],   digits=1)
            ne = length(coords[i])
            print(row("$ax (inner)",    @sprintf("%.1f : %.2f : %.1f  (n=%d)", i0, dx, i1, ni)))
            print(row("$ax (extended)", @sprintf("%.1f : %.2f : %.1f  (n=%d)", e0, dx, e1, ne)))
        end
        print(row("Total nodes", prod(length.(coords))))
        print(row("PML nodes",   n_pml))
        print(row("PML sides",   isempty(pml_sides) ? "none" : join(pml_sides, ", ")))
        println(div)
        println(head("Time"))
        print(row("Interval", @sprintf("%.4f : %.2e : %.4f", t.t0, t.dt, t.tend)))
        print(row("Timesteps", t.nt))
        print(row("CFL number", @sprintf("%.3f %s", courant, courant > 1 ? "⚠ UNSTABLE" : "✓")))
        print(row("Points/wavelength", @sprintf("%.1f %s", ppw, ppw < 5 ? "⚠ LOW" : "✓")))
        println(div)
        println(head("Source"))
        if N == 2
            print(row("Position (x,z)", @sprintf("(%.1f, %.1f)", s.x, s.z)))
        else
            print(row("Position (x,y,z)", @sprintf("(%.1f, %.1f, %.1f)", s.x, s.y, s.z)))
        end
        print(row("fdom / t0", @sprintf("%.1f Hz / %.4f s", s.fdom, c.dict["source"]["wavelet_center"])))
        println(div)
        println(head("Receivers"))
        print(row("Geophones", fdsg.geophones.n))
        print(row("DAS fibers", sum(f.n for f in fdsg.das.fibers; init=0)))
        print(row("Snapshots", fdsg.snapshots.n))
        println(bot)
    end