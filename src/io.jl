const _DAS_AXIS_NAMES_2D = ("x_aligned", "z_aligned")
const _DAS_AXIS_NAMES_3D = ("x_aligned", "y_aligned", "z_aligned")

function save_results(fdsg::FDSG)

    path = get(fdsg.config.dict["settings"], "output_file", nothing)
    (isnothing(path) || path == "null") && return

    d  = fdsg.domain
    t  = fdsg.time
    s  = fdsg.source
    N  = fdsg.config.dim

    # inner-axis coordinate ranges (StepRangeLen → Vector for HDF5)
    inner_coords = [collect(d.coordinates[k][d.inner_ids[k]]) for k in 1:N]

    HDF5.h5open(path, "w") do file

        # grid
        g1 = HDF5.create_group(file, "grid")
        ax = N == 2 ? ("x", "z") : ("x", "y", "z")
        for (k, name) in enumerate(ax)
            g1["$(name)_coordinates"] = inner_coords[k]
        end

        # time
        g2 = HDF5.create_group(file, "time")
        g2["t0"]   = t.t0
        g2["tend"] = t.tend
        g2["dt"]   = t.dt
        g2["time"] = collect(t.t)

        # source
        g3 = HDF5.create_group(file, "source")
        g3["dominant_frequency"]    = s.fdom
        g3["source_time_function"]  = collect(s.stf)
        g3["stf_d1"]                = collect(s.stf_d1)
        g3["Mxx"] = s.Mxx
        g3["Mxz"] = s.Mxz
        g3["Mzz"] = s.Mzz
        if N == 2
            g3["location"] = [s.x, s.z]
        else
            g3["Mxy"] = s.Mxy
            g3["Myy"] = s.Myy
            g3["Myz"] = s.Myz
            g3["location"] = [s.x, s.y, s.z]
        end

        # geophones 
        geo = fdsg.geophones
        if geo.n > 0
            g4 = HDF5.create_group(file, "geophones")
            for i in 1:geo.n
                gi = HDF5.create_group(g4, "geophone_$i")
                gi["data"]     = geo.data[i, :, :]   # (ncomp, nt)
                gi["location"] = geo.coords[:, i]     # (ndim,)
            end
        end

        # DAS
        axis_names = N == 2 ? _DAS_AXIS_NAMES_2D : _DAS_AXIS_NAMES_3D
        das_has_data = any(fg -> fg.n > 0, fdsg.das.fibers)
        if das_has_data
            g5 = HDF5.create_group(file, "das")
            for (fg, aname) in zip(fdsg.das.fibers, axis_names)
                fg.n == 0 && continue
                ga = HDF5.create_group(g5, aname)
                for i in 1:fg.n
                    gfi = HDF5.create_group(ga, "fiber_$i")
                    gfi["data"]     = fg.data[i]       # (nch, nt)
                    gfi["location"] = fg.coords[i]     # (ndim, nch)
                end
            end
        end

        # snapshots
        snap = fdsg.snapshots
        if snap.n > 0
            g6 = HDF5.create_group(file, "snapshots")
            g6["fields"] = snap.fieldnames
            if N == 2
                # data: (ntime, nfields, nx, nz) — crop to inner indices
                ix = d.inner_ids[1]
                iz = d.inner_ids[2]
                g6["XZ"] = snap.data[:, :, ix, iz]
                g6["times"] = [collect(t.t)[ti] for ti in sort(collect(keys(snap.tid_map)))]
            else
                ix = d.inner_ids[1]
                iy = d.inner_ids[2]
                iz = d.inner_ids[3]
                g6["times"] = [collect(t.t)[ti] for ti in sort(collect(keys(snap.tid_map)))]
                for (p, gids) in enumerate(snap.grid_ids)
                    gp = HDF5.create_group(g6, "plane_$p")
                    gp["XY_data"]   = snap.XY[p, :, :, ix, iy]   # (ntime, nfields, nx_inner, ny_inner)
                    gp["XZ_data"]   = snap.XZ[p, :, :, ix, iz]
                    gp["YZ_data"]   = snap.YZ[p, :, :, iy, iz]
                    gp["XY_plane_Z"] = d.coordinates[3][gids[3]]
                    gp["XZ_plane_Y"] = d.coordinates[2][gids[2]]
                    gp["YZ_plane_X"] = d.coordinates[1][gids[1]]
                end
            end
        end
    end

    if fdsg.config.dict["settings"]["verbose"]
         println("Results saved.")
    end
    return path
end

"""
    load_results(path) -> Dict{String, Any}

Load simulation results from an HDF5 file into a nested Julia dictionary.

# Example
```julia
data = load_results("output.h5")
vx_snap = data["snapshots"]["XZ"]   # 2D wavefield snapshots
geo1    = data["geophones"]["geophone_1"]["data"]
```
"""
function load_results(path::String)
    @assert isfile(path) "Results file not found: $path"
    HDF5.h5open(path, "r") do file
        return _h5_to_dict(file)
    end
end

function _h5_to_dict(obj)
    if isa(obj, HDF5.File) || isa(obj, HDF5.Group)
        return Dict{String, Any}(k => _h5_to_dict(obj[k]) for k in keys(obj))
    elseif isa(obj, HDF5.Dataset)
        return HDF5.read(obj)
    else
        return obj
    end
end