# ------ some checks ------

function check_dir(dirpath)
    if ispath(dirpath) != true
        error("Output directory $dirpath not found")
    end
end;

function check_format(format)
    supported_formats = ["h5"] 
    if format ∉ supported_formats
        error("Output format $format not supported. Supported formats: $supported_formats")
    end
end

function check_filename(filename, Log)
    if isfile(filename)
        add_message!(Log,"Warning: Output filename: $filename already exists. Once the calculations are complete, this file will be overwritten.")
    end
end;

# --------- save ---------

function save_results_h5(fdsg3d::FDSG3D)
  
    file = h5open(fdsg3d.settings.filename, "w")
    # coordinates 

    write(file, "x_coordinates", collect(fdsg3d.domain.xcoords[fdsg3d.domain.inner_id[3]]))
    write(file, "y_coordinates", collect(fdsg3d.domain.ycoords[fdsg3d.domain.inner_id[2]]))
    write(file, "z_coordinates", collect(fdsg3d.domain.zcoords[fdsg3d.domain.inner_id[1]]))
    # time 
    write(file, "time", collect(fdsg3d.time.t))
    # geophones 
    write(file, "geophones_data", fdsg3d.geophones.data)
    write(file, "geophones_coords", Array(fdsg3d.geophones.coords))
    write(file, "geophones_x_indices", fdsg3d.geophones.xids)
    write(file, "geophones_y_indices", fdsg3d.geophones.yids)
    write(file, "geophones_z_indices", fdsg3d.geophones.zids)

    # das x_aligned
    for i in 1:fdsg3d.das.x_aligned.n
        write(file, "das_x_aligned_data_$i", fdsg3d.das.x_aligned.data[i])

        xids, yids, zids = fdsg3d.das.x_aligned.ids[i]
        ids = zeros(3,size(xids)[1])
        ids[1,:] .= xids
        ids[2,:] .= yids
        ids[3,:] .= zids
        write(file, "das_x_aligned_indices_$i", ids)

        xcoords, ycoords, zcoords = fdsg3d.das.x_aligned.coords[i]
        coords = zeros(3,size(xcoords)[1])
        coords[1,:] .= xcoords
        coords[2,:] .= ycoords
        coords[3,:] .= zcoords
        write(file, "das_x_aligned_coords_$i", coords)
    end

    # das y_aligned
    for i in 1:fdsg3d.das.y_aligned.n
        write(file, "das_y_aligned_data_$i", fdsg3d.das.y_aligned.data[i])

        xids, yids, zids = fdsg3d.das.y_aligned.ids[i]
        ids = zeros(3,size(yids)[1])
        ids[1,:] .= xids
        ids[2,:] .= yids
        ids[3,:] .= zids
        write(file, "das_y_aligned_indices_$i", ids)

        xcoords, ycoords, zcoords = fdsg3d.das.y_aligned.coords[i]
        coords = zeros(3,size(ycoords)[1])
        coords[1,:] .= xcoords
        coords[2,:] .= ycoords
        coords[3,:] .= zcoords
        write(file, "das_y_aligned_coords_$i", coords)
    end

    # das z_aligned
    for i in 1:fdsg3d.das.z_aligned.n
        write(file, "das_z_aligned_data_$i", fdsg3d.das.z_aligned.data[i])

        xids, yids, zids = fdsg3d.das.z_aligned.ids[i]
        ids = zeros(3,size(zids)[1])
        ids[1,:] .= xids
        ids[2,:] .= yids
        ids[3,:] .= zids
        write(file, "das_z_aligned_indices_$i", ids)

        xcoords, ycoords, zcoords = fdsg3d.das.z_aligned.coords[i]
        coords = zeros(3,size(zcoords)[1])
        coords[1,:] .= xcoords
        coords[2,:] .= ycoords
        coords[3,:] .= zcoords
        write(file, "das_z_aligned_coords_$i", coords)
    end

    # snapshots (shrinked to physical domain)
    write(file, "snapshots_XYdata", fdsg3d.snapshots.XY[:,:,fdsg3d.domain.inner_id[2], fdsg3d.domain.inner_id[3]])
    write(file, "snapshots_XZdata", fdsg3d.snapshots.XZ[:,:,fdsg3d.domain.inner_id[1], fdsg3d.domain.inner_id[3]])
    write(file, "snapshots_YZdata", fdsg3d.snapshots.YZ[:,:,fdsg3d.domain.inner_id[1], fdsg3d.domain.inner_id[2]])

    write(file, "snapshots_fieldnames",fdsg3d.snapshots.fieldnames)
    write(file, "snapshots_times", fdsg3d.snapshots.t)
    write(file, "snapshots_times_indices",fdsg3d.snapshots.tid)

    close(file)

end


function save_results(fdsg3d::FDSG3D)
    if fdsg3d.settings.save == true
        # h5 format
        if occursin(".h5", fdsg3d.settings.filename)
            save_results_h5(fdsg3d)
            println("Results saved.")
        end
    end
end;


# ------- template --------

"""
Elastic_FDSG.dim3.configtemplate(path)

Creates and empty template configuration file for a 3D ElasticFDSG simulation. 
The user must fill the template afterwards. 

# Arguments
- `path::String`: Path where the template is saved.

# Returns
- `Nothing`: The function saves the template as specified in path.
"""
function configtemplate(path::String)

    template = """
    # This is a template configuration.yaml file for a 3D ElasticFDSG simulation.
    # Any other configuration file can be prepared in the same manner.
    # All keywords are case sensitive, and must be named as shown below.
    # Users must fill the whole template before running a simulation.  
    # The velocity model is prepared in another file.

    settings:
        device:                          # cpu / gpu-cuda / gpu-metal 
        precision:                       # Float64 / Float32
        spatial_derivative_order:        # (1-10)
        show_summary_in_console:         # true / false  
        show_progress_in_console:        # true / false
        save_results:                    # true / false  
        output:
            destination_folder: 
            file_name:  
            format: h5                   # h5 (only supported format yet) 

    time:
        start:                # [sec]
        end:                  # [sec]
        timestep:             # [sec]

    source:
        dominant_frequency:          # [Hz]
        wavelet_type:                # ricker / gauss1d 
        wavelet_center:              # (should be ≥ 1/fdom) [sec]
        amplitude:                        
        point_source:          
            use:                     # true / false  (if true, double_couple.use should be false)
            phi:                     # [°] Elevation
            theta:                   # [°] Azimuth 
        double_couple: 
            use:                     # true / false  (if true, point_source.use should be false)
            strike:                  # [°] ∈ [0, 360]
            dip:                     # [°] ∈ [0, 90]
            rake:                    # [°] ∈ [0, 360]
            anisotropic_moment_tensor:      # true / false
        location:
            x:              # [m]
            y:              # [m]
            z:              # [m]

    boundaries: # absorbing / else
        xstart:       
        xend:        
        ystart:        
        yend:          
        zstart:         
        zend:           

    pml:
        nlayer:                      # 5-15 recommended
        reflection_coefficient:      # 1e-10 recommended

    receivers:
        geophones: 
            #- { x: 100, y: 200, z: 300 }    #[m] Example geophone

        das:
            x_aligned: 
            #- { x_range: "0:2.5:500", y: 100, z: 100 }  #[m] Example das 
           
            y_aligned: 
            #- { x: 200, y_range: "100:2.5:400", z: 100 } #[m] Example das 

            z_aligned:
            #- { x: 300, y: 250, z_range: "0:2.5:400" } #[m] Example das 
                    
                    
        snapshots: # XY-, XZ-, YZ-planes only
            times: 
                # [0.25, 0.5, 0.75, 1.0] # [sec] Example
            fields: 
                # [vx, vy, vz, sxx, sxy, sxz, syy, syz, szz] # Example 
            origins:
                # { x: 250, y: 250, z: 250 } # [m] Coordinates of snapshot origin. Example: XY-plane snapshots will be at z=250m


    """


    open(path, "w") do io
        write(io, template)
    end
    println("Template saved at $path")
end