# ---------- checks  ----------

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

# ---------- save ----------

function save_results(fdsg2d::FDSG2D)
    if fdsg2d.settings.save == true
        
        # h5 format
        if occursin(".h5", fdsg2d.settings.filename)
            save_results_h5(fdsg2d)
            println("Results saved.")
        end
    end
end;


function save_results_h5(fdsg2d::FDSG2D)
  
    file = h5open(fdsg2d.settings.filename, "w")
    # coordinates 
    write(file, "x_coordinates", collect(fdsg2d.domain.xcoords[fdsg2d.domain.inner_id[2]]))
    write(file, "y_coordinates", collect(fdsg2d.domain.ycoords[fdsg2d.domain.inner_id[1]]))
    # time 
    write(file, "time", collect(fdsg2d.time.t))
    # geophones 
    write(file, "geophones_data", fdsg2d.geophones.data)
    write(file, "geophones_coords", Array(fdsg2d.geophones.coords))
    write(file, "geophones_x_indices", fdsg2d.geophones.xids)
    write(file, "geophones_y_indices", fdsg2d.geophones.yids)
    # das x_aligned
    for i in 1:fdsg2d.das.x_aligned.n
        write(file, "das_x_aligned_data_$i", fdsg2d.das.x_aligned.data[i])
        xids, yids = fdsg2d.das.x_aligned.ids[i]
        ids = zeros(2,size(xids)[1])
        ids[1,:] .= xids
        ids[2,:] .= yids
        write(file, "das_x_aligned_indices_$i", ids)

        xcoords, ycoords = fdsg2d.das.x_aligned.coords[i]
        coords = zeros(2,size(xcoords)[1])
        coords[1,:] .= xcoords
        coords[2,:] .= ycoords
        write(file, "das_x_aligned_coords_$i", coords)
    end

    # das y_aligned
    for i in 1:fdsg2d.das.y_aligned.n
        write(file, "das_y_aligned_data_$i", fdsg2d.das.y_aligned.data[i])

        xids, yids = fdsg2d.das.y_aligned.ids[i]
        ids = zeros(2,size(yids)[1])
        ids[1,:] .= xids
        ids[2,:] .= yids
        write(file, "das_y_aligned_indices_$i", ids)

        xcoords, ycoords = fdsg2d.das.y_aligned.coords[i]
        coords = zeros(2,size(ycoords)[1])
        coords[1,:] .= xcoords
        coords[2,:] .= ycoords
        write(file, "das_y_aligned_coords_$i", coords)
    end
    
    # snapshots (shrinked to physical domain)
    write(file, "snapshots_data",fdsg2d.snapshots.XY[:,:,fdsg2d.domain.inner_id[1], fdsg2d.domain.inner_id[2]])
    write(file, "snapshots_fieldnames",fdsg2d.snapshots.fieldnames)
    write(file, "snapshots_times", fdsg2d.snapshots.t)
    write(file, "snapshots_times_indices",fdsg2d.snapshots.tid)

    close(file)

end



# ------- template --------

"""
Elastic_FDSG.dim2.configtemplate(path)

Creates and empty template configuration file for a 2D ElasticFDSG simulation. 
The user must fill the template afterwards. 

# Arguments
- `path::String`: Path where the template is saved.

# Returns
- `Nothing`: The function saves the template as specified in path.
"""
function configtemplate(path::String)

    template = """
    # This is a template configuration.yaml file for a 2D ElasticFDSG simulation.
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
        dominant_frequency:        # [Hz]
        wavelet_type:              # ricker / gauss1d 
        wavelet_center:            # (should be ≥ 1/fdom) [sec]
        amplitude:                
        force_angle:               # [°]            
        point_source:              # true / false  (if true, double_couple should be false)
        double_couple:             # true / false  (if true, point_source should be false)
        location:
            x:                     # [m]
            y:                     # [m]

    boundaries: # absorbing / else
        left:         # [x-start] 
        right:        # [x-end]    
        top:          # [y-start]  
        bottom:       # [y-end]  

    pml:
        nlayer:                       # 5-15 recommended
        reflection_coefficient:       # 1e-10 recommended

    receivers:
        geophones:
            # - { x: 3000, y: 3000 }        # [m] Example geophone

        das:
            x_aligned:
            #- { x_range: "0:25:10000", y: 4000 }  # [m] Example das
          
            y_aligned:
            #- { x: 1000, y_range: "1000:25:9000"} # [m] Example das
            
            
    snapshots:
        times: 
            # [0.25, 0.5, 0.75, 1.0] # [sec] Example
        fields: 
            # [vx, vy, sxx, sxy, syy] # Example 

    """


    open(path, "w") do io
        write(io, template)
    end
    println("Template saved at $path")
end