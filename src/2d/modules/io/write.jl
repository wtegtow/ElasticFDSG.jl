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

            if fdsg2d.settings.showinfo println("Results saved.") end
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
    ElasticFDSG.dim2.configtemplate(path)

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
    # The templates is filled with some example placeholder values.
    # The user must fill them before running a simulation. 
    # The velocity model is prepared in another file.

    settings:
        device: cpu                         # cpu / gpu-cuda / gpu-metal 
        precision: Float64                  # Float64 / Float32
        spatial_derivative_order: 10        # (1-10)
        show_summary_in_console: true       # true / false  
        show_progress_in_console: true      # true / false
        save_results: true                  # true / false  
        output:
            destination_folder: path/to/destination/folder
            file_name: cool_simulation 
            format: h5                   # h5 (only supported format yet) 

    time:
        start: 0              # [sec]
        end: 1                # [sec]
        timestep: 0.005       # [sec] (will be checked and changed if unstable)

    source:
        dominant_frequency: 15        # [Hz]
        wavelet_type: gauss1d         # ricker / gauss1d (first derivative of a gaussian)
        wavelet_center: 0.067         # (should be ≥ 1/fdom) [sec]
        amplitude: 1e5                # wavelet amplifier
        location: 
            x: 2500                   # [m]
            y: 2500                   # [m]             

        point_source:             
            use: true               # Enables point source (if true, double_couple.use should be false)
            act_on:
                on_vx: true         # Apply only on vx-component
                on_vy: false        # Apply only on vy-component
                on_vx_and_vy: 
                    force_angle: 0   # [°] Only used if both on_vx and on_vy are enabled. ∈[0°,360°] 0° -> in y negative direction, 180° -> in x negative direction
    
        double_couple:            
            use: false               # Enables double couple source (if true, point_source.use should be false)
            strike: 0                # ∈[0°,360°] 0° -> in y negative direction, 180° -> in x negative direction

    boundaries: # (absorbing / else)
        xstart: absorbing      
        xend: absorbing       
        ystart: absorbing       
        yend: absorbing

    pml:
        nlayer: 10                      # about 10-20 works reasonable
        reflection_coefficient: 1e-5    # about 1e-5 works reasonable

    receivers:
        geophones:
            - { x: 3000, y: 3000 }        # [m] Example geophone

        das:
            x_aligned:
                - { x_range: "0:25:10000", y: 4000 }  # [m] Example das
          
            y_aligned:
                - { x: 1000, y_range: "1000:25:9000"} # [m] Example das
            
        snapshots:
            times: 
                [0.25, 0.5, 0.75, 1.0] # [sec] Example
            fields: 
                [vx, vy, sxx, sxy, syy] # Example 

    """

    open(path, "w") do io
        write(io, template)
    end
    println("Template saved at $path")
end