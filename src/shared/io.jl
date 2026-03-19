function print_h5_tree_rekursive(obj; indent=0, last=true, prefix="", show_header=true)
    if indent == 0 && show_header
    print("""\
    ********************************************************************
    H5-File Tree
    ********************************************************************
    """)
    end

    connector = last ? "└─ " : "├─ "
    child_prefix = last ? "   " : "│  "

    if isa(obj, HDF5.File) || isa(obj, HDF5.Group)
        names = collect(keys(obj))
        for (i, name) in enumerate(names)
            item = obj[name]
            is_last = i == length(names)
            if isa(item, HDF5.Group)
                println(prefix, connector, "Group: ", name)
                print_h5_tree_rekursive(item; indent=indent+1, last=is_last, prefix=prefix*child_prefix, show_header=false)
            elseif isa(item, HDF5.Dataset)
                println(prefix, connector, "Dataset: ", name, "  size=", size(item), "  eltype=", eltype(read(item)))
            end
        end
    end
end

"""
ElasticFDSG.print_h5_tree(filename)
Prints tree structure of a result.h5 file.
# Arguments
- `filename::Stringg`: Path to the `results.h5` file.
"""
function print_h5_tree(filename::String)
    h5open(filename, "r") do file
        print_h5_tree_rekursive(file)
    end
end

function h5_to_dict(obj)
    if isa(obj, HDF5.File) || isa(obj, HDF5.Group)
        d = Dict{String, Any}()
        for k in keys(obj)
            child = obj[k]
            d[string(k)] = h5_to_dict(child)
        end
        return d
    elseif isa(obj, HDF5.Dataset)
        return read(obj)
    else
        return obj
    end
end

"""
ElasticFDSG.load_results(filename)
Loads a result.h5 file into a dictionary.
# Arguments
- `filename::Stringg`: Path to the `results.h5` file.
"""
function load_results(filename::String)
    h5open(filename, "r") do file
        return h5_to_dict(file)
    end
end

"""
ElasticFDSG.config_template(path, dim)

Creates and saves an empty template configuration.yaml file for a ElasticFDSG simulation. 
Users can fill the template afterwards. 

# Arguments
- `path::String`: Path where the template is saved.
- `dim::Real`: dim = 2 => 2D-template, dim=3 => 3D-template

"""
function config_template(path::String; dim::Real=3)

    @assert dim == 2 || dim == 3 "dim != 1 or 2"
    if dim == 2
        config_template2D(path)
    else  
        config_template3D(path)
    end

end;


function config_template2D(path::String)

    template = """
    # This is a template configuration.yaml file for a 2D ElasticFDSG simulation.
    # Any other configuration file can be prepared in the same manner.
    # The user must fill them before running a simulation. 
    # The velocity model is prepared in another file.

    settings:
        device: cpu                         # cpu / cuda / metal / intel / amd  
        precision: Float64                  # Float64 / Float32
        spatial_derivative_order: 4         # 1-10, but 4 recommended 
        show_progress_in_console: true      # true / false 
        output_file: path/to/my/output/file 
    time:
        start: 0              
        end: 1              
        timestep: 0.005  # (will be checked and changed if unstable)

    source:
        dominant_frequency:                
        wavelet_type: ricker               # ricker / gauss1d 
        wavelet_center:                    # (should be ≥ 1.25/fdom) 
        seismic_moment: 
        location:
            x: 0                                                        
            z: 0                          
        moment_tensor:
            Mxx: 0
            Mxz: 0
            Mzz: 0
            anisotropic: false       # true / false 

    boundaries: # (absorbing / else)
        xstart: absorbing      
        xend:   absorbing         
        zstart: absorbing        
        zend:   absorbing
        pml_layer: 10          

    receivers:
        geophones:
            - { x: 0, z: 0 }
            - { x: 0, z: 0 }

        das:
            x_aligned:
                - { x: { start: 0, step: 5, end: 100 }, z: 0 }
                - { x: { start: 0, step: 5, end: 100 }, z: 0 }

            z_aligned:
                - { x: 0, z: { start: 0, step: 5, end: 100 } }

        snapshots:
            # 2D plane snapshots XZ-plane
            fields: ["vx", "vz", "sxx", "sxz", "szz"]
            times: [0, 0.1, 1]
    """

    open(path, "w") do io
        write(io, template)
    end
    println("Configuration template saved at $path")
end


function config_template3D(path::String)

    template = """
    # This is a template configuration.yaml file for a 3D ElasticFDSG simulation.
    # Any other configuration file can be prepared in the same manner.
    # The user must fill them before running a simulation. 
    # The velocity model is prepared in another file.

    settings:
        device: cpu                         # cpu / cuda / metal / intel / amd  
        precision: Float64                  # Float64 / Float32
        spatial_derivative_order: 4         # 1-10, but 4 recommended 
        show_progress_in_console: true      # true / false 
        output_file: path/to/my/output/file 
    time:
        start: 0              
        end: 1              
        timestep: 0.005  # (will be checked and changed if unstable)

    source:
        dominant_frequency:                
        wavelet_type: ricker               # ricker / gauss1d 
        wavelet_center:                    # (should be ≥ 1.25/fdom) 
        seismic_moment: 
        location:
            x: 0                             
            y: 0                            
            z: 0                          
        moment_tensor:
            Mxx: 0
            Mxy: 0
            Mxz: 0
            Myy: 0
            Myz: 0
            Mzz: 0
            anisotropic: false       # true / false 

    boundaries: # (absorbing / else)
        xstart: absorbing      
        xend:   absorbing       
        ystart: absorbing       
        yend:   absorbing         
        zstart: absorbing        
        zend:   absorbing
        pml_layer: 10          

    receivers:
        geophones:
            - { x: 0, y: 0, z: 0 }
            - { x: 0, y: 0, z: 0 }

        das:
            x_aligned:
                - { x: { start: 0, step: 5, end: 100 }, y: 0, z: 0 }
                - { x: { start: 0, step: 5, end: 100 }, y: 0, z: 0 }

            y_aligned: # empty, if no receiver is required

            z_aligned:
                - { x: 0, y: 0, z: { start: 0, step: 5, end: 100 } }

        snapshots:
            # 2D plane snapshots | XY, XZ, YZ - planes centered at plane_positions
            fields: ["vx", "vy", "vz", "sxx", "sxy", "syy", "szz"]
            times: [0, 0.1, 1]
            plane_positions:
                - { x: 100, y: 0, z: 0 }
                - { x: 0, y:100, z: 0 }
        
    """

    open(path, "w") do io
        write(io, template)
    end
    println("Configuration template saved at $path")
end