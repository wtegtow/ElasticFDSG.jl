# Configurations

Empty template configuration.yaml files can be generated with:
- [ElasticFDSG.config_template()](@ref ElasticFDSG.config_template).
The contents of these files must then be fully completed by the user, either manually or with self-written scripts (see [examples](https://github.com/wtegtow/ElasticFDSG.jl/tree/main/examples)).



## 2D 

```julia 
using ElasticFDSG 

CONFIGFILE_PATH = "path/to/my/config_file.yaml" 
ElasticFDSG.config_template(CONFIGFILE_PATH; dim=2)

```

This saves an empty .yaml file, with hopefully self-explaining configurations:

```yaml
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

```


## 3D 

```julia 
using ElasticFDSG 

CONFIGFILE_PATH = "path/to/my/config_file.yaml"  # add your actual path here 
ElasticFDSG.config_template(configuration_file_path; dim=3)

```

This saves an empty .yaml file, with hopefully self-explaining configurations:

```yaml
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

```