# Configurations

Empty template configuration.yaml files can be generated with:
- [ElasticFDSG.dim2.configtemplate()](@ref ElasticFDSG.dim2.configtemplate) for 2D simulations 
- [ElasticFDSG.dim3.configtemplate()](@ref ElasticFDSG.dim3.configtemplate) for 3D simulations.
The contents of these files must then be fully completed by the user either manually or with self-written scripts.

!!! note

    The configuration file reader yet dont perform any type or spell checks. User should fill the form like indicated below to avoid mysterious errors further downstream in the application. For instance, ```settings.spatial_derivative_order: 10``` will work, while ```"10"``` would throw an error.


## 2D 

```julia 
using ElasticFDSG 

configuration_file_path = "path/to/my/config_file.yaml" 
ElasticFDSG.dim2.configtemplate(configuration_file_path)

```

This creates:

```yaml
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

```


## 3D 

```julia 
using ElasticFDSG 

configuration_file_path = "path/to/my/config_file.yaml"  # add your actual path here 
ElasticFDSG.dim3.configtemplate(configuration_file_path)

```

This creates:

```yaml
# This is a template configuration.yaml file for a 3D ElasticFDSG simulation.
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
    dominant_frequency: 25                 # [Hz]
    wavelet_type: ricker                   # ricker / gauss1d 
    wavelet_center: 0.04                   # (should be ≥ 1/fdom) [sec]
    amplitude: 1 
    location:
        x: 200                             # [m]
        y: 200                             # [m]
        z: 100                             # [m]      
    point_source:                 
        use: true                          # true / false (if true, double_couple.use should be false)
        act_on:
            on_vx: true                    # Apply only on vx-component 
            on_vy: false                   # Apply only on vy-component
            on_vz: false                   # Apply only on vz-component
            on_all:                        # Only activates, if all v-components are enabled true
                phi: 0                     # [°] Elevation angle
                theta: 0                   # [°] Azimuth angle (x-y plane)
    double_couple: 
        use: false                         # true / false  (if true, point_source.use should be false)
        strike: 0                          # [°] ∈ [0, 360]
        dip: 0                             # [°] ∈ [0, 90]
        rake: 0                            # [°] ∈ [0, 360]
        anisotropic_moment_tensor: true    # true / false

boundaries: # (absorbing / else)
    xstart: absorbing      
    xend:   absorbing       
    ystart: absorbing       
    yend:   absorbing         
    zstart: absorbing        
    zend:   absorbing          

pml:
    nlayer: 10                       # about 10-20 works reasonable
    reflection_coefficient: 1e-15    # about 1e-5 works reasonable

receivers:
    geophones: 
        - { x: 100, y: 200, z: 300 }    #[m] Example geophone

    das:
        x_aligned: 
            - { x_range: "0:2.5:500", y: 100, z: 100 }  #[m] Example das 
        
        y_aligned: 
            - { x: 200, y_range: "100:2.5:400", z: 100 } #[m] Example das 

        z_aligned:
            - { x: 300, y: 250, z_range: "0:2.5:400" } #[m] Example das 
                
                
    snapshots: # 2D plane-snapshots: XY, XZ, YZ planes only
        times: 
             [0.25, 0.5, 0.75, 1.0] # [sec] Example
        fields: 
             [vx, vy, vz, sxx, sxy, sxz, syy, syz, szz] # Example 
        origins:
             { x: 250, y: 250, z: 250 } # [m] Coordinates of snapshot origin. Example: XY-plane snapshots will be at z=250m

```