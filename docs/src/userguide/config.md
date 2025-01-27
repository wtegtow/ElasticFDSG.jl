# Configurations

Empty template configuration.yaml files can be generated with:
- [ElasticFDSG.dim2.configtemplate()](@ref ElasticFDSG.dim2.configtemplate) for 2D simulations 
- [ElasticFDSG.dim3.configtemplate()](@ref ElasticFDSG.dim3.configtemplate) for 3D simulations.
The contents of these files must then be fully completed by the user either manually or with self-written scripts.

**Note:** The configuration file reader yet dont perform any type or spell checks. User should fill the form like indicated below. For instance, ```settings.spatial_derivative_order: 10``` will work, while ```"10"``` would throw an error.


## 2D 

```julia 
using ElasticFDSG 

configuration_file_path = joinpath(@__DIR__, "my_2D_configurations.yaml") # add your actual path here 
ElasticFDSG.dim2.configtemplate(configuration_file_path)

```

This creates:

```yaml
# This is a template configuration.yaml file for a 2D ElasticFDSG simulation.
# Any other configuration file can be prepared in the same manner.
# All keywords are case sensitive, and must be named as shown below.
# Users must fill the whole template before running a simulation. 
# The velocity model is prepared in another file.

settings:
    device:                          # cpu / gpu-cuda / gpu-metal 
    precision:                       # Float64 / Float32
    spatial_derivative_order:        # 1-10
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

boundaries: # absorbing / else (reflecting)
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

```


## 3D 

```julia 
using ElasticFDSG 

configuration_file_path = joinpath(@__DIR__, "my_3D_configurations.yaml") # add your actual path here 
ElasticFDSG.dim3.configtemplate(configuration_file_path)

```

This creates:

```yaml
# This is a template configuration.yaml file for a 3D ElasticFDSG simulation.
# Any other configuration file can be prepared in the same manner.
# All keywords are case sensitive, and must be named as shown below.
# Users must fill the whole template before running a simulation.  
# The velocity model is prepared in another file.

settings:
    device:                          # cpu / gpu-cuda / gpu-metal 
    precision:                       # Float64 / Float32
    spatial_derivative_order:        # 1-10
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

boundaries: # absorbing / else (reflecting)
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


```
