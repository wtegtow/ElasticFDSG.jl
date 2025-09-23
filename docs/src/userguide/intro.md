# General Usage

The core functionalities of the application are:
- [ElasticFDSG.dim2.runsim()](@ref ElasticFDSG.dim2.runsim) for 2D simulations 
- [ElasticFDSG.dim3.runsim()](@ref ElasticFDSG.dim3.runsim) for 3D simulations.

To run a simulation, the following inputs are required:

(1) a velocity model stored in a julia (.jld2) or numpy array (.npy/npz)

(2) a configuration.yaml file specifying simulation parameters.

Once the inputs are ready, a simulation can be initiated like:

```julia 

using ElasticFDSG

# paths
path_to_configfile = "path/to/my/config_file.yaml"
path_to_velmodfile = "path/to/my/velocity_model.jld2"

# run a 2D simulation
ElasticFDSG.dim2.runsim(path_to_configfile, path_to_velmodfile)

# run a 3D simulation
ElasticFDSG.dim3.runsim(path_to_configfile, path_to_velmodfile)

```


After the calculations are completed, the results will be saved as a .h5 file at the specified location. These results can then be processed using any tool of choice that supports HDF5. 


```julia 
ElasticFDSG.load_results("path_to_my_resulth5_file")
```

In the next sections, the process for creating required velocity models and configuration files will be demonstrated.