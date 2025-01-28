# Velocity Models

Velocity models can be prepared in .jld2 (Julia) or .npy/.npz (Python) formats. 
This gives users the freedom to design their velocity models in the tool of their choice and adapt them to the required format. 
The following explains the required structure for 2D and 3D models and shows how to create simple velocity models in Julia.

!!!Note
    To simplify the creation of velocity models, all elastic parameters are defined on full-integer grid points in the staggered grid.


## 2D 
The 2D solver expects an (7, ny, nx) array with:

- 1: X - 2D meshgrid coordinates
- 2: Y - 2D meshgrid coordinates 
- 3: P-wave velocities [m/s]
- 4: S-wave velocities [m/s]
- 5: Densities [kg/m$^{3}$]
- 6: 2D Thomsen Parameter $\epsilon$ 
- 7: 2D Thomsen Parameter $\delta$ 

For an isotropic medium, set 6 & 7 zero.

The following script shows how to create a simple 2D velocity model.
```julia

using JLD2

# spatial extends
x_start = 0
x_end = 10000

y_start = 0 
y_end = 10000

dx = 10 # cell size x-direction 
dy = 10 # cell size y-direction 

xcoords = x_start:dx:x_end # x-coordinates
ycoords = y_start:dy:y_end # y-coordinates

nx = length(xcoords) # number of grid points x-direction
ny = length(ycoords) # number of grid points y-direction

dim = (ny, nx) # model dimensions 

# 2D meshgrid
X = repeat(xcoords', ny, 1)
Y = repeat(ycoords,  1, nx)

vp = zeros(dim);   # P-wave velocity
vs = zeros(dim);   # S-wave velocity
rho = zeros(dim);  # Density
eps0 = zeros(dim); # 2D Thomson parameter epsilon
del0 = zeros(dim); # 2D Thomson parameter delta

# fill arrays with values
vp[:,:]  .= 5000;  
vs[:,:]  .= 2500;
rho[:,:] .= 2800;

# thomson parameter for vti media (zero for isotropic medium)
eps0[:,:] .= 0.2
del0[:,:] .= -0.15

# velocity model array
veldim = (7, ny, nx)
velmod = zeros(veldim)
# fill with elastic properties
velmod[1,:,:] .= X
velmod[2,:,:] .= Y
velmod[3,:,:] .= vp
velmod[4,:,:] .= vs
velmod[5,:,:] .= rho
velmod[6,:,:] .= eps0
velmod[7,:,:] .= del0

# save the velocity model
path = joinpath(@__DIR__,"velmod.jld2") # add your actual path here 
jldsave(path; velmod)
```

## 3D 

The 3D solver expects an (13, nx, ny, nz) array with:

- 1: X - 3D meshgrid coordinates 
- 2: Y - 3D meshgrid coordinates 
- 3: Z - 3D meshgrid coordinates 

- 4: P-wave velocities [m/s] 
- 5: S-wave velocities [m/s] 
- 6: Densities [kg/m$^{3}$] 

Tsvankin Parameter 
- 7: $\varepsilon_1$ 
- 8: $\varepsilon_2$ 
- 9: $\gamma_1$ 
- 10: $\gamma_2$ 
- 11: $\delta_1$ 
- 12: $\delta_2$
- 13: $\delta_3$ 

For an isotropic medium, set 7-13 to zero.

For a VTI medium, the relationship between Tsvankin and Thomsen Parameters can be used:

| Tsvankin Parameter      | `<=>` | Thomsen Parameter |
|-------------------------|--------------|-------------------|
| $\varepsilon_1 = \varepsilon_2$ | `<=>`         | $\varepsilon$       |
| $\gamma_1 = \gamma_2$            | `<=>`         | $\gamma$            |
| $\delta_1 = \delta_2$, $\delta_3 = 0$             | `<=>`         | $\delta$            |




The following script shows how to create a simple 3D velocity model.

```julia

using JLD2

# spatial extends
x_start = 0
x_end = 500

y_start = 0 
y_end = 500

z_start = 0
z_end = 500

dx = 5 # cell size x-direction 
dy = 5 # cell size y-direction 
dz = 5 # cell size z-direction 

xcoords = x_start:dx:x_end # x-coordinates
ycoords = y_start:dy:y_end # y-coordinates
zcoords = z_start:dz:z_end # z-coordinates

nx = length(xcoords) # number of grid points x-direction
ny = length(ycoords) # number of grid points y-direction
nz = length(zcoords) # number of grid points z-direction

# 3D Meshgrid
X = getindex.(Iterators.product(xcoords, ycoords, zcoords), 1)
Y = getindex.(Iterators.product(xcoords, ycoords, zcoords), 2)
Z = getindex.(Iterators.product(xcoords, ycoords, zcoords), 3)

dim = (nx, ny, nz) # model dimension

vp = zeros(dim);    # P-wave velocity
vs = zeros(dim);    # S-wave velocity
rho = zeros(dim);   # Density velocity
eps1 = zeros(dim);  # Tsvankin parameter epsilon 1 
eps2 = zeros(dim);  # Tsvankin parameter epsilon 2
gam1 = zeros(dim);  # Tsvankin parameter gamma 1 
gam2 = zeros(dim);  # Tsvankin parameter gamma 2 
del1 = zeros(dim);  # Tsvankin parameter delta 1 
del2 = zeros(dim);  # Tsvankin parameter delta 2 
del3 = zeros(dim);  # Tsvankin parameter delta 3 

# fill arrays with values
vp[:,:,:] .= 5000;
vs[:,:,:]  .= 2500;
rho[:,:,:] .= 2800;

# Tsvankin Parameter
eps1[:,:,:]  .= 0.05
eps2[:,:,:]  .= 0.1
gam1[:,:,:]  .= 0.1
gam2[:,:,:]  .= -0.05
del1[:,:,:]  .= 0.05
del2[:,:,:]  .= 0.025
del3[:,:,:]  .= -0.1

# velocity model array
veldim = (13, nx, ny, nz)
velmod = zeros(veldim);
# fill with elastic properties
velmod[1,:,:,:] .= X
velmod[2,:,:,:] .= Y
velmod[3,:,:,:] .= Z

velmod[4,:,:,:] .= vp
velmod[5,:,:,:] .= vs
velmod[6,:,:,:] .= rho

velmod[7,:,:,:] .= eps1
velmod[8,:,:,:] .= eps2
velmod[9,:,:,:] .= gam1
velmod[10,:,:,:] .= gam2
velmod[11,:,:,:] .= del1
velmod[12,:,:,:] .= del2
velmod[13,:,:,:] .= del3

# save velocity model in jld2 file
using JLD2
path = joinpath(@__DIR__,"velmod.jld2") # add your actual path here 
jldsave(path; velmod)
```