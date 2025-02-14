# spatial extends
x0 = 0
xend = 5000
y0 = 0 
yend = 5000
dx = 10
dy = 10

xcoords = x0:dx:xend # x-coordinates
ycoords = y0:dy:yend # y-coordinates

dim = (length(ycoords),(length(xcoords)))

# 2D meshgrid
X = repeat(xcoords', length(ycoords), 1)
Y = repeat(ycoords, 1, length(xcoords))

# velocity
vp = zeros(dim);
vs = zeros(dim);
rho = zeros(dim);
eps0 = zeros(dim);
del0 = zeros(dim);

vp[:,:]  .= 5000;
vs[:,:]  .= 2500;
rho[:,:] .= 2800;

veldim = (7,length(ycoords),(length(xcoords)))

# velocity model array
velmod = zeros(veldim)
velmod[1,:,:] .= X
velmod[2,:,:] .= Y
velmod[3,:,:] .= vp
velmod[4,:,:] .= vs
velmod[5,:,:] .= rho
velmod[6,:,:] .= eps0
velmod[7,:,:] .= del0

# save velocity model
using JLD2
path = joinpath(@__DIR__,"velmod.jld2")
jldsave(path;velmod)