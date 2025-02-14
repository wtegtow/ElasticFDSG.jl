using JLD2

function create_3dtest_velmod(comp, vp0, vs0, rho0)

x0 = nothing 
xend = nothing
y0 = nothing 
yend = nothing
z0 = nothing 
zend = nothing

# narrow velocity model to save memory and time  
if comp == "on_vx" || comp == "on_vy"

    x0 = 0
    xend = 1250
    y0 = 0 
    yend = 1250
    z0 = 0
    zend = 200

elseif comp == "on_vz"

    x0 = 0
    xend = 1250
    y0 = 0 
    yend = 200
    z0 = 0
    zend = 1250

end

dx = 10
dy = 10
dz = 10

xcoords = x0:dx:xend
ycoords = y0:dy:yend 
zcoords = z0:dz:zend

nx = length(xcoords)
ny = length(ycoords)
nz = length(zcoords)

# 3D Meshgrid
X = getindex.(Iterators.product(xcoords, ycoords, zcoords), 1)
Y = getindex.(Iterators.product(xcoords, ycoords, zcoords), 2)
Z = getindex.(Iterators.product(xcoords, ycoords, zcoords), 3)

# model dimension
dim = (length(xcoords),length(ycoords),(length(zcoords)))

# allocate empty arrays
vp = zeros(dim);
vs = zeros(dim);
rho = zeros(dim);
eps1 = zeros(dim);
eps2 = zeros(dim);
gam1 = zeros(dim);
gam2 = zeros(dim);
del1 = zeros(dim);
del2 = zeros(dim);
del3 = zeros(dim);


vp[:,:,:]  .= vp0
vs[:,:,:]  .= vs0
rho[:,:,:] .= rho0

# Tsvankin Parameter
#eps1[:,:,:]  .= 0.05
#eps2[:,:,:]  .= 0.1
#gam1[:,:,:]  .= 0.1
#gam2[:,:,:]  .= -0.05
#del1[:,:,:]  .= 0.05
#del2[:,:,:]  .= 0.025
#del3[:,:,:]  .= -0.1

#println("NoN: ", prod(dim))

# velocity model array
velmod = zeros(Float32, 13, nx, ny ,nz);
velmod[1,:,:,:] .= X;
velmod[2,:,:,:] .= Y;
velmod[3,:,:,:] .= Z;
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
path = joinpath(@__DIR__,"velmod.jld2")
jldsave(path;velmod)

end;