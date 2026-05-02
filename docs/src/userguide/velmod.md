# Velocity Models

The velocity model encodes all material properties on a regular Cartesian grid.
It must be provided as a multi-dimensional array — either passed directly as a Julia `AbstractArray`,
or loaded from a `.jld2` (Julia) or `.npy` / `.npz` (NumPy) file.

!!! note
    All elastic parameters are defined at **full integer grid points** (co-located with $v_x$).
    The solver interpolates effective values at staggered positions internally.

---

## 2D velocity model

A 2D model is a **3-dimensional array of shape `(7, nx, nz)`**:

| Index | Field | Unit |
|-------|-------|------|
| `[1,:,:]` | X-coordinate grid | m |
| `[2,:,:]` | Z-coordinate grid | m |
| `[3,:,:]` | P-wave velocity $V_{P0}$ | m/s |
| `[4,:,:]` | S-wave velocity $V_{S0}$ | m/s |
| `[5,:,:]` | Density $\rho$ | kg/m³ |
| `[6,:,:]` | Thomsen parameter $\epsilon$ | — |
| `[7,:,:]` | Thomsen parameter $\delta$ | — |

For an **isotropic** medium set $\epsilon = \delta = 0$.
For a **VTI** medium use the standard Thomsen parameterisation.

### Example — 2D homogeneous VTI model

```julia
using JLD2

h  = 10.0          # grid spacing [m]
xc = range(0.0, step=h, length=300)
zc = range(0.0, step=h, length=300)
nx, nz = length(xc), length(zc)

# 2D coordinate grids  (nx × nz)
X = repeat(xc,          1, nz)
Z = repeat(reshape(zc, 1, :), nx, 1)

vp  = fill(4000.0, nx, nz)   # P-wave velocity [m/s]
vs  = fill(2000.0, nx, nz)   # S-wave velocity [m/s]
rho = fill(2600.0, nx, nz)   # density [kg/m³]
eps = fill(0.1,    nx, nz)   # Thomsen ε
del = fill(0.05,   nx, nz)   # Thomsen δ

velmod = zeros(7, nx, nz)
velmod[1,:,:] .= X
velmod[2,:,:] .= Z
velmod[3,:,:] .= vp
velmod[4,:,:] .= vs
velmod[5,:,:] .= rho
velmod[6,:,:] .= eps
velmod[7,:,:] .= del

# save — the array must be stored under the key "velmod"
jldsave("velmod2d.jld2"; velmod)
```

### Example — 2D layered model

```julia
using JLD2

h = 5.0
xc = range(0.0, step=h, length=500)
zc = range(0.0, step=h, length=400)
nx, nz = length(xc), length(zc)

X = repeat(xc, 1, nz)
Z = repeat(reshape(zc, 1, :), nx, 1)

vp  = zeros(nx, nz)
vs  = zeros(nx, nz)
rho = zeros(nx, nz)

# layer 1: z < 1000 m
mask1 = Z .< 1000
vp[mask1]  .= 3000;  vs[mask1]  .= 1700;  rho[mask1]  .= 2400

# layer 2: z ≥ 1000 m
mask2 = .!mask1
vp[mask2]  .= 4500;  vs[mask2]  .= 2500;  rho[mask2]  .= 2700

velmod = zeros(7, nx, nz)
velmod[1,:,:] .= X;    velmod[2,:,:] .= Z
velmod[3,:,:] .= vp;   velmod[4,:,:] .= vs;  velmod[5,:,:] .= rho
# ε and δ are zero → isotropic

jldsave("velmod2d_layered.jld2"; velmod)
```

---

## 3D velocity model

A 3D model is a **4-dimensional array of shape `(13, nx, ny, nz)`**:

| Index | Field | Unit |
|-------|-------|------|
| `[1,:,:,:]` | X-coordinate grid | m |
| `[2,:,:,:]` | Y-coordinate grid | m |
| `[3,:,:,:]` | Z-coordinate grid | m |
| `[4,:,:,:]` | P-wave velocity $V_{P0}$ | m/s |
| `[5,:,:,:]` | S-wave velocity $V_{S0}$ | m/s |
| `[6,:,:,:]` | Density $\rho$ | kg/m³ |
| `[7,:,:,:]` | Tsvankin $\varepsilon_1$ | — |
| `[8,:,:,:]` | Tsvankin $\varepsilon_2$ | — |
| `[9,:,:,:]` | Tsvankin $\gamma_1$ | — |
| `[10,:,:,:]` | Tsvankin $\gamma_2$ | — |
| `[11,:,:,:]` | Tsvankin $\delta_1$ | — |
| `[12,:,:,:]` | Tsvankin $\delta_2$ | — |
| `[13,:,:,:]` | Tsvankin $\delta_3$ | — |

For an **isotropic** medium set all Tsvankin parameters to zero.
For a **VTI** medium use the equivalences:

| Tsvankin | Thomsen |
|----------|---------|
| $\varepsilon_1 = \varepsilon_2$ | $\varepsilon$ |
| $\gamma_1 = \gamma_2$ | $\gamma$ |
| $\delta_1 = \delta_2$, $\delta_3 = 0$ | $\delta$ |

### Example — 3D homogeneous isotropic model

```julia
using JLD2

h  = 10.0
xc = range(0.0, step=h, length=150)
yc = range(0.0, step=h, length=100)
zc = range(0.0, step=h, length=150)
nx, ny, nz = length(xc), length(yc), length(zc)

# 3D coordinate grids  (nx × ny × nz)
X = getindex.(Iterators.product(xc, yc, zc), 1)
Y = getindex.(Iterators.product(xc, yc, zc), 2)
Z = getindex.(Iterators.product(xc, yc, zc), 3)

vp  = fill(3500.0, nx, ny, nz)
vs  = fill(2000.0, nx, ny, nz)
rho = fill(2500.0, nx, ny, nz)

velmod = zeros(13, nx, ny, nz)
velmod[1,:,:,:] .= X
velmod[2,:,:,:] .= Y
velmod[3,:,:,:] .= Z
velmod[4,:,:,:] .= vp
velmod[5,:,:,:] .= vs
velmod[6,:,:,:] .= rho
# indices 7–13 (Tsvankin parameters) remain zero → isotropic

jldsave("velmod3d.jld2"; velmod)
```

---

## Using NumPy files

The same array can be saved from Python using NumPy:

```python
import numpy as np

# 2D example: shape (7, nx, nz)
np.save("velmod2d.npy", velmod)

# or with named key
np.savez("velmod2d.npz", velmod=velmod)
```

Both `.npy` and `.npz` (with the key `"velmod"`) are supported.
