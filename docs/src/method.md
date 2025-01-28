## Theory
This section provides a brief overview of the theoretical basis of the application.
For more comprehensive descriptions, readers are referred to literature listed below.

Elastic wave propagation in a 3D anisotropic linear elastic material can be described by the following set of first order partial differential equations:

```math
\begin{aligned}
\partial_{t} \sigma_{xx} &=  c_{11} \partial_{x} v_{x} + c_{12} \partial_{y} v_{y} + c_{13} \partial_{z} v_{z}  \\ 
\partial_{t} \sigma_{yy} &=  c_{12} \partial_{x} v_{x} + c_{22} \partial_{y} v_{y} + c_{23} \partial_{z} v_{z}  \\
\partial_{t} \sigma_{zz} &=  c_{13} \partial_{x} v_{x} + c_{23} \partial_{y} v_{y} + c_{33} \partial_{z} v_{z}  \\
\partial_{t} \sigma_{xy} &=  c_{66} (\partial_{x} v_{y} + \partial_{y} v_{x})  \\
\partial_{t} \sigma_{xz} &=  c_{55} (\partial_{x} v_{z} + \partial_{z} v_{x})  \\
\partial_{t} \sigma_{yz} &=  c_{44} (\partial_{y} v_{z} + \partial_{z} v_{y}) \tag{1} \\

\end{aligned}
```

```math
\begin{aligned}
\rho \partial_{t} v_{x} &= \partial_{x} \sigma_{xx} + \partial_{y} \sigma_{xy} + \partial_{z} \sigma_{xz} + f_{x} \\
\rho \partial_{t} v_{y} &= \partial_{x} \sigma_{xy} + \partial_{y} \sigma_{yy} + \partial_{z} \sigma_{yz} + f_{y} \\
\rho \partial_{t} v_{z} &= \partial_{x} \sigma_{xz} + \partial_{y} \sigma_{xz} + \partial_{z} \sigma_{zz} + f_{z} \tag{2} \\
\end{aligned}
```

In these equations:
- are $\sigma_{xx}, \sigma_{yy}, \sigma_{zz}$ the normal stress components.
- are $\sigma_{xy}, \sigma_{xz}, \sigma_{yz}$ the shear stress components.
- are $c_{11}, c_{12}, c_{13}, c_{22}, c_{23}, c_{33}, c_{44}, c_{55}, c_{66}$ the 9 elastic constants needed to describe an orthorhombic material.
- is $\rho$ is the density of the material.
- are $v_{x}, v_{y}, v_{z}$ the particle velocities in the $x$, $y$, and $z$ directions, respectively.
- are $f_{x}, f_{y}, f_{z}$ the external body forces acting in the $x$, $y$, and $z$ directions, respectively.

Using the Tsvankin notation, isotropic, vertical transversal isotropic (VTI) and ortorhombic materials (ORT) can be characterized by two vertical velocities and 7 dimensionless paramters $\epsilon_{1}, \epsilon_{2}, \gamma_{1}, \gamma_{2}, \delta_{1}, \delta_{2}, \delta_{3}$ :

```math
\begin{aligned}
c_{33} &= v_{p0}^2 \cdot \rho \\
c_{55} &= v_{s0}^2 \cdot \rho \\
c_{11} &= (2 \epsilon_{2} + 1) \cdot c_{33} \\
c_{22} &= c_{33} \cdot (2 \epsilon_{1} + 1) \\
c_{66} &= c_{55} \cdot (2 \gamma_{1} + 1) \\
c_{44} &= \frac{c_{66}}{1 + \gamma_{2}} \\
c_{13} &= \sqrt{2 c_{33} \cdot (c_{33} - c_{55}) \cdot \delta_{2} + (c_{33} - c_{55})^2} - c_{55} \\
c_{23} &= \sqrt{2 c_{33} \cdot (c_{33} - c_{44}) \cdot \delta_{1} + (c_{33} - c_{44})^2} - c_{44} \\
c_{12} &= \sqrt{2 c_{11} \cdot (c_{11} - c_{66}) \cdot \delta_{3} + (c_{11} - c_{66})^2} - c_{66}
\end{aligned}
```



## Numerical Scheme 

The set of equations (1) and (2) can be solved using a finite-difference staggered-grid scheme. In a staggered grid, field quantities are not co-located on same grid points but are distributed across predefined grid cells. The primary advantage of a staggered grid is an enhanced accuracy of spatial central difference operators, as well as staggered temporal finite difference operators, which also improves the accuracy of the time marching.

However, the distribution of field quantities within a grid cell also has drawbacks, especially during the processing of the results. For example, users working with geophone data must consider that the velocity components $v_{x}, v_{y}, v_{z}$ are not measured at the same location.

The staggered grid scheme used in the application is illustrated below:


          ↑ y                           │      y              △             ■ vₓ  i,j,k
            │         i+½,j+½           │  z↑ /              /              ◆ vᵧ  i+½,j+½,k
      i,j+½ │         │vᵧ               │   │/              /               ● v𝑧  i+½,j,k+½ 
       σₓᵧ ─◇─────────◆───              │   ○ ............ ●
            │         │                 │   ╎              ╎                □ σₓₓ,σᵧᵧ,σ𝑧𝑧 i+½,j,k
            │         │                 │   ╎              ╎                ○ σₓ𝑧 i,j,k+½
            │         │                 │   ╎   ◇          ╎  ◆             ◇ σₓᵧ i,j+½,k
            ■─────────□──── ──→ x       │   ╎  /           ╎ /              △ σᵧ𝑧 i+½,j+½,k+½
            vₓ      σₓₓ,σᵧᵧ             │   ╎ /            ╎/  
           i,j       i+½,j              │   ■ ............ □ ─→ x 
                                        │ 

Using this scheme, the discrete form of (1) and (2) are given by:

```math
\begin{aligned}

\sigma_{xx \; (i+\frac{1}{2},j,k)}^{t_{n}+\frac{1}{2}} &= \sigma_{xx \; (i+\frac{1}{2},j,k)}^{t_{n}} \Delta t (c_{11} \mathcal{D}_{x} v_{x} + c_{12} \mathcal{D}_{y} v_{y} + c_{13} \mathcal{D}_{z} v_{z} ) \bigg|_{(i+\frac{1}{2},j,k)} \\ 

\sigma_{yy \; (i+\frac{1}{2},j,k)}^{t_{n}+\frac{1}{2}} &= \sigma_{yy \; (i+\frac{1}{2},j,k)}^{t_{n}} \Delta t (c_{12} \mathcal{D}_{x} v_{x} + c_{22} \mathcal{D}_{y} v_{y} + c_{23} \mathcal{D}_{z} v_{z} ) \bigg|_{(i+\frac{1}{2},j,k)} \\ 

\sigma_{zz \; (i+\frac{1}{2},j,k)}^{t_{n}+\frac{1}{2}} &= \sigma_{zz \; (i+\frac{1}{2},j,k)}^{t_{n}} \Delta t (c_{13} \mathcal{D}_{x} v_{x} + c_{23} \mathcal{D}_{y} v_{y} + c_{33} \mathcal{D}_{z} v_{z} ) \bigg|_{(i+\frac{1}{2},j,k)} \\ 

\sigma_{xy \; (i,j+\frac{1}{2},k)}^{t_{n}+\frac{1}{2}} &= \sigma_{xy \; (i,j+\frac{1}{2},k)}^{t_{n}} \Delta t \, c_{66} ( \mathcal{D}_{x} v_{y} +  \mathcal{D}_{y} v_{x})\bigg|_{(i,j+\frac{1}{2},k)} \\

\sigma_{xz \; (i,j,k+\frac{1}{2})}^{t_{n}+\frac{1}{2}} &= \sigma_{xz \; (i,j,k+\frac{1}{2})}^{t_{n}} \Delta t \, c_{55} ( \mathcal{D}_{x} v_{z} + \mathcal{D}_{z} v_{x})\bigg|_{(i,j,k+\frac{1}{2})} \\

\sigma_{yz \; (i+\frac{1}{2},j+\frac{1}{2},k+\frac{1}{2})}^{t_{n}+\frac{1}{2}} &= \sigma_{yz \; (i+\frac{1}{2},j+\frac{1}{2},k+\frac{1}{2})}^{t_{n}} \Delta t \, c_{44} ( \mathcal{D}_{y} v_{z} +  \mathcal{D}_{y} v_{z})\bigg|_{(i+\frac{1}{2},j+\frac{1}{2},k+\frac{1}{2})} \\
            
            
\end{aligned}
```      

```math
\begin{aligned}
v_{x \; (i,j,k)}^{t_{n}+1} &= v_{x \; (i,j,k)}^{t_{n}} \frac{\Delta t}{\rho} (\mathcal{D}_{x} \sigma_{xx} + \mathcal{D}_{y} \sigma_{xy} + \mathcal{D}_{z} \sigma_{xz} + f_{x}) \bigg|_{(i,j,k)} \\
v_{y \; (i+\frac{1}{2},j+\frac{1}{2},k)}^{t_{n}+1} &= v_{y \; (i+\frac{1}{2},j+\frac{1}{2},k)}^{t_{n}} \frac{\Delta t}{\rho} (\mathcal{D}_{x} \sigma_{xy} + \mathcal{D}_{y} \sigma_{yy} + \mathcal{D}_{z} \sigma_{yz} + f_{y}) \bigg|_{(i+\frac{1}{2},j+\frac{1}{2},k)} \\
v_{z \; (i+\frac{1}{2},j,k+\frac{1}{2})}^{t_{n}+1} &= v_{z \; (i+\frac{1}{2},j,k+\frac{1}{2})}^{t_{n}} \frac{\Delta t}{\rho} (\mathcal{D}_{x} \sigma_{xz} + \mathcal{D}_{y} \sigma_{xz} + \mathcal{D}_{z} \sigma_{zz} + f_{z}) \bigg|_{(i+\frac{1}{2},j,k+\frac{1}{2})} \tag{3}

\end{aligned}
```
Here $(i,j,k)$ represent the grid points in the x,y,z-directions, respectively; 
$t_{n}$ denotes the $n$-th time step; 
and $\Delta t$ the time increment.

The differential operators $\mathcal{D}$ are given by central difference approximations:

```math
\begin{aligned}

\mathcal{D_{x}} f = \sum_{n=1}^{N} \frac{c_{n}}{\Delta x} ( f(x_{i+n}, y_{j}, z_{k}) - f(x_{i- n}, y_{j}, z_{k})) \\
\mathcal{D_{y}} f = \sum_{n=1}^{N} \frac{c_{n}}{\Delta y} ( f(x_{i}, y_{j+n}, z_{k}) - f(x_{i}, y_{j-n}, z_{k})) \\
\mathcal{D_{z}} f = \sum_{n=1}^{N} \frac{c_{n}}{\Delta z} ( f(x_{i}, y_{j}, z_{k+n}) - f(x_{i}, y_{j}, z_{k-n}))
        
\end{aligned}
```      
with order $N$ and coefficients $c_{n}$.

!!! note

  To compute a spatial derivative of order $N$ at a specific location, at least $N$ neighboring nodes are required on both sides. As a result, for edge nodes no spatial derivatives can be calculated. These edge nodes (ghost nodes) are effective model boundaries and should lie outside the physical (user-defined) domain. The application automatically extends the user-defined domain with $N$ ghost node layers.

  In staggered grids, it is often beneficial to assign certain elastic properties to specific points within the grid cell. However, requiring users to define such grids can become difficult to manage, particularly for complex media. To simplify this, we define all elastic properties at full integer grid points (e.g., on $v_{x}$) and effective properties are calculated by interpolating values from neighboring nodes.


## Sources 
Earthquake simulations require the excitation of point or double-couple forces using the body force term in the equation of motion.

- Point Sources:
Point sources can directly be applied to the velocity components at the desired source location $(s_{x},s_{y},s_{y})$, using elevation and azimuth angles combined with a source-time function that contains a wavelet. However, due to the positions of velocity components in the grid cell, the excitation occurs at slightly different locations.

- Double couple sources
Double-couple sources involve applying force pairs $(p,q)$ with strength $M_{pq}$ to the velocity field at a desired source location. In the used staggered grid scheme, double-couple sources are centered around normal stresses grid points, i.e., $(i+\frac{1}{2},j,k)$.
As a result, a defined double couple source location is centered shifted by half a grid point in the x-direction! 
Using moment tensors to represent the equivalent distribution of force pairs, at least 30 force contributions need to be applied to the surounding velocity field.
The 10 force components contributing to $f_{x}$ are:

```math
\begin{aligned}
f_{x \; i,j,k} &= - f_{x \; i+1,j,k} &= \frac{M_{xx}(t)}{\Delta x^{2} \Delta y \Delta z} \\\\

f_{x \; i,j+1,k} &= - f_{x \; i,j-1,k} &= \frac{M_{xy}(t)}{4 \Delta x \Delta y^{2} \Delta z} \\
f_{x \; i+1,j+1,k} &= - f_{x \; i+1,j-1,k} &= \frac{M_{xy}(t)}{4 \Delta x  \Delta y^{2}  \Delta z} \\\\

f_{x \; i,j,k+1} &= - f_{x \; i,j,k-1} &= \frac{M_{xz}(t)}{4 \Delta x \Delta y \Delta z^{2}} \\
f_{x \; i+1,j,k+1} &= - f_{x \;i+1,j,k-1} &= \frac{M_{xz}(t)}{4 \Delta x \Delta y \Delta z^{2}} \\

\end{aligned}
```  

The same logic applies to the body force components contributing to $f_{y}$ and $f_{z}$.

The moment tensor components are computed by user defined dip, strike and rake values for the following coordinate system:
  
         x                  
       /                     
      /                        δ = Dip, λ = Rake, Φ = Strike                
      ---- y                   Coordinate System: x -> north,            
     |                                            y -> east,
     |                                            z -> positive downward
     z  

```math
\begin{aligned}
M_{xx} &= -(sin(\delta) \, cos(\lambda) \, sin(2\Phi) + sin(2\delta) \, sin(\lambda) \, sin(\Phi)^{2}) \\
M_{xy} &= \sin(\delta) \, \cos(\lambda) \, \cos(2\Phi) + \frac{1}{2} \sin(2\delta) \, \sin(\lambda) \, \sin(2\Phi) \\
M_{xz} &= -(\cos(\delta) \, \cos(\lambda) \, \cos(\Phi) + \cos(2\delta) \, \sin(\lambda) \, \sin(\Phi)) \\
M_{yy} &= \sin(\delta) \, \cos(\lambda) \, \sin(2\Phi) - \sin(2\delta) \, \sin(\lambda) \, \cos(\Phi)^2 \\
M_{yz} &= -(\cos(\delta) \, \cos(\lambda) \, \sin(\Phi) - \cos(2\delta) \, \sin(\lambda) \, \cos(\Phi)) \\
M_{zz} &= \sin(2\delta) \, \sin(\lambda)
\end{aligned}
```  

    
## C-PML 

Earthquake simulations often require modeling wave propagation in unbounded media. Perfectly Matched Layer (PML) is a very effective method to prevent artificial reflections at model boundaries.
The underlying idea is to manipulate the wave equation to obtain exponentially decaying plane wave solutions for complex arguments.
Assuming that the computational domain (real arguments) is surrounded by a complex region (the PML region), the amplitudes of incident waves decay and cause negligible reflections at the model boundaries.


Komatitsch & Martin (2007) introduced a memory efficient convolution-based unsplit PML formulation (C-PML). 
This approach requires storing one additional memory variable for each spatial derivative, but only in the PML region.
For the PML region, spatial derivative operators are replaced by:

```math
\begin{aligned}
\mathcal{D}_{\hat{x}} &= \frac{\mathcal{D_{x}}}{\kappa_{x}} + \Psi_{x} \tag{4}
\end{aligned}
```
Here, $\Psi_{x}$ represents the memory variable associated with the field from which the derivative is taken. The memory variable $\Psi_{x}$ is updated at each time step by:

```math
\begin{aligned}
\Psi_{x}^{t_{n}} &= b_{x} \Psi_{x}^{t_{n}-1} + a_{x} (\mathcal{D_{x}})^{t_{n} + \frac{1}{2}} \tag{5}
\end{aligned}
```

while $a_{x}$, $b_{x}$ $\kappa_{x}$ are precomputed PML related parameter. 

Accordingly, differential operators $\mathcal{D_{y}}$ and $\mathcal{D_{z}}$ are replaced by $\mathcal{D}_{\hat{y}}$ and $\mathcal{D}_{\hat{z}}$ in the PML-region.


## References 

Komatitsch, D., & Martin, R. (2007). An unsplit convolutional perfectly matched layer improved at grazing incidence for the seismic wave equation. Geophysics, 72(5), SM155-SM167.

Moczo, P., Kristek, J., & Gális, M. (2014). The finite-difference modelling of earthquake motions: Waves and ruptures. Cambridge University Press.

Tsvankin, I. (1997). Anisotropic parameters and P-wave velocity for orthorhombic media. Geophysics, 62(4), 1292-1309.

Virieux, J. (1984). SH-wave propagation in heterogeneous media: Velocity-stress finite-difference method. Geophysics, 49(11), 1933-1942. 