## Validation

To validate the correctness of the application, the numerically derived seismometer data is compared with analytical solutions for a given homogeneous medium.

Solutions for inhomogeneous partial differential equations can be obtained using Green's functions 
$G(\mathbf{x}, t; \mathbf{x_{0}}, t_{0})$ with $\delta$-distributions as source terms acting on $(\mathbf{x}, t)$ and activated on $(\mathbf{x_{0}}, t_{0})$.

The solution to that problem leads to the practical scenario where the displacement field can be described by a convolution of the Green's function with the source-time function:

```math

u_{i} = G * S, \quad i \in \{x, y, z\}.

```

The velocity can then be derived as:

```math

v_{i} = \frac{\partial u_{i}}{\partial t}, \quad i \in \{x, y, z\}.

```

In the 2D case, the Green's function is given by:

```math
G_{2D}(\mathbf{x},t) = \frac{1}{2\pi \rho c^2} \frac{H\biggl((t-t_{0})-\frac{r}{c}\biggr)}{\sqrt{(t-t_{0})^2-\frac{r^2}{c^2}}},

```

and in the 3D case:


```math

G_{3D}(\mathbf{x},t) = \frac{1}{4 \pi \rho c^2 r} \delta(t - \frac{r}{c})


```

In the equations above, $\rho$ represents the density, $c$ is the wave speed, $H$ is the Heaviside step function, and $\delta$ is the Dirac delta function and $r$ the Euclidian distance for 2 and 3D, respectivly.

The following images show a comparison of numerical seismograms with corresponding analytical solutions (red dashed lines). 


2D:
![comp](assets/comparision.png)


3D: coming soon



The examples were selected to illustrate grid-dispersion.
Grid dispersion can be counteracted with higher order spatial derivative operators $N$. However, users should set the grid-spacing sufficient fine to avoid grid-dispersion.

For heterogeneous media, analytical solutions are not so easily obtainable. 
Furthermore, due to different schemes, results from different numerical methods are also more difficult to compare. 
For that reason, the "validation" of heterogeneous media is yet left completely on a visible basis:

![2danim](assets/2dvalid.gif)