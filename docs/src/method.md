# Method

This section gives a concise overview of the numerical method implemented in ElasticFDSG.
For a thorough mathematical treatment readers are referred to
Virieux (1986), Levander (1988), and Komatitsch & Martin (2007).

---

## Governing equations

Elastic wave propagation in a linear anisotropic medium is governed by the first-order
velocity–stress system:

```math
\begin{aligned}
\rho \, \partial_t v_x &= \partial_x \sigma_{xx} + \partial_y \sigma_{xy} + \partial_z \sigma_{xz} \\
\rho \, \partial_t v_y &= \partial_x \sigma_{xy} + \partial_y \sigma_{yy} + \partial_z \sigma_{yz} \\
\rho \, \partial_t v_z &= \partial_x \sigma_{xz} + \partial_y \sigma_{yz} + \partial_z \sigma_{zz}
\end{aligned}
```

```math
\begin{aligned}
\partial_t \sigma_{xx} &= C_{11}\,\partial_x v_x + C_{12}\,\partial_y v_y + C_{13}\,\partial_z v_z \\
\partial_t \sigma_{yy} &= C_{12}\,\partial_x v_x + C_{22}\,\partial_y v_y + C_{23}\,\partial_z v_z \\
\partial_t \sigma_{zz} &= C_{13}\,\partial_x v_x + C_{23}\,\partial_y v_y + C_{33}\,\partial_z v_z \\
\partial_t \sigma_{xy} &= C_{66}\!\left(\partial_x v_y + \partial_y v_x\right) \\
\partial_t \sigma_{xz} &= C_{55}\!\left(\partial_x v_z + \partial_z v_x\right) \\
\partial_t \sigma_{yz} &= C_{44}\!\left(\partial_y v_z + \partial_z v_y\right)
\end{aligned}
```

The nine independent stiffness coefficients $C_{11}, C_{12}, C_{13}, C_{22}, C_{23}, C_{33}, C_{44}, C_{55}, C_{66}$
describe the most general symmetry supported: **orthorhombic (ORT)**.
VTI and isotropic models are special cases.

### Stiffness from Tsvankin parameters

For 3D models, the stiffness tensor is assembled from the Tsvankin parameterisation:

```math
\begin{aligned}
C_{33} &= V_{P0}^2 \, \rho \\
C_{55} &= V_{S0}^2 \, \rho \\
C_{11} &= (2\varepsilon_2 + 1)\,C_{33} \\
C_{22} &= (2\varepsilon_1 + 1)\,C_{33} \\
C_{66} &= (2\gamma_1   + 1)\,C_{55} \\
C_{44} &= C_{66} / (1 + \gamma_2) \\
C_{13} &= \sqrt{2C_{33}(C_{33}-C_{55})\,\delta_2 + (C_{33}-C_{55})^2} - C_{55} \\
C_{23} &= \sqrt{2C_{33}(C_{33}-C_{44})\,\delta_1 + (C_{33}-C_{44})^2} - C_{44} \\
C_{12} &= \sqrt{2C_{11}(C_{11}-C_{66})\,\delta_3 + (C_{11}-C_{66})^2} - C_{66}
\end{aligned}
```

For 2D models the standard Thomsen parameterisation is used instead ($\epsilon$, $\delta$ — see fwm.tex for details).

---

## Staggered-grid discretisation

Field components are not co-located but distributed across the unit cell following the
standard staggered-grid layout of Virieux (1986):

```
3D cell (one octant shown):

      ↑ z               field           grid position
      │                 ─────────────────────────────────────
   k  ● ─────── ■       σxx, σyy, σzz   (i,   j,   k)
      ╎         ╎       vx              (i-½, j,   k)
      ╎    ◆    ╎       vy              (i,   j-½, k)
      ╎         ╎       vz              (i,   j,   k-½)
      □ ─────── ○ ── → x   σxy          (i-½, j-½, k)
     /                  σxz             (i-½, j,   k-½)
    / ── y               σyz            (i,   j-½, k-½)
  i,j
```

The discrete **leapfrog** update scheme advances velocities and stresses alternately
at half-integer time steps ($t^{n}$ and $t^{n+\frac{1}{2}}$):

**Velocity update** (at time $t^{n+\frac{1}{2}}$):

```math
v_{x(I)}^{n+\frac{1}{2}} = v_{x(I)}^{n-\frac{1}{2}}
    + \frac{\Delta t}{\rho_{(I)}}
      \!\left(\partial_x \sigma_{xx} + \partial_y \sigma_{xy} + \partial_z \sigma_{xz}\right)\!\bigg|_{(I)}^{n},
\quad I = \!\left(i{-}\tfrac{1}{2},j,k\right)
```

(analogously for $v_y$ and $v_z$).

**Stress update** (at time $t^{n+1}$):

```math
\sigma_{xx(I)}^{n+1} = \sigma_{xx(I)}^{n}
    + \Delta t \!\left(C_{11}\,\partial_x v_x + C_{12}\,\partial_y v_y + C_{13}\,\partial_z v_z\right)\!\bigg|_{(I)}^{n+\frac{1}{2}},
\quad I = (i,j,k)
```

(analogously for the remaining five stress components).

Spatial derivatives are approximated by central differences of order $N$:

```math
\partial_x f \big|_{(i,j,k)} \approx \sum_{n=1}^{N} \frac{k_n}{\Delta x}
    \Bigl(f(x_{i+n},y_j,z_k) - f(x_{i-n},y_j,z_k)\Bigr)
```

with coefficients $k_n$ optimised for maximum wavenumber accuracy.
The combination of staggered spatial and temporal grids yields **second-order accuracy in both space and time**, $\mathcal{O}(\Delta t^2, \Delta x^{2N})$, without requiring intermediate storage states.

Since density and stiffness are defined at integer nodes (not staggered), effective values
at the staggered positions are obtained by averaging the values from neighboring grid points.

---

## Domain extension

Starting from the user-defined physical domain of size $n_x \times n_y \times n_z$, the grid is
extended in two stages:

1. **PML layers** — $N_\mathrm{PML}$ absorbing cells appended at each absorbing boundary.
2. **Ghost nodes** — $N_g = N$ additional cells beyond the PML to support the finite-difference
   stencil at every interior point.

Material properties in the extended region are filled by nearest-neighbour replication from
the physical domain, avoiding artificial impedance contrasts.

---

## Absorbing boundaries — C-PML

Without absorbing boundaries, waves reflect off the computational domain edges.
ElasticFDSG uses the **convolutional perfectly matched layer (C-PML)** of
Komatitsch & Martin (2007), which replaces each spatial derivative $\partial_x$
with the modified operator:

```math
\tilde{\partial}_x f = \partial_x f + \psi_x
```

where $\psi_x$ is a memory variable updated each time step:

```math
\psi_x^{n+1} = b_x \, \psi_x^{n} + a_x \, \partial_x f^{n+1}
```

The precomputed coefficients are:

```math
b_x = \exp\!\bigl(-(d_x + \alpha_x)\Delta t\bigr), \qquad
a_x = \frac{d_x}{d_x + \alpha_x}(b_x - 1)
```

The damping profile $d_x$ increases polynomially from zero at the physical/PML interface
to a maximum at the outer boundary, while the frequency-shift parameter $\alpha_x$ is set to
$\pi f_\mathrm{dom}$.
This formulation requires **18 additional memory variables** in 3D (6 per spatial direction),
stored only inside the PML region.

!!! warning "PML stability for anisotropic media"
    The C-PML is intrinsically unstable for anisotropic materials whose quasi-shear slowness
    surface is non-convex (Bécache et al. 2003).
    This includes many shale formations with high anellipticity.
    In practice, sufficiently thick PML regions typically attenuate such modes before
    they re-enter the physical domain, but instabilities may still appear in extreme cases.

---

## Moment tensor source

Sources are injected as incremental stresses following Shi et al. (2018).
The source time function derivative $\dot{s}(t)$ is scaled by the scalar seismic moment $M_0$
and distributed over the stress components:

```math
\sigma_{xx(i,j,k)}^{n+1} \mathrel{-}= \frac{\Delta t}{\Delta x \Delta y \Delta z} M_0 \hat{M}_{xx} \dot{s}(t_{n+1})
```

The normal stress components are applied at a single node $(i,j,k)$.
Each shear component is symmetrically distributed over **four** neighboring nodes:

```math
\sigma_{xy(I)}^{n+1} \mathrel{-}=
\frac{\Delta t}{4 \Delta x \Delta y \Delta z} M_0 \hat{M}_{xy} \dot{s}(t_{n+1}),
\quad I = \!\left(i{\pm}\tfrac{1}{2},\, j{\pm}\tfrac{1}{2},\, k\right)
```

(and analogously for $\sigma_{xz}$ and $\sigma_{yz}$).
This symmetric distribution produces the correct radiation pattern for each moment tensor configuration.

### Moment tensor from fault geometry

For shear-slip (double-couple) sources, the moment tensor components can be computed
from strike $\Phi$, dip $\delta$, and rake $\lambda$ (coordinate system: $x$→ North, $y$→ East, $z$↓):

```math
\begin{aligned}
M_{xx} &= -\!\bigl(\sin\delta\cos\lambda\sin 2\Phi + \sin 2\delta\sin\lambda\sin^2\!\Phi\bigr) \\
M_{yy} &= \phantom{-}\sin\delta\cos\lambda\sin 2\Phi - \sin 2\delta\sin\lambda\cos^2\!\Phi \\
M_{zz} &= \phantom{-}\sin 2\delta\sin\lambda \\
M_{xy} &= \phantom{-}\sin\delta\cos\lambda\cos 2\Phi + \tfrac{1}{2}\sin 2\delta\sin\lambda\sin 2\Phi \\
M_{xz} &= -\!\bigl(\cos\delta\cos\lambda\cos\Phi + \cos 2\delta\sin\lambda\sin\Phi\bigr) \\
M_{yz} &= -\!\bigl(\cos\delta\cos\lambda\sin\Phi - \cos 2\delta\sin\lambda\cos\Phi\bigr)
\end{aligned}
```

---

## Receivers

### Geophones

Point receivers record the three particle velocity components $v_x$, $v_y$, $v_z$ at each time step.
Due to the staggered arrangement, each component is sampled at its natural staggered grid location
(offset by half a cell from the nominal receiver position).

### DAS

Distributed acoustic sensing (DAS) receivers record **axial strain** along coordinate-aligned profiles.
Within the staggered-grid framework, strain can be reconstructed from the co-located normal stress
components via the compliance relation:

```math
\begin{bmatrix} \varepsilon_{xx} \\ \varepsilon_{yy} \\ \varepsilon_{zz} \end{bmatrix}
= \begin{bmatrix} C_{11} & C_{12} & C_{13} \\ C_{12} & C_{22} & C_{23} \\ C_{13} & C_{23} & C_{33} \end{bmatrix}^{-1}
\begin{bmatrix} \sigma_{xx} \\ \sigma_{yy} \\ \sigma_{zz} \end{bmatrix}
```

The axial strain component along the fiber orientation is extracted from the resulting vector.
Gauge-length integration or conversion to strain rate can be performed in post-processing.
Only fibers aligned with the $x$-, $y$-, or $z$-axis are supported.

---

## References

- Virieux, J. (1986). P-SV wave propagation in heterogeneous media: Velocity-stress finite-difference method. *Geophysics*, 51(4), 889–901.
- Levander, A. R. (1988). Fourth-order finite-difference P-SV seismograms. *Geophysics*, 53(11), 1425–1436.
- Komatitsch, D., & Martin, R. (2007). An unsplit convolutional perfectly matched layer improved at grazing incidence for the seismic wave equation. *Geophysics*, 72(5), SM155–SM167.
- Shi, P., Angus, D., Nowacki, A., Yuan, S., & Wang, Y. (2018). Microseismic full waveform modeling in anisotropic media with moment tensor implementation. *Surveys in Geophysics*, 39, 567–611.
- Bécache, E., Fauqueux, S., & Joly, P. (2003). Stability of perfectly matched layers, group velocities and anisotropic waves. *Journal of Computational Physics*, 188(2), 399–433.
