var documenterSearchIndex = {"docs":
[{"location":"userguide/config/#Configurations","page":"Configurations","title":"Configurations","text":"","category":"section"},{"location":"userguide/config/","page":"Configurations","title":"Configurations","text":"Empty template configuration.yaml files can be generated with:","category":"page"},{"location":"userguide/config/","page":"Configurations","title":"Configurations","text":"ElasticFDSG.dim2.configtemplate() for 2D simulations \nElasticFDSG.dim3.configtemplate() for 3D simulations.","category":"page"},{"location":"userguide/config/","page":"Configurations","title":"Configurations","text":"The contents of these files must then be fully completed by the user either manually or with self-written scripts.","category":"page"},{"location":"userguide/config/","page":"Configurations","title":"Configurations","text":"note: Note\nThe configuration file reader yet dont perform any type or spell checks. User should fill the form like indicated below. For instance, settings.spatial_derivative_order: 10 will work, while \"10\" would throw an error.","category":"page"},{"location":"userguide/config/#2D","page":"Configurations","title":"2D","text":"","category":"section"},{"location":"userguide/config/","page":"Configurations","title":"Configurations","text":"using ElasticFDSG \n\nconfiguration_file_path = joinpath(@__DIR__, \"my_2D_configurations.yaml\") # add your actual path here \nElasticFDSG.dim2.configtemplate(configuration_file_path)\n","category":"page"},{"location":"userguide/config/","page":"Configurations","title":"Configurations","text":"This creates:","category":"page"},{"location":"userguide/config/","page":"Configurations","title":"Configurations","text":"# This is a template configuration.yaml file for a 2D ElasticFDSG simulation.\n# Any other configuration file can be prepared in the same manner.\n# All keywords are case sensitive, and must be named as shown below.\n# Users must fill the whole template before running a simulation. \n# The velocity model is prepared in another file.\n\nsettings:\n    device:                          # cpu / gpu-cuda / gpu-metal \n    precision:                       # Float64 / Float32\n    spatial_derivative_order:        # 1-10\n    show_summary_in_console:         # true / false  \n    show_progress_in_console:        # true / false\n    save_results:                    # true / false  \n    output:\n        destination_folder: \n        file_name:  \n        format: h5                   # h5 (only supported format yet) \n\ntime:\n    start:                # [sec]\n    end:                  # [sec]\n    timestep:             # [sec]\n\nsource:\n    dominant_frequency:        # [Hz]\n    wavelet_type:              # ricker / gauss1d \n    wavelet_center:            # (should be ≥ 1/fdom) [sec]\n    amplitude:                \n    force_angle:               # [°]            \n    point_source:              # true / false  (if true, double_couple should be false)\n    double_couple:             # true / false  (if true, point_source should be false)\n    location:\n        x:                     # [m]\n        y:                     # [m]\n\nboundaries: # absorbing / else (reflecting)\n    left:         # [x-start] \n    right:        # [x-end]    \n    top:          # [y-start]  \n    bottom:       # [y-end]  \n\npml:\n    nlayer:                       # 5-15 recommended\n    reflection_coefficient:       # 1e-10 recommended\n\nreceivers:\n    geophones:\n        # - { x: 3000, y: 3000 }        # [m] Example geophone\n\n    das:\n        x_aligned:\n        #- { x_range: \"0:25:10000\", y: 4000 }  # [m] Example das\n      \n        y_aligned:\n        #- { x: 1000, y_range: \"1000:25:9000\"} # [m] Example das\n        \n        \nsnapshots:\n    times: \n        # [0.25, 0.5, 0.75, 1.0] # [sec] Example\n    fields: \n        # [vx, vy, sxx, sxy, syy] # Example \n","category":"page"},{"location":"userguide/config/#3D","page":"Configurations","title":"3D","text":"","category":"section"},{"location":"userguide/config/","page":"Configurations","title":"Configurations","text":"using ElasticFDSG \n\nconfiguration_file_path = joinpath(@__DIR__, \"my_3D_configurations.yaml\") # add your actual path here \nElasticFDSG.dim3.configtemplate(configuration_file_path)\n","category":"page"},{"location":"userguide/config/","page":"Configurations","title":"Configurations","text":"This creates:","category":"page"},{"location":"userguide/config/","page":"Configurations","title":"Configurations","text":"# This is a template configuration.yaml file for a 3D ElasticFDSG simulation.\n# Any other configuration file can be prepared in the same manner.\n# All keywords are case sensitive, and must be named as shown below.\n# Users must fill the whole template before running a simulation.  \n# The velocity model is prepared in another file.\n\nsettings:\n    device:                          # cpu / gpu-cuda / gpu-metal \n    precision:                       # Float64 / Float32\n    spatial_derivative_order:        # 1-10\n    show_summary_in_console:         # true / false  \n    show_progress_in_console:        # true / false\n    save_results:                    # true / false  \n    output:\n        destination_folder: \n        file_name:  \n        format: h5                   # h5 (only supported format yet) \n\ntime:\n    start:                # [sec]\n    end:                  # [sec]\n    timestep:             # [sec]\n\nsource:\n    dominant_frequency:          # [Hz]\n    wavelet_type:                # ricker / gauss1d \n    wavelet_center:              # (should be ≥ 1/fdom) [sec]\n    amplitude:                        \n    point_source:          \n        use:                     # true / false  (if true, double_couple.use should be false)\n        phi:                     # [°] Elevation\n        theta:                   # [°] Azimuth \n    double_couple: \n        use:                     # true / false  (if true, point_source.use should be false)\n        strike:                  # [°] ∈ [0, 360]\n        dip:                     # [°] ∈ [0, 90]\n        rake:                    # [°] ∈ [0, 360]\n        anisotropic_moment_tensor:      # true / false\n    location:\n        x:              # [m]\n        y:              # [m]\n        z:              # [m]\n\nboundaries: # absorbing / else (reflecting)\n    xstart:       \n    xend:        \n    ystart:        \n    yend:          \n    zstart:         \n    zend:           \n\npml:\n    nlayer:                      # 5-15 recommended\n    reflection_coefficient:      # 1e-10 recommended\n\nreceivers:\n    geophones: \n        #- { x: 100, y: 200, z: 300 }    #[m] Example geophone\n\n    das:\n        x_aligned: \n        #- { x_range: \"0:2.5:500\", y: 100, z: 100 }  #[m] Example das \n       \n        y_aligned: \n        #- { x: 200, y_range: \"100:2.5:400\", z: 100 } #[m] Example das \n\n        z_aligned:\n        #- { x: 300, y: 250, z_range: \"0:2.5:400\" } #[m] Example das \n                \n                \n    snapshots: # XY-, XZ-, YZ-planes only\n        times: \n            # [0.25, 0.5, 0.75, 1.0] # [sec] Example\n        fields: \n            # [vx, vy, vz, sxx, sxy, sxz, syy, syz, szz] # Example \n        origins:\n            # { x: 250, y: 250, z: 250 } # [m] Coordinates of snapshot origin. Example: XY-plane snapshots will be at z=250m\n\n","category":"page"},{"location":"method/#Theory","page":"Method","title":"Theory","text":"","category":"section"},{"location":"method/","page":"Method","title":"Method","text":"This section provides a brief overview of the theoretical basis of the application. For more comprehensive descriptions, readers are referred to literature listed below.","category":"page"},{"location":"method/","page":"Method","title":"Method","text":"Elastic wave propagation in a 3D anisotropic linear elastic material can be described by the following set of first order partial differential equations:","category":"page"},{"location":"method/","page":"Method","title":"Method","text":"beginaligned\npartial_t sigma_xx =  c_11 partial_x v_x + c_12 partial_y v_y + c_13 partial_z v_z   \npartial_t sigma_yy =  c_12 partial_x v_x + c_22 partial_y v_y + c_23 partial_z v_z  \npartial_t sigma_zz =  c_13 partial_x v_x + c_23 partial_y v_y + c_33 partial_z v_z  \npartial_t sigma_xy =  c_66 (partial_x v_y + partial_y v_x)  \npartial_t sigma_xz =  c_55 (partial_x v_z + partial_z v_x)  \npartial_t sigma_yz =  c_44 (partial_y v_z + partial_z v_y) tag1 \n\nendaligned","category":"page"},{"location":"method/","page":"Method","title":"Method","text":"beginaligned\nrho partial_t v_x = partial_x sigma_xx + partial_y sigma_xy + partial_z sigma_xz + f_x \nrho partial_t v_y = partial_x sigma_xy + partial_y sigma_yy + partial_z sigma_yz + f_y \nrho partial_t v_z = partial_x sigma_xz + partial_y sigma_xz + partial_z sigma_zz + f_z tag2 \nendaligned","category":"page"},{"location":"method/","page":"Method","title":"Method","text":"In these equations:","category":"page"},{"location":"method/","page":"Method","title":"Method","text":"are sigma_xx sigma_yy sigma_zz the normal stress components.\nare sigma_xy sigma_xz sigma_yz the shear stress components.\nare c_11 c_12 c_13 c_22 c_23 c_33 c_44 c_55 c_66 the 9 elastic constants needed to describe an orthorhombic material.\nis rho is the density of the material.\nare v_x v_y v_z the particle velocities in the x, y, and z directions, respectively.\nare f_x f_y f_z the external body forces acting in the x, y, and z directions, respectively.","category":"page"},{"location":"method/","page":"Method","title":"Method","text":"Using the Tsvankin notation, isotropic, vertical transversal isotropic (VTI) and ortorhombic materials (ORT) can be characterized by two vertical velocities and 7 dimensionless paramters epsilon_1 epsilon_2 gamma_1 gamma_2 delta_1 delta_2 delta_3 :","category":"page"},{"location":"method/","page":"Method","title":"Method","text":"beginaligned\nc_33 = v_p0^2 cdot rho \nc_55 = v_s0^2 cdot rho \nc_11 = (2 epsilon_2 + 1) cdot c_33 \nc_22 = c_33 cdot (2 epsilon_1 + 1) \nc_66 = c_55 cdot (2 gamma_1 + 1) \nc_44 = fracc_661 + gamma_2 \nc_13 = sqrt2 c_33 cdot (c_33 - c_55) cdot delta_2 + (c_33 - c_55)^2 - c_55 \nc_23 = sqrt2 c_33 cdot (c_33 - c_44) cdot delta_1 + (c_33 - c_44)^2 - c_44 \nc_12 = sqrt2 c_11 cdot (c_11 - c_66) cdot delta_3 + (c_11 - c_66)^2 - c_66\nendaligned","category":"page"},{"location":"method/#Numerical-Scheme","page":"Method","title":"Numerical Scheme","text":"","category":"section"},{"location":"method/","page":"Method","title":"Method","text":"The set of equations (1) and (2) can be solved using a finite-difference staggered-grid scheme. In a staggered grid, field quantities are not co-located on same grid points but are distributed across predefined grid cells. The primary advantage of a staggered grid is an enhanced accuracy of spatial central difference operators, as well as staggered temporal finite difference operators, which also improves the accuracy of the time marching.","category":"page"},{"location":"method/","page":"Method","title":"Method","text":"However, the distribution of field quantities within a grid cell also has drawbacks, especially during the processing of the results. For example, users working with geophone data must consider that the velocity components v_x v_y v_z are not measured at the same location.","category":"page"},{"location":"method/","page":"Method","title":"Method","text":"The staggered grid scheme used in the application is illustrated below:","category":"page"},{"location":"method/","page":"Method","title":"Method","text":"      ↑ y                           │      y              △             ■ vₓ  i,j,k\n        │         i+½,j+½           │  z↑ /              /              ◆ vᵧ  i+½,j+½,k\n  i,j+½ │         │vᵧ               │   │/              /               ● v𝑧  i+½,j,k+½ \n   σₓᵧ ─◇─────────◆───              │   ○ ............ ●\n        │         │                 │   ╎              ╎                □ σₓₓ,σᵧᵧ,σ𝑧𝑧 i+½,j,k\n        │         │                 │   ╎              ╎                ○ σₓ𝑧 i,j,k+½\n        │         │                 │   ╎   ◇          ╎  ◆             ◇ σₓᵧ i,j+½,k\n        ■─────────□──── ──→ x       │   ╎  /           ╎ /              △ σᵧ𝑧 i+½,j+½,k+½\n        vₓ      σₓₓ,σᵧᵧ             │   ╎ /            ╎/  \n       i,j       i+½,j              │   ■ ............ □ ─→ x \n                                    │","category":"page"},{"location":"method/","page":"Method","title":"Method","text":"Using this scheme, the discrete form of (1) and (2) are given by:","category":"page"},{"location":"method/","page":"Method","title":"Method","text":"beginaligned\n\nsigma_xx  (i+frac12jk)^t_n+frac12 = sigma_xx  (i+frac12jk)^t_n-frac12 Delta t (c_11 mathcalD_x v_x + c_12 mathcalD_y v_y + c_13 mathcalD_z v_z ) bigg_(i+frac12jk)^t_n  \n\nsigma_yy  (i+frac12jk)^t_n+frac12 = sigma_yy  (i+frac12jk)^t_n-frac12 Delta t (c_12 mathcalD_x v_x + c_22 mathcalD_y v_y + c_23 mathcalD_z v_z ) bigg_(i+frac12jk)^t_n  \n\nsigma_zz  (i+frac12jk)^t_n+frac12 = sigma_zz  (i+frac12jk)^t_n-frac12 Delta t (c_13 mathcalD_x v_x + c_23 mathcalD_y v_y + c_33 mathcalD_z v_z ) bigg_(i+frac12jk)^t_n  \n\nsigma_xy  (ij+frac12k)^t_n+frac12 = sigma_xy  (ij+frac12k)^t_n-frac12 Delta t  c_66 ( mathcalD_x v_y +  mathcalD_y v_x)bigg_(ij+frac12k)^t_n \n\nsigma_xz  (ijk+frac12)^t_n+frac12 = sigma_xz  (ijk+frac12)^t_n-frac12 Delta t  c_55 ( mathcalD_x v_z + mathcalD_z v_x)bigg_(ijk+frac12)^t_n \n\nsigma_yz  (i+frac12j+frac12k+frac12)^t_n+frac12 = sigma_yz  (i+frac12j+frac12k+frac12)^t_n-frac12 Delta t  c_44 ( mathcalD_y v_z +  mathcalD_y v_z)bigg_(i+frac12j+frac12k+frac12)^t_n \n            \n            \nendaligned","category":"page"},{"location":"method/","page":"Method","title":"Method","text":"beginaligned\nv_x  (ijk)^t_n+1 = v_x  (ijk)^t_n fracDelta trho (mathcalD_x sigma_xx + mathcalD_y sigma_xy + mathcalD_z sigma_xz + f_x) bigg_(ijk)^t_n+frac12 \nv_y  (i+frac12j+frac12k)^t_n+1 = v_y  (i+frac12j+frac12k)^t_n fracDelta trho (mathcalD_x sigma_xy + mathcalD_y sigma_yy + mathcalD_z sigma_yz + f_y) bigg_(i+frac12j+frac12k)^t_n+frac12 \nv_z  (i+frac12jk+frac12)^t_n+1 = v_z  (i+frac12jk+frac12)^t_n fracDelta trho (mathcalD_x sigma_xz + mathcalD_y sigma_xz + mathcalD_z sigma_zz + f_z) bigg_(i+frac12jk+frac12)^t_n+frac12 tag3\n\nendaligned","category":"page"},{"location":"method/","page":"Method","title":"Method","text":"Here (ijk) represent the grid points in the x,y,z-directions, respectively;  t_n denotes the n-th time step;  and Delta t the time increment.","category":"page"},{"location":"method/","page":"Method","title":"Method","text":"The differential operators mathcalD are given by central difference approximations:","category":"page"},{"location":"method/","page":"Method","title":"Method","text":"beginaligned\n\nmathcalD_x f = sum_n=1^N fracc_nDelta x ( f(x_i+n y_j z_k) - f(x_i- n y_j z_k)) \nmathcalD_y f = sum_n=1^N fracc_nDelta y ( f(x_i y_j+n z_k) - f(x_i y_j-n z_k)) \nmathcalD_z f = sum_n=1^N fracc_nDelta z ( f(x_i y_j z_k+n) - f(x_i y_j z_k-n))\n        \nendaligned","category":"page"},{"location":"method/","page":"Method","title":"Method","text":"with order N and coefficients c_n.","category":"page"},{"location":"method/","page":"Method","title":"Method","text":"To compute a spatial derivative of order N at a specific location, at least N neighboring nodes are required on both sides. As a result, for edge nodes no spatial derivatives can be calculated. These edge nodes (ghost nodes) are effective model boundaries and should lie outside the physical (user-defined) domain. The application automatically extends the user-defined domain with N ghost node layers.","category":"page"},{"location":"method/","page":"Method","title":"Method","text":"In staggered grids, it is often beneficial to assign certain elastic properties to specific points within the grid cell. However, requiring users to define such grids can become difficult to manage, particularly for complex media. To simplify this, we define all elastic properties at full integer grid points (e.g., on v_x) and effective properties are calculated by interpolating values from neighboring nodes.","category":"page"},{"location":"method/#Sources","page":"Method","title":"Sources","text":"","category":"section"},{"location":"method/","page":"Method","title":"Method","text":"Earthquake simulations require the excitation of point or double-couple forces using the body force term in the equation of motion.","category":"page"},{"location":"method/","page":"Method","title":"Method","text":"Point Sources:","category":"page"},{"location":"method/","page":"Method","title":"Method","text":"Point sources can directly be applied to the velocity components at the desired source location (s_xs_ys_y), using elevation and azimuth angles combined with a source-time function that contains a wavelet. However, due to the positions of velocity components in the grid cell, the excitation occurs at slightly different locations.","category":"page"},{"location":"method/","page":"Method","title":"Method","text":"Double couple sources","category":"page"},{"location":"method/","page":"Method","title":"Method","text":"Double-couple sources involve applying force pairs (pq) with strength M_pq to the velocity field at a desired source location. In the used staggered grid scheme, double-couple sources are centered around normal stresses grid points, i.e., (i+frac12jk). As a result, a defined double couple source location is centered shifted by half a grid point in the x-direction!  Using moment tensors to represent the equivalent distribution of force pairs, at least 30 force contributions need to be applied to the surounding velocity field. The 10 force components contributing to f_x are:","category":"page"},{"location":"method/","page":"Method","title":"Method","text":"beginaligned\nf_x  ijk = - f_x  i+1jk = fracM_xx(t)Delta x^2 Delta y Delta z \n\nf_x  ij+1k = - f_x  ij-1k = fracM_xy(t)4 Delta x Delta y^2 Delta z \nf_x  i+1j+1k = - f_x  i+1j-1k = fracM_xy(t)4 Delta x  Delta y^2  Delta z \n\nf_x  ijk+1 = - f_x  ijk-1 = fracM_xz(t)4 Delta x Delta y Delta z^2 \nf_x  i+1jk+1 = - f_x i+1jk-1 = fracM_xz(t)4 Delta x Delta y Delta z^2 \n\nendaligned","category":"page"},{"location":"method/","page":"Method","title":"Method","text":"The same logic applies to the body force components contributing to f_y and f_z.","category":"page"},{"location":"method/","page":"Method","title":"Method","text":"The moment tensor components are computed by user defined dip, strike and rake values for the following coordinate system:","category":"page"},{"location":"method/","page":"Method","title":"Method","text":"     x                  \n   /                     \n  /                        δ = Dip, λ = Rake, Φ = Strike                \n  ---- y                   Coordinate System: x -> north,            \n |                                            y -> east,\n |                                            z -> positive downward\n z","category":"page"},{"location":"method/","page":"Method","title":"Method","text":"beginaligned\nM_xx = -(sin(delta)  cos(lambda)  sin(2Phi) + sin(2delta)  sin(lambda)  sin(Phi)^2) \nM_xy = sin(delta)  cos(lambda)  cos(2Phi) + frac12 sin(2delta)  sin(lambda)  sin(2Phi) \nM_xz = -(cos(delta)  cos(lambda)  cos(Phi) + cos(2delta)  sin(lambda)  sin(Phi)) \nM_yy = sin(delta)  cos(lambda)  sin(2Phi) - sin(2delta)  sin(lambda)  cos(Phi)^2 \nM_yz = -(cos(delta)  cos(lambda)  sin(Phi) - cos(2delta)  sin(lambda)  cos(Phi)) \nM_zz = sin(2delta)  sin(lambda)\nendaligned","category":"page"},{"location":"method/#C-PML","page":"Method","title":"C-PML","text":"","category":"section"},{"location":"method/","page":"Method","title":"Method","text":"Earthquake simulations often require modeling wave propagation in unbounded media. Perfectly Matched Layer (PML) is a very effective method to prevent artificial reflections at model boundaries. The underlying idea is to manipulate the wave equation to obtain exponentially decaying plane wave solutions for complex arguments. Assuming that the computational domain (real arguments) is surrounded by a complex region (the PML region), the amplitudes of incident waves decay and cause negligible reflections at the model boundaries.","category":"page"},{"location":"method/","page":"Method","title":"Method","text":"Komatitsch & Martin (2007) introduced a memory efficient convolution-based unsplit PML formulation (C-PML).  This approach requires storing one additional memory variable for each spatial derivative, but only in the PML region. For the PML region, spatial derivative operators are replaced by:","category":"page"},{"location":"method/","page":"Method","title":"Method","text":"beginaligned\nmathcalD_hatx = fracmathcalD_xkappa_x + Psi_x tag4\nendaligned","category":"page"},{"location":"method/","page":"Method","title":"Method","text":"Here, Psi_x represents the memory variable associated with the field from which the derivative is taken. The memory variable Psi_x is updated at each time step by:","category":"page"},{"location":"method/","page":"Method","title":"Method","text":"beginaligned\nPsi_x^t_n = b_x Psi_x^t_n-1 + a_x (mathcalD_x)^t_n + frac12 tag5\nendaligned","category":"page"},{"location":"method/","page":"Method","title":"Method","text":"while a_x, b_x kappa_x are precomputed PML related parameter. ","category":"page"},{"location":"method/","page":"Method","title":"Method","text":"Accordingly, differential operators mathcalD_y and mathcalD_z are replaced by mathcalD_haty and mathcalD_hatz in the PML-region.","category":"page"},{"location":"method/#References","page":"Method","title":"References","text":"","category":"section"},{"location":"method/","page":"Method","title":"Method","text":"Komatitsch, D., & Martin, R. (2007). An unsplit convolutional perfectly matched layer improved at grazing incidence for the seismic wave equation. Geophysics, 72(5), SM155-SM167.","category":"page"},{"location":"method/","page":"Method","title":"Method","text":"Moczo, P., Kristek, J., & Gális, M. (2014). The finite-difference modelling of earthquake motions: Waves and ruptures. Cambridge University Press.","category":"page"},{"location":"method/","page":"Method","title":"Method","text":"Tsvankin, I. (1997). Anisotropic parameters and P-wave velocity for orthorhombic media. Geophysics, 62(4), 1292-1309.","category":"page"},{"location":"method/","page":"Method","title":"Method","text":"Virieux, J. (1984). SH-wave propagation in heterogeneous media: Velocity-stress finite-difference method. Geophysics, 49(11), 1933-1942. ","category":"page"},{"location":"userguide/velmod/#Velocity-Models","page":"Velocity Models","title":"Velocity Models","text":"","category":"section"},{"location":"userguide/velmod/","page":"Velocity Models","title":"Velocity Models","text":"Velocity models can be prepared in .jld2 (Julia) or .npy/.npz (Python) formats.  This gives users the freedom to design their velocity models in the tool of their choice and adapt them to the required format.  The following explains the required structure for 2D and 3D models and shows how to create simple velocity models in Julia.","category":"page"},{"location":"userguide/velmod/","page":"Velocity Models","title":"Velocity Models","text":"note: Note\nTo simplify the creation of velocity models, all elastic parameters are defined on full-integer grid points in the staggered grid.","category":"page"},{"location":"userguide/velmod/#2D","page":"Velocity Models","title":"2D","text":"","category":"section"},{"location":"userguide/velmod/","page":"Velocity Models","title":"Velocity Models","text":"The 2D solver expects an (7, ny, nx) array with:","category":"page"},{"location":"userguide/velmod/","page":"Velocity Models","title":"Velocity Models","text":"1: X - 2D meshgrid coordinates\n2: Y - 2D meshgrid coordinates \n3: P-wave velocities [m/s]\n4: S-wave velocities [m/s]\n5: Densities [kg/m^3]\n6: 2D Thomsen Parameter epsilon \n7: 2D Thomsen Parameter delta ","category":"page"},{"location":"userguide/velmod/","page":"Velocity Models","title":"Velocity Models","text":"For an isotropic medium, set 6 & 7 zero.","category":"page"},{"location":"userguide/velmod/","page":"Velocity Models","title":"Velocity Models","text":"The following script shows how to create a simple 2D velocity model.","category":"page"},{"location":"userguide/velmod/","page":"Velocity Models","title":"Velocity Models","text":"\nusing JLD2\n\n# spatial extends\nx_start = 0\nx_end = 10000\n\ny_start = 0 \ny_end = 10000\n\ndx = 10 # cell size x-direction \ndy = 10 # cell size y-direction \n\nxcoords = x_start:dx:x_end # x-coordinates\nycoords = y_start:dy:y_end # y-coordinates\n\nnx = length(xcoords) # number of grid points x-direction\nny = length(ycoords) # number of grid points y-direction\n\ndim = (ny, nx) # model dimensions \n\n# 2D meshgrid\nX = repeat(xcoords', ny, 1)\nY = repeat(ycoords,  1, nx)\n\nvp = zeros(dim);   # P-wave velocity\nvs = zeros(dim);   # S-wave velocity\nrho = zeros(dim);  # Density\neps0 = zeros(dim); # 2D Thomson parameter epsilon\ndel0 = zeros(dim); # 2D Thomson parameter delta\n\n# fill arrays with values\nvp[:,:]  .= 5000;  \nvs[:,:]  .= 2500;\nrho[:,:] .= 2800;\n\n# thomson parameter for vti media (zero for isotropic medium)\neps0[:,:] .= 0.2\ndel0[:,:] .= -0.15\n\n# velocity model array\nveldim = (7, ny, nx)\nvelmod = zeros(veldim)\n# fill with elastic properties\nvelmod[1,:,:] .= X\nvelmod[2,:,:] .= Y\nvelmod[3,:,:] .= vp\nvelmod[4,:,:] .= vs\nvelmod[5,:,:] .= rho\nvelmod[6,:,:] .= eps0\nvelmod[7,:,:] .= del0\n\n# save the velocity model\npath = joinpath(@__DIR__,\"velmod.jld2\") # add your actual path here \njldsave(path; velmod)","category":"page"},{"location":"userguide/velmod/#3D","page":"Velocity Models","title":"3D","text":"","category":"section"},{"location":"userguide/velmod/","page":"Velocity Models","title":"Velocity Models","text":"The 3D solver expects an (13, nx, ny, nz) array with:","category":"page"},{"location":"userguide/velmod/","page":"Velocity Models","title":"Velocity Models","text":"1: X - 3D meshgrid coordinates \n2: Y - 3D meshgrid coordinates \n3: Z - 3D meshgrid coordinates \n4: P-wave velocities [m/s] \n5: S-wave velocities [m/s] \n6: Densities [kg/m^3] ","category":"page"},{"location":"userguide/velmod/","page":"Velocity Models","title":"Velocity Models","text":"Tsvankin Parameter ","category":"page"},{"location":"userguide/velmod/","page":"Velocity Models","title":"Velocity Models","text":"7: varepsilon_1 \n8: varepsilon_2 \n9: gamma_1 \n10: gamma_2 \n11: delta_1 \n12: delta_2\n13: delta_3 ","category":"page"},{"location":"userguide/velmod/","page":"Velocity Models","title":"Velocity Models","text":"For an isotropic medium, set 7-13 to zero.","category":"page"},{"location":"userguide/velmod/","page":"Velocity Models","title":"Velocity Models","text":"For a VTI medium, the relationship between Tsvankin and Thomsen Parameters can be used:","category":"page"},{"location":"userguide/velmod/","page":"Velocity Models","title":"Velocity Models","text":"Tsvankin Parameter <=> Thomsen Parameter\nvarepsilon_1 = varepsilon_2 <=> varepsilon\ngamma_1 = gamma_2 <=> gamma\ndelta_1 = delta_2, delta_3 = 0 <=> delta","category":"page"},{"location":"userguide/velmod/","page":"Velocity Models","title":"Velocity Models","text":"The following script shows how to create a simple 3D velocity model.","category":"page"},{"location":"userguide/velmod/","page":"Velocity Models","title":"Velocity Models","text":"\nusing JLD2\n\n# spatial extends\nx_start = 0\nx_end = 500\n\ny_start = 0 \ny_end = 500\n\nz_start = 0\nz_end = 500\n\ndx = 5 # cell size x-direction \ndy = 5 # cell size y-direction \ndz = 5 # cell size z-direction \n\nxcoords = x_start:dx:x_end # x-coordinates\nycoords = y_start:dy:y_end # y-coordinates\nzcoords = z_start:dz:z_end # z-coordinates\n\nnx = length(xcoords) # number of grid points x-direction\nny = length(ycoords) # number of grid points y-direction\nnz = length(zcoords) # number of grid points z-direction\n\n# 3D Meshgrid\nX = getindex.(Iterators.product(xcoords, ycoords, zcoords), 1)\nY = getindex.(Iterators.product(xcoords, ycoords, zcoords), 2)\nZ = getindex.(Iterators.product(xcoords, ycoords, zcoords), 3)\n\ndim = (nx, ny, nz) # model dimension\n\nvp = zeros(dim);    # P-wave velocity\nvs = zeros(dim);    # S-wave velocity\nrho = zeros(dim);   # Density velocity\neps1 = zeros(dim);  # Tsvankin parameter epsilon 1 \neps2 = zeros(dim);  # Tsvankin parameter epsilon 2\ngam1 = zeros(dim);  # Tsvankin parameter gamma 1 \ngam2 = zeros(dim);  # Tsvankin parameter gamma 2 \ndel1 = zeros(dim);  # Tsvankin parameter delta 1 \ndel2 = zeros(dim);  # Tsvankin parameter delta 2 \ndel3 = zeros(dim);  # Tsvankin parameter delta 3 \n\n# fill arrays with values\nvp[:,:,:] .= 5000;\nvs[:,:,:]  .= 2500;\nrho[:,:,:] .= 2800;\n\n# Tsvankin Parameter\neps1[:,:,:]  .= 0.05\neps2[:,:,:]  .= 0.1\ngam1[:,:,:]  .= 0.1\ngam2[:,:,:]  .= -0.05\ndel1[:,:,:]  .= 0.05\ndel2[:,:,:]  .= 0.025\ndel3[:,:,:]  .= -0.1\n\n# velocity model array\nveldim = (13, nx, ny, nz)\nvelmod = zeros(veldim);\n# fill with elastic properties\nvelmod[1,:,:,:] .= X\nvelmod[2,:,:,:] .= Y\nvelmod[3,:,:,:] .= Z\n\nvelmod[4,:,:,:] .= vp\nvelmod[5,:,:,:] .= vs\nvelmod[6,:,:,:] .= rho\n\nvelmod[7,:,:,:] .= eps1\nvelmod[8,:,:,:] .= eps2\nvelmod[9,:,:,:] .= gam1\nvelmod[10,:,:,:] .= gam2\nvelmod[11,:,:,:] .= del1\nvelmod[12,:,:,:] .= del2\nvelmod[13,:,:,:] .= del3\n\n# save velocity model in jld2 file\nusing JLD2\npath = joinpath(@__DIR__,\"velmod.jld2\") # add your actual path here \njldsave(path; velmod)","category":"page"},{"location":"userguide/intro/#General-Usage","page":"General Usage","title":"General Usage","text":"","category":"section"},{"location":"userguide/intro/","page":"General Usage","title":"General Usage","text":"The core functionalities of the application are:","category":"page"},{"location":"userguide/intro/","page":"General Usage","title":"General Usage","text":"ElasticFDSG.dim2.runsim() for 2D simulations \nElasticFDSG.dim3.runsim() for 3D simulations.","category":"page"},{"location":"userguide/intro/","page":"General Usage","title":"General Usage","text":"To run a simulation, the following inputs are required:","category":"page"},{"location":"userguide/intro/","page":"General Usage","title":"General Usage","text":"(1) a velocity model stored in a julia (.jld2) or numpy array (.npy/npz)","category":"page"},{"location":"userguide/intro/","page":"General Usage","title":"General Usage","text":"(2) a configuration.yaml file specifying simulation parameters.","category":"page"},{"location":"userguide/intro/","page":"General Usage","title":"General Usage","text":"Once the inputs are ready, a simulation can be initiated like:","category":"page"},{"location":"userguide/intro/","page":"General Usage","title":"General Usage","text":"\nusing ElasticFDSG\n\n# paths\npath_to_configfile = joinpath(@__DIR__, \"config.yaml\")\npath_to_velmodfile = joinpath(@__DIR__, \"velmod.jld2\")\n\n# run a 2D simulation\nElasticFDSG.dim2.runsim(path_to_configfile, path_to_velmodfile)\n\n# run a 3D simulation\nElasticFDSG.dim3.runsim(path_to_configfile, path_to_velmodfile)\n","category":"page"},{"location":"userguide/intro/","page":"General Usage","title":"General Usage","text":"After the calculations are completed, the results will be saved as a .h5 file at the specified location. These results can then be processed using any tool of choice that supports HDF5. Below is an example of how to access the content of the results file with Julia:","category":"page"},{"location":"userguide/intro/","page":"General Usage","title":"General Usage","text":"using HDF5\n\n# extract content from result file \nfunction extract_hdf5_content(file_path::String)\n    file = h5open(file_path, \"r\")\n    result = Dict()\n    \n    function process_group(group, prefix=\"\")\n        for name in keys(group)\n            path = joinpath(prefix, name)\n            obj = group[name]\n            if obj isa HDF5.Group\n                process_group(obj, path)\n            elseif obj isa HDF5.Dataset\n                result[path] = read(obj)\n            else\n                ;\n            end\n        end\n    end\n    \n    process_group(file)\n    close(file)\n    return result\nend;\n\nfile_path = joinpath(@__DIR__, \"results.h5\"); # add the actual path here\ncontent = extract_hdf5_content(file_path);\n\nprintln(\"Content keys\")\nfor (key, value) in content\n    println(\"$key\")\nend","category":"page"},{"location":"userguide/intro/","page":"General Usage","title":"General Usage","text":"In the next sections, the process for creating required velocity models and configuration files will be demonstrated.","category":"page"},{"location":"reference/#API-Reference","page":"API Reference","title":"API Reference","text":"","category":"section"},{"location":"reference/","page":"API Reference","title":"API Reference","text":"Pages = [\"reference.md\"]","category":"page"},{"location":"reference/","page":"API Reference","title":"API Reference","text":"Pages = [\"reference.md\"]","category":"page"},{"location":"reference/","page":"API Reference","title":"API Reference","text":"Modules = [ElasticFDSG.dim2, ElasticFDSG.dim3]","category":"page"},{"location":"reference/#ElasticFDSG.dim2.configtemplate-Tuple{String}","page":"API Reference","title":"ElasticFDSG.dim2.configtemplate","text":"Elastic_FDSG.dim2.configtemplate(path)\n\nCreates and empty template configuration file for a 2D ElasticFDSG simulation.  The user must fill the template afterwards. \n\nArguments\n\npath::String: Path where the template is saved.\n\nReturns\n\nNothing: The function saves the template as specified in path.\n\n\n\n\n\n","category":"method"},{"location":"reference/#ElasticFDSG.dim2.runsim-Tuple{String, String}","page":"API Reference","title":"ElasticFDSG.dim2.runsim","text":"Elastic_FDSG.dim2.runsim(CONFIGPATH, VELMODPATH)\n\nRun the 2D elastic forward simulation using the specified configuration and velocity model.\n\nArguments\n\nCONFIGPATH::String: Path to the configuration.yaml file containing simulation settings.\nVELMODPATH::String: Path to the velocity model file. Supported formats include .jld2, .npy, and .npz.\n\nReturns\n\nNothing: The function runs the simulation and saves the results as specified in the configuration file.\n\n\n\n\n\n","category":"method"},{"location":"reference/#ElasticFDSG.dim3.configtemplate-Tuple{String}","page":"API Reference","title":"ElasticFDSG.dim3.configtemplate","text":"Elastic_FDSG.dim3.configtemplate(path)\n\nCreates and empty template configuration file for a 3D ElasticFDSG simulation.  The user must fill the template afterwards. \n\nArguments\n\npath::String: Path where the template is saved.\n\nReturns\n\nNothing: The function saves the template as specified in path.\n\n\n\n\n\n","category":"method"},{"location":"reference/#ElasticFDSG.dim3.runsim-Tuple{String, String}","page":"API Reference","title":"ElasticFDSG.dim3.runsim","text":"Elastic_FDSG.dim3.runsim(CONFIGPATH, VELMODPATH)\n\nRun the 3D elastic forward simulation using the specified configuration and velocity model.\n\nArguments\n\nCONFIGPATH::String: Path to the configuration.yaml file containing simulation settings.\nVELMODPATH::String: Path to the velocity model file. Supported formats include .jld2, .npy, and .npz.\n\nReturns\n\nNothing: The function runs the simulation and saves the results as specified in the configuration file.\n\n\n\n\n\n","category":"method"},{"location":"#ElasticFDSG.jl","page":"Home","title":"ElasticFDSG.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"ElasticFDSG.jl solves the elastic wave equation using finite differences on a staggered grid (FDSG) in the velocity-stress formulation [Virieux (1986)].","category":"page"},{"location":"","page":"Home","title":"Home","text":"But why yet another implementation? In many current applications, the installation of various dependencies, setup and definition of model parameters are often already challenging tasks, especially for inexperienced users seeking a quick and straightforward workflow. ElasticFDSG was developed with ease of use in mind.  It aims to offer a user-friendly experience while also maintaining flexibility to be applied to a wide variety of simulation scenarios. Users can easily customize their simulations by creating configuration files and defining velocity models in a straightforward manner, that can then be passed directly to the solvers. Furthermore, ElasticFDSG is implemented primarily in base Julia, minimizing dependencies on external packages to ensure the package remains robust over time.","category":"page"},{"location":"#Features","page":"Home","title":"Features","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"User friendly usage.\n2D and 3D elastic forward modelling on regular grids.\nMultithreaded CPU and parallized GPU (Cuda & Metal) solver.\nSpatial derivatives of order 1 to 10.\nFirst order time stepping.\nElastic isotropic or vertical transversal isotrop (VTI) 2D models using 2D Thomson parameter.\nElastic isotropic, VTI or orthorhombic (ORT) 3D models using Tsvankin parameter.\nAll elastic properties are defined on full integer grid points to simplify the creation of velocity models.\nSolver can handle fully heterogeneous media.\nPoint and anisotropic double couple sources. \nAbsorbing boundaries implemented with Convolutional-Perfectly-Matched-Layer.\nSave geophone receiver (velocity point sensors). \nSave snapshots at specified time steps.\nSave Distributed Acoustic Sensing (DAS) receiver, aligned with model coordinate axis (strain-profiles).\nEasy-to-read source code.","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"using Pkg\nPkg.add(ElasticFDSG)","category":"page"},{"location":"#Citing","page":"Home","title":"Citing","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"If you find this package helpful for your research, please consider citing:","category":"page"},{"location":"","page":"Home","title":"Home","text":"@misc{ElasticFDSG,\n  author       = {William Tegtow},\n  title        = {ElasticFDSG.jl: Simulating elastic wave propagation in 2D and 3D anisotropic media.},\n  year         = {2025},\n  publisher    = {GitHub},\n  journal      = {GitHub repository},\n  howpublished = {\\url{https://github.com/wtegtow/ElasticFDSG.jl}},\n  note         = {Version 1.0.0},\n  doi          = {https://doi.org/10.5281/zenodo.14752931}\n}\n","category":"page"}]
}
