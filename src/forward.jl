# ─────────────────────────────────────────────────────────────────────────────
# Staggered-grid finite-difference update kernels for elastic wave propagation.
#
# Array index convention
#   2D: field[i_x, i_z]         @index returns (x, z)
#   3D: field[i_x, i_y, i_z]    @index returns (x, y, z)
#
# Stiffness flat-matrix column layout:
#   2D: c11, c13, c33, c44, rho
#   3D: c11, c12, c13, c22, c23, c33, c44, c55, c66, rho
#
# CPML coefficient table layout (3 × n matrix):
#   row 1: a  
#   row 2: b  
#   row 3: K  
# ─────────────────────────────────────────────────────────────────────────────

# 2D: columns = (c11, c13, c33, c44, rho)
function _flatten_stiffness(tensors::Vector{Stiffness}, ::Val{2}, fp::DataType)
    n   = length(tensors)
    mat = Matrix{fp}(undef, n, 5)
    for (i, s) in enumerate(tensors)
        f = s.fields
        mat[i, 1] = fp(f.c11)
        mat[i, 2] = fp(f.c13)
        mat[i, 3] = fp(f.c33)
        mat[i, 4] = fp(f.c44)
        mat[i, 5] = fp(f.rho)
    end
    return mat
end

# 3D: columns = (c11, c12, c13, c22, c23, c33, c44, c55, c66, rho)
function _flatten_stiffness(tensors::Vector{Stiffness}, ::Val{3}, fp::DataType)
    n   = length(tensors)
    mat = Matrix{fp}(undef, n, 10)
    for (i, s) in enumerate(tensors)
        f = s.fields
        mat[i,  1] = fp(f.c11)
        mat[i,  2] = fp(f.c12)
        mat[i,  3] = fp(f.c13)
        mat[i,  4] = fp(f.c22)
        mat[i,  5] = fp(f.c23)
        mat[i,  6] = fp(f.c33)
        mat[i,  7] = fp(f.c44)
        mat[i,  8] = fp(f.c55)
        mat[i,  9] = fp(f.c66)
        mat[i, 10] = fp(f.rho)
    end
    return mat
end

function diff_coeff(M::Int)
    coeffs = if M == 1
        [1.0]
    elseif M == 2
        [0.1129136e+1, -0.4304542e-1]
    elseif M == 3
        [0.1186247e+1, -0.7266808e-1, 0.6351497e-2]
    elseif M == 4
        [0.1218159e+1, -0.9397218e-1, 0.1519043e-1, -0.1742128e-2]
    elseif M == 5
        [0.1236425e+1, -0.1081130e+0, 0.2339911e-1, -0.5061550e-2, 0.7054313e-3]
    elseif M == 6
        [0.1247576e+1, -0.1174969e+0, 0.2997288e-1, -0.8741572e-2, 0.2262285e-2, -0.3745306e-3]
    elseif M == 7
        [0.1254380e+1, -0.1235307e+0, 0.3467231e-1, -0.1192915e-1, 0.4057090e-2, -0.1191005e-2, 0.2263204e-3]
    elseif M == 8
        [0.1259012e+1, -0.1277647e+0, 0.3820715e-1, -0.1458251e-1, 0.5845385e-2, -0.2213861e-2, 0.7243880e-3, -0.1566173e-3]
    elseif M == 9
        [0.1262147e+1, -0.1306967e+0, 0.4075792e-1, -0.1665221e-1, 0.7377057e-2, -0.3258150e-2, 0.1336259e-2, -0.4775830e-3, 0.1151664e-3]
    elseif M == 10
        [0.1264362e+1, -0.1327958e+0, 0.4264687e-1, -0.1824918e-1, 0.8656223e-2, -0.4200034e-2, 0.1989180e-2, -0.8686637e-3, 0.3342741e-3, -0.8854090e-4]
    else
        error("diff_coeff: M must be between 1 and 10, got $M")
    end
    return coeffs
end

struct SimParams{N, T, Cfd<:AbstractVector{T}, Cdat<:AbstractMatrix{T}}
    N_fd   :: Int
    c_fd   :: Cfd
    c_data :: Cdat
end

function init_simparams(fdsg::FDSG)
    fp     = eval(Symbol(fdsg.config.dict["settings"]["precision"]))
    N      = Int(fdsg.config.dict["settings"]["spatial_derivative_order"])
    dim    = fdsg.config.dim
    AT     = fdsg.device.array
    c_fd   = AT(fp.(diff_coeff(N)))
    c_data = AT(_flatten_stiffness(fdsg.elastic.c_tensors, Val(dim), fp))
    return SimParams{dim, fp, typeof(c_fd), typeof(c_data)}(N, c_fd, c_data)
end

# ══════════════════════════════════════════════════════════════════════════════
# 2D KERNELS  
# ══════════════════════════════════════════════════════════════════════════════

# ──────────────────────────────────────────────────────────────────────────────
# sxx, szz  (normal stresses)
# guard:  N ≤ x ≤ nx−N,   N+1 ≤ z ≤ nz−N
# c_eff:  average along x  ([x,z] + [x+1,z])
# stencils:
#   vx_x → x_odd   (forward in x)
#   vz_z → z_evn   (backward in z)
# ──────────────────────────────────────────────────────────────────────────────
@kernel inbounds=true function _sxx_szz_2d!(
        sxx, szz, vx, vz,
        c_data, c_lookup, pml_lookup,
        vx_x_pml, vz_z_pml,
        x_odd, z_evn,
        c_fd, dx::T, dz::T, dt::T, N::Int, nx::Int, nz::Int) where T

    x, z = @index(Global, NTuple)

    if N <= x <= nx-N && N+1 <= z <= nz-N

        cidx    = c_lookup[x,   z]
        cidx_xp = c_lookup[x+1, z]
        c11eff  = T(0.5) * (c_data[cidx, 1] + c_data[cidx_xp, 1])
        c13eff  = T(0.5) * (c_data[cidx, 2] + c_data[cidx_xp, 2])
        c33eff  = T(0.5) * (c_data[cidx, 3] + c_data[cidx_xp, 3])

        vx_x = zero(T)
        vz_z = zero(T)
        for i in 1:N
            vx_x += c_fd[i]/dx * (vx[x+i, z] - vx[x-(i-1), z])
            vz_z += c_fd[i]/dz * (vz[x, z+(i-1)] - vz[x, z-i])
        end

        if pml_lookup[x, z] > 0
            pml_idx = pml_lookup[x, z]
            vx_x_pml[pml_idx] = x_odd[2, x] * vx_x_pml[pml_idx] + x_odd[1, x] * vx_x
            vz_z_pml[pml_idx] = z_evn[2, z] * vz_z_pml[pml_idx] + z_evn[1, z] * vz_z
            vx_x = vx_x / x_odd[3, x] + vx_x_pml[pml_idx]
            vz_z = vz_z / z_evn[3, z] + vz_z_pml[pml_idx]
        end

        sxx[x, z] += dt * (c11eff * vx_x + c13eff * vz_z)
        szz[x, z] += dt * (c13eff * vx_x + c33eff * vz_z)
    end
end

# ──────────────────────────────────────────────────────────────────────────────
# sxz  (shear stress)
# guard:  N+1 ≤ x ≤ nx−N+1,   N ≤ z ≤ nz−N
# c_eff:  average along z  ([x,z] + [x,z+1])
# stencils:
#   vz_x → x_evn   (backward in x)
#   vx_z → z_odd   (forward in z)
# ──────────────────────────────────────────────────────────────────────────────
@kernel inbounds=true function _sxz_2d!(
        sxz, vx, vz,
        c_data, c_lookup, pml_lookup,
        vz_x_pml, vx_z_pml,
        x_evn, z_odd,
        c_fd, dx::T, dz::T, dt::T, N::Int, nx::Int, nz::Int) where T

    x, z = @index(Global, NTuple)

    if N+1 <= x <= nx-N+1 && N <= z <= nz-N

        cidx    = c_lookup[x, z  ]
        cidx_zp = c_lookup[x, z+1]
        c44eff  = T(0.5) * (c_data[cidx, 4] + c_data[cidx_zp, 4])

        vz_x = zero(T)
        vx_z = zero(T)
        for i in 1:N
            vz_x += c_fd[i]/dx * (vz[x+(i-1), z] - vz[x-i, z])
            vx_z += c_fd[i]/dz * (vx[x, z+i] - vx[x, z-(i-1)])
        end

        if pml_lookup[x, z] > 0
            pml_idx = pml_lookup[x, z]
            vz_x_pml[pml_idx] = x_evn[2, x] * vz_x_pml[pml_idx] + x_evn[1, x] * vz_x
            vx_z_pml[pml_idx] = z_odd[2, z] * vx_z_pml[pml_idx] + z_odd[1, z] * vx_z
            vz_x = vz_x / x_evn[3, x] + vz_x_pml[pml_idx]
            vx_z = vx_z / z_odd[3, z] + vx_z_pml[pml_idx]
        end

        sxz[x, z] += dt * c44eff * (vz_x + vx_z)
    end
end

# ──────────────────────────────────────────────────────────────────────────────
# vx  (x-velocity)
# guard:  N+1 ≤ x ≤ nx−N+1,   N+1 ≤ z ≤ nz−N+1
# rho:    direct lookup, no averaging
# stencils:
#   sxx_x → x_evn   (backward in x)
#   sxz_z → z_evn   (backward in z)
# ──────────────────────────────────────────────────────────────────────────────
@kernel inbounds=true function _vx_2d!(
        vx, sxx, sxz,
        c_data, c_lookup, pml_lookup,
        sxx_x_pml, sxz_z_pml,
        x_evn, z_evn,
        c_fd, dx::T, dz::T, dt::T, N::Int, nx::Int, nz::Int) where T

    x, z = @index(Global, NTuple)

    if N+1 <= x <= nx-N+1 && N+1 <= z <= nz-N+1

        cidx = c_lookup[x, z]
        rho  = c_data[cidx, 5]

        sxx_x = zero(T)
        sxz_z = zero(T)
        for i in 1:N
            sxx_x += c_fd[i]/dx * (sxx[x+(i-1), z] - sxx[x-i, z])
            sxz_z += c_fd[i]/dz * (sxz[x, z+(i-1)] - sxz[x, z-i])
        end

        if pml_lookup[x, z] > 0
            pml_idx = pml_lookup[x, z]
            sxx_x_pml[pml_idx] = x_evn[2, x] * sxx_x_pml[pml_idx] + x_evn[1, x] * sxx_x
            sxz_z_pml[pml_idx] = z_evn[2, z] * sxz_z_pml[pml_idx] + z_evn[1, z] * sxz_z
            sxx_x = sxx_x / x_evn[3, x] + sxx_x_pml[pml_idx]
            sxz_z = sxz_z / z_evn[3, z] + sxz_z_pml[pml_idx]
        end

        vx[x, z] += dt/rho * (sxx_x + sxz_z)
    end
end

# ──────────────────────────────────────────────────────────────────────────────
# vz  (z-velocity)
# guard:  N ≤ x ≤ nx−N,   N ≤ z ≤ nz−N
# rho:    4-point average [x,z],[x+1,z],[x,z+1],[x+1,z+1]
# stencils:
#   sxz_x → x_odd   (forward in x)
#   szz_z → z_odd   (forward in z)
# ──────────────────────────────────────────────────────────────────────────────
@kernel inbounds=true function _vz_2d!(
        vz, sxz, szz,
        c_data, c_lookup, pml_lookup,
        sxz_x_pml, szz_z_pml,
        x_odd, z_odd,
        c_fd, dx::T, dz::T, dt::T, N::Int, nx::Int, nz::Int) where T

    x, z = @index(Global, NTuple)

    if N <= x <= nx-N && N <= z <= nz-N

        cidx      = c_lookup[x,   z  ]
        cidx_xp   = c_lookup[x+1, z  ]
        cidx_zp   = c_lookup[x,   z+1]
        cidx_xpzp = c_lookup[x+1, z+1]
        rhoeff = T(0.25) * (c_data[cidx, 5] + c_data[cidx_xp, 5] +
                             c_data[cidx_zp, 5] + c_data[cidx_xpzp, 5])

        sxz_x = zero(T)
        szz_z = zero(T)
        for i in 1:N
            sxz_x += c_fd[i]/dx * (sxz[x+i, z] - sxz[x-(i-1), z])
            szz_z += c_fd[i]/dz * (szz[x, z+i] - szz[x, z-(i-1)])
        end

        if pml_lookup[x, z] > 0
            pml_idx = pml_lookup[x, z]
            sxz_x_pml[pml_idx] = x_odd[2, x] * sxz_x_pml[pml_idx] + x_odd[1, x] * sxz_x
            szz_z_pml[pml_idx] = z_odd[2, z] * szz_z_pml[pml_idx] + z_odd[1, z] * szz_z
            sxz_x = sxz_x / x_odd[3, x] + sxz_x_pml[pml_idx]
            szz_z = szz_z / z_odd[3, z] + szz_z_pml[pml_idx]
        end

        vz[x, z] += dt/rhoeff * (sxz_x + szz_z)
    end
end


# ══════════════════════════════════════════════════════════════════════════════
# 3D KERNELS  
# ══════════════════════════════════════════════════════════════════════════════

# ──────────────────────────────────────────────────────────────────────────────
# sxx, syy, szz  (normal stresses)
# guard:  N ≤ x ≤ nx−N,   N+1 ≤ y ≤ ny−N,   N+1 ≤ z ≤ nz−N
# c_eff:  average along x  ([x,y,z] + [x+1,y,z])
# stencils:
#   vx_x → x_odd   (forward in x)
#   vy_y → y_evn   (backward in y)
#   vz_z → z_evn   (backward in z)
# ──────────────────────────────────────────────────────────────────────────────
@kernel inbounds=true function _sxx_syy_szz_3d!(
        sxx, syy, szz, vx, vy, vz,
        c_data, c_lookup, pml_lookup,
        vx_x_pml, vy_y_pml, vz_z_pml,
        x_odd, y_evn, z_evn,
        c_fd, dx::T, dy::T, dz::T, dt::T, N::Int,
        nx::Int, ny::Int, nz::Int) where T

    x, y, z = @index(Global, NTuple)

    if N <= x <= nx-N && N+1 <= y <= ny-N && N+1 <= z <= nz-N

        cidx    = c_lookup[x,   y, z]
        cidx_xp = c_lookup[x+1, y, z]
        c11eff  = T(0.5) * (c_data[cidx, 1] + c_data[cidx_xp, 1])
        c12eff  = T(0.5) * (c_data[cidx, 2] + c_data[cidx_xp, 2])
        c13eff  = T(0.5) * (c_data[cidx, 3] + c_data[cidx_xp, 3])
        c22eff  = T(0.5) * (c_data[cidx, 4] + c_data[cidx_xp, 4])
        c23eff  = T(0.5) * (c_data[cidx, 5] + c_data[cidx_xp, 5])
        c33eff  = T(0.5) * (c_data[cidx, 6] + c_data[cidx_xp, 6])

        vx_x = zero(T);  vy_y = zero(T);  vz_z = zero(T)
        for i in 1:N
            vx_x += c_fd[i]/dx * (vx[x+i, y, z] - vx[x-(i-1), y, z])
            vy_y += c_fd[i]/dy * (vy[x, y+(i-1), z] - vy[x, y-i, z])
            vz_z += c_fd[i]/dz * (vz[x, y, z+(i-1)] - vz[x, y, z-i])
        end

        if pml_lookup[x, y, z] > 0
            pml_idx = pml_lookup[x, y, z]
            vx_x_pml[pml_idx] = x_odd[2, x] * vx_x_pml[pml_idx] + x_odd[1, x] * vx_x
            vy_y_pml[pml_idx] = y_evn[2, y] * vy_y_pml[pml_idx] + y_evn[1, y] * vy_y
            vz_z_pml[pml_idx] = z_evn[2, z] * vz_z_pml[pml_idx] + z_evn[1, z] * vz_z
            vx_x = vx_x / x_odd[3, x] + vx_x_pml[pml_idx]
            vy_y = vy_y / y_evn[3, y] + vy_y_pml[pml_idx]
            vz_z = vz_z / z_evn[3, z] + vz_z_pml[pml_idx]
        end

        sxx[x, y, z] += dt * (c11eff * vx_x + c12eff * vy_y + c13eff * vz_z)
        syy[x, y, z] += dt * (c12eff * vx_x + c22eff * vy_y + c23eff * vz_z)
        szz[x, y, z] += dt * (c13eff * vx_x + c23eff * vy_y + c33eff * vz_z)
    end
end

# ──────────────────────────────────────────────────────────────────────────────
# sxy  (XY-shear stress)
# guard:  N+1 ≤ x ≤ nx−N+1,   N ≤ y ≤ ny−N,   1 ≤ z ≤ nz
# c_eff:  average along y  ([x,y,z] + [x,y+1,z])
# stencils:
#   vy_x → x_evn   (backward in x)
#   vx_y → y_odd   (forward in y)
# ──────────────────────────────────────────────────────────────────────────────
@kernel inbounds=true function _sxy_3d!(
        sxy, vx, vy,
        c_data, c_lookup, pml_lookup,
        vy_x_pml, vx_y_pml,
        x_evn, y_odd,
        c_fd, dx::T, dy::T, dt::T, N::Int,
        nx::Int, ny::Int, nz::Int) where T

    x, y, z = @index(Global, NTuple)

    if N+1 <= x <= nx-N+1 && N <= y <= ny-N && 1 <= z <= nz

        cidx    = c_lookup[x, y,   z]
        cidx_yp = c_lookup[x, y+1, z]
        c66eff  = T(0.5) * (c_data[cidx, 9] + c_data[cidx_yp, 9])

        vy_x = zero(T);  vx_y = zero(T)
        for i in 1:N
            vy_x += c_fd[i]/dx * (vy[x+(i-1), y, z] - vy[x-i, y, z])
            vx_y += c_fd[i]/dy * (vx[x, y+i, z] - vx[x, y-(i-1), z])
        end

        if pml_lookup[x, y, z] > 0
            pml_idx = pml_lookup[x, y, z]
            vy_x_pml[pml_idx] = x_evn[2, x] * vy_x_pml[pml_idx] + x_evn[1, x] * vy_x
            vx_y_pml[pml_idx] = y_odd[2, y] * vx_y_pml[pml_idx] + y_odd[1, y] * vx_y
            vy_x = vy_x / x_evn[3, x] + vy_x_pml[pml_idx]
            vx_y = vx_y / y_odd[3, y] + vx_y_pml[pml_idx]
        end

        sxy[x, y, z] += dt * c66eff * (vy_x + vx_y)
    end
end

# ──────────────────────────────────────────────────────────────────────────────
# sxz  (XZ-shear stress)
# guard:  N+1 ≤ x ≤ nx−N+1,   1 ≤ y ≤ ny,   N ≤ z ≤ nz−N
# c_eff:  average along z  ([x,y,z] + [x,y,z+1])
# stencils:
#   vz_x → x_evn   (backward in x)
#   vx_z → z_odd   (forward in z)
# ──────────────────────────────────────────────────────────────────────────────
@kernel inbounds=true function _sxz_3d!(
        sxz, vx, vz,
        c_data, c_lookup, pml_lookup,
        vz_x_pml, vx_z_pml,
        x_evn, z_odd,
        c_fd, dx::T, dz::T, dt::T, N::Int,
        nx::Int, ny::Int, nz::Int) where T

    x, y, z = @index(Global, NTuple)

    if N+1 <= x <= nx-N+1 && 1 <= y <= ny && N <= z <= nz-N

        cidx    = c_lookup[x, y, z  ]
        cidx_zp = c_lookup[x, y, z+1]
        c55eff  = T(0.5) * (c_data[cidx, 8] + c_data[cidx_zp, 8])

        vz_x = zero(T);  vx_z = zero(T)
        for i in 1:N
            vz_x += c_fd[i]/dx * (vz[x+(i-1), y, z] - vz[x-i, y, z])
            vx_z += c_fd[i]/dz * (vx[x, y, z+i] - vx[x, y, z-(i-1)])
        end

        if pml_lookup[x, y, z] > 0
            pml_idx = pml_lookup[x, y, z]
            vz_x_pml[pml_idx] = x_evn[2, x] * vz_x_pml[pml_idx] + x_evn[1, x] * vz_x
            vx_z_pml[pml_idx] = z_odd[2, z] * vx_z_pml[pml_idx] + z_odd[1, z] * vx_z
            vz_x = vz_x / x_evn[3, x] + vz_x_pml[pml_idx]
            vx_z = vx_z / z_odd[3, z] + vx_z_pml[pml_idx]
        end

        sxz[x, y, z] += dt * c55eff * (vz_x + vx_z)
    end
end

# ──────────────────────────────────────────────────────────────────────────────
# syz  (YZ-shear stress)
# guard:  1 ≤ x ≤ nx,   N ≤ y ≤ ny−N,   N ≤ z ≤ nz−N
# c_eff:  4-point c44 average [x,y,z],[x,y+1,z],[x,y,z+1],[x,y+1,z+1]
# stencils:
#   vz_y → y_odd   (forward in y)
#   vy_z → z_odd   (forward in z)
# ──────────────────────────────────────────────────────────────────────────────
@kernel inbounds=true function _syz_3d!(
        syz, vy, vz,
        c_data, c_lookup, pml_lookup,
        vz_y_pml, vy_z_pml,
        y_odd, z_odd,
        c_fd, dy::T, dz::T, dt::T, N::Int,
        nx::Int, ny::Int, nz::Int) where T

    x, y, z = @index(Global, NTuple)

    if 1 <= x <= nx && N <= y <= ny-N && N <= z <= nz-N

        cidx      = c_lookup[x, y, z]
        cidx_yp   = c_lookup[x, y+1, z]
        cidx_zp   = c_lookup[x, y, z+1]
        cidx_ypzp = c_lookup[x, y+1, z+1]
        c44eff = T(0.25) * (c_data[cidx, 7] + c_data[cidx_yp, 7] +
                             c_data[cidx_zp, 7] + c_data[cidx_ypzp, 7])

        vz_y = zero(T);  vy_z = zero(T)
        for i in 1:N
            vz_y += c_fd[i]/dy * (vz[x, y+i, z] - vz[x, y-(i-1), z])
            vy_z += c_fd[i]/dz * (vy[x, y, z+i] - vy[x, y, z-(i-1)])
        end

        if pml_lookup[x, y, z] > 0
            pml_idx = pml_lookup[x, y, z]
            vz_y_pml[pml_idx] = y_odd[2, y] * vz_y_pml[pml_idx] + y_odd[1, y] * vz_y
            vy_z_pml[pml_idx] = z_odd[2, z] * vy_z_pml[pml_idx] + z_odd[1, z] * vy_z
            vz_y = vz_y / y_odd[3, y] + vz_y_pml[pml_idx]
            vy_z = vy_z / z_odd[3, z] + vy_z_pml[pml_idx]
        end

        syz[x, y, z] += dt * c44eff * (vz_y + vy_z)
    end
end

# ──────────────────────────────────────────────────────────────────────────────
# vx  (x-velocity in 3D)
# guard:  N+1 ≤ x ≤ nx−N+1,   N+1 ≤ y ≤ ny−N+1,   N+1 ≤ z ≤ nz−N+1
# rho:    direct lookup
# stencils:
#   sxx_x → x_evn   (backward in x)
#   sxy_y → y_evn   (backward in y)
#   sxz_z → z_evn   (backward in z)
# ──────────────────────────────────────────────────────────────────────────────
@kernel inbounds=true function _vx_3d!(
        vx, sxx, sxy, sxz,
        c_data, c_lookup, pml_lookup,
        sxx_x_pml, sxy_y_pml, sxz_z_pml,
        x_evn, y_evn, z_evn,
        c_fd, dx::T, dy::T, dz::T, dt::T, N::Int,
        nx::Int, ny::Int, nz::Int) where T

    x, y, z = @index(Global, NTuple)

    if N+1 <= x <= nx-N+1 && N+1 <= y <= ny-N+1 && N+1 <= z <= nz-N+1

        cidx = c_lookup[x, y, z]
        rho  = c_data[cidx, 10]

        sxx_x = zero(T);  sxy_y = zero(T);  sxz_z = zero(T)
        for i in 1:N
            sxx_x += c_fd[i]/dx * (sxx[x+(i-1), y, z] - sxx[x-i, y, z])
            sxy_y += c_fd[i]/dy * (sxy[x, y+(i-1),  z] - sxy[x, y-i, z])
            sxz_z += c_fd[i]/dz * (sxz[x, y, z+(i-1) ] - sxz[x, y, z-i])
        end

        if pml_lookup[x, y, z] > 0
            pml_idx = pml_lookup[x, y, z]
            sxx_x_pml[pml_idx] = x_evn[2, x] * sxx_x_pml[pml_idx] + x_evn[1, x] * sxx_x
            sxy_y_pml[pml_idx] = y_evn[2, y] * sxy_y_pml[pml_idx] + y_evn[1, y] * sxy_y
            sxz_z_pml[pml_idx] = z_evn[2, z] * sxz_z_pml[pml_idx] + z_evn[1, z] * sxz_z
            sxx_x = sxx_x / x_evn[3, x] + sxx_x_pml[pml_idx]
            sxy_y = sxy_y / y_evn[3, y] + sxy_y_pml[pml_idx]
            sxz_z = sxz_z / z_evn[3, z] + sxz_z_pml[pml_idx]
        end

        vx[x, y, z] += dt/rho * (sxx_x + sxy_y + sxz_z)
    end
end

# ──────────────────────────────────────────────────────────────────────────────
# vy  (y-velocity in 3D)
# guard:  N ≤ x ≤ nx−N,   N ≤ y ≤ ny−N,   N+1 ≤ z ≤ nz−N+1
# rho:    4-point average [x,y,z],[x,y+1,z],[x+1,y,z],[x+1,y+1,z]
# stencils:
#   sxy_x → x_odd   (forward in x)
#   syy_y → y_odd   (forward in y)
#   syz_z → z_evn   (backward in z)
# ──────────────────────────────────────────────────────────────────────────────
@kernel inbounds=true function _vy_3d!(
        vy, sxy, syy, syz,
        c_data, c_lookup, pml_lookup,
        sxy_x_pml, syy_y_pml, syz_z_pml,
        x_odd, y_odd, z_evn,
        c_fd, dx::T, dy::T, dz::T, dt::T, N::Int,
        nx::Int, ny::Int, nz::Int) where T

    x, y, z = @index(Global, NTuple)

    if N <= x <= nx-N && N <= y <= ny-N && N+1 <= z <= nz-N+1

        cidx      = c_lookup[x,   y,   z]
        cidx_yp   = c_lookup[x,   y+1, z]
        cidx_xp   = c_lookup[x+1, y,   z]
        cidx_xpyp = c_lookup[x+1, y+1, z]
        rhoeff = T(0.25) * (c_data[cidx, 10] + c_data[cidx_yp, 10] +
                             c_data[cidx_xp, 10] + c_data[cidx_xpyp, 10])

        sxy_x = zero(T);  syy_y = zero(T);  syz_z = zero(T)
        for i in 1:N
            sxy_x += c_fd[i]/dx * (sxy[x+i, y, z] - sxy[x-(i-1), y, z])
            syy_y += c_fd[i]/dy * (syy[x, y+i, z] - syy[x, y-(i-1), z])
            syz_z += c_fd[i]/dz * (syz[x, y, z+(i-1)] - syz[x, y, z-i])
        end

        if pml_lookup[x, y, z] > 0
            pml_idx = pml_lookup[x, y, z]
            sxy_x_pml[pml_idx] = x_odd[2, x] * sxy_x_pml[pml_idx] + x_odd[1, x] * sxy_x
            syy_y_pml[pml_idx] = y_odd[2, y] * syy_y_pml[pml_idx] + y_odd[1, y] * syy_y
            syz_z_pml[pml_idx] = z_evn[2, z] * syz_z_pml[pml_idx] + z_evn[1, z] * syz_z
            sxy_x = sxy_x / x_odd[3, x] + sxy_x_pml[pml_idx]
            syy_y = syy_y / y_odd[3, y] + syy_y_pml[pml_idx]
            syz_z = syz_z / z_evn[3, z] + syz_z_pml[pml_idx]
        end

        vy[x, y, z] += dt/rhoeff * (sxy_x + syy_y + syz_z)
    end
end

# ──────────────────────────────────────────────────────────────────────────────
# vz  (z-velocity in 3D)
# guard:  N ≤ x ≤ nx−N,   N+1 ≤ y ≤ ny−N+1,   N ≤ z ≤ nz−N
# rho:    4-point average [x,y,z],[x+1,y,z],[x,y,z+1],[x+1,y,z+1]
# stencils:
#   sxz_x → x_odd   (forward in x)
#   syz_y → y_evn   (backward in y)
#   szz_z → z_odd   (forward in z)
# ──────────────────────────────────────────────────────────────────────────────
@kernel inbounds=true function _vz_3d!(
        vz, sxz, syz, szz,
        c_data, c_lookup, pml_lookup,
        sxz_x_pml, syz_y_pml, szz_z_pml,
        x_odd, y_evn, z_odd,
        c_fd, dx::T, dy::T, dz::T, dt::T, N::Int,
        nx::Int, ny::Int, nz::Int) where T

    x, y, z = @index(Global, NTuple)

    if N <= x <= nx-N && N+1 <= y <= ny-N+1 && N <= z <= nz-N

        cidx      = c_lookup[x,   y, z  ]
        cidx_xp   = c_lookup[x+1, y, z  ]
        cidx_zp   = c_lookup[x,   y, z+1]
        cidx_xpzp = c_lookup[x+1, y, z+1]
        rhoeff = T(0.25) * (c_data[cidx, 10] + c_data[cidx_xp, 10] +
                             c_data[cidx_zp, 10] + c_data[cidx_xpzp, 10])

        sxz_x = zero(T);  syz_y = zero(T);  szz_z = zero(T)
        for i in 1:N
            sxz_x += c_fd[i]/dx * (sxz[x+i, y, z] - sxz[x-(i-1), y, z])
            syz_y += c_fd[i]/dy * (syz[x, y+(i-1), z] - syz[x, y-i, z])
            szz_z += c_fd[i]/dz * (szz[x, y, z+i] - szz[x, y, z-(i-1)])
        end

        if pml_lookup[x, y, z] > 0
            pml_idx = pml_lookup[x, y, z]
            sxz_x_pml[pml_idx] = x_odd[2, x] * sxz_x_pml[pml_idx] + x_odd[1, x] * sxz_x
            syz_y_pml[pml_idx] = y_evn[2, y] * syz_y_pml[pml_idx] + y_evn[1, y] * syz_y
            szz_z_pml[pml_idx] = z_odd[2, z] * szz_z_pml[pml_idx] + z_odd[1, z] * szz_z
            sxz_x = sxz_x / x_odd[3, x] + sxz_x_pml[pml_idx]
            syz_y = syz_y / y_evn[3, y] + syz_y_pml[pml_idx]
            szz_z = szz_z / z_odd[3, z] + szz_z_pml[pml_idx]
        end

        vz[x, y, z] += dt/rhoeff * (sxz_x + syz_y + szz_z)
    end
end


# ══════════════════════════════════════════════════════════════════════════════
# DISPATCHER 
# ══════════════════════════════════════════════════════════════════════════════
T_spacing(coords::AbstractVector) = eltype(coords)(step(coords))

function _launch(kernel_fn, backend, block_size, ndrange, args...)
    k = kernel_fn(backend, block_size)
    k(args...; ndrange=ndrange)
end

function update_velocities!(fields::Fields2D, pml::CPML2D,
                             elastic::Elastic{2}, domain::Domain{2},
                             time::SimTime, params::SimParams{2},
                             backend, block_size)
    nx, nz = domain.shape
    dx  = T_spacing(domain.coordinates[1])
    dz  = T_spacing(domain.coordinates[2])
    dt  = time.dt
    N   = params.N_fd
    c_fd  = params.c_fd
    c_dat = params.c_data
    clu   = elastic.c_lookup
    plu   = domain.pml_lookup
    nr    = (nx, nz)

    _launch(_vx_2d!, backend, block_size, nr,
            fields.vx, fields.sxx, fields.sxz,
            c_dat, clu, plu,
            pml.sxx_x, pml.sxz_z,
            pml.x_evn, pml.z_evn,
            c_fd, dx, dz, dt, N, nx, nz)

    _launch(_vz_2d!, backend, block_size, nr,
            fields.vz, fields.sxz, fields.szz,
            c_dat, clu, plu,
            pml.sxz_x, pml.szz_z,
            pml.x_odd, pml.z_odd,
            c_fd, dx, dz, dt, N, nx, nz)
end


function update_stresses!(fields::Fields2D, pml::CPML2D,
                           elastic::Elastic{2}, domain::Domain{2},
                           time::SimTime, params::SimParams{2},
                           backend, block_size)
    nx, nz = domain.shape
    dx  = T_spacing(domain.coordinates[1])
    dz  = T_spacing(domain.coordinates[2])
    dt  = time.dt
    N   = params.N_fd
    c_fd  = params.c_fd
    c_dat = params.c_data
    clu   = elastic.c_lookup
    plu   = domain.pml_lookup
    nr    = (nx, nz)

    _launch(_sxx_szz_2d!, backend, block_size, nr,
            fields.sxx, fields.szz, fields.vx, fields.vz,
            c_dat, clu, plu,
            pml.vx_x, pml.vz_z,
            pml.x_odd, pml.z_evn,
            c_fd, dx, dz, dt, N, nx, nz)

    _launch(_sxz_2d!, backend, block_size, nr,
            fields.sxz, fields.vx, fields.vz,
            c_dat, clu, plu,
            pml.vz_x, pml.vx_z,
            pml.x_evn, pml.z_odd,
            c_fd, dx, dz, dt, N, nx, nz)
end


function update_velocities!(fields::Fields3D, pml::CPML3D,
                             elastic::Elastic{3}, domain::Domain{3},
                             time::SimTime, params::SimParams{3},
                             backend, block_size)
    nx, ny, nz = domain.shape
    dx  = T_spacing(domain.coordinates[1])
    dy  = T_spacing(domain.coordinates[2])
    dz  = T_spacing(domain.coordinates[3])
    dt  = time.dt
    N   = params.N_fd
    c_fd  = params.c_fd
    c_dat = params.c_data
    clu   = elastic.c_lookup
    plu   = domain.pml_lookup
    nr    = (nx, ny, nz)

    _launch(_vx_3d!, backend, block_size, nr,
            fields.vx, fields.sxx, fields.sxy, fields.sxz,
            c_dat, clu, plu,
            pml.sxx_x, pml.sxy_y, pml.sxz_z,
            pml.x_evn, pml.y_evn, pml.z_evn,
            c_fd, dx, dy, dz, dt, N, nx, ny, nz)

    _launch(_vy_3d!, backend, block_size, nr,
            fields.vy, fields.sxy, fields.syy, fields.syz,
            c_dat, clu, plu,
            pml.sxy_x, pml.syy_y, pml.syz_z,
            pml.x_odd, pml.y_odd, pml.z_evn,
            c_fd, dx, dy, dz, dt, N, nx, ny, nz)

    _launch(_vz_3d!, backend, block_size, nr,
            fields.vz, fields.sxz, fields.syz, fields.szz,
            c_dat, clu, plu,
            pml.sxz_x, pml.syz_y, pml.szz_z,
            pml.x_odd, pml.y_evn, pml.z_odd,
            c_fd, dx, dy, dz, dt, N, nx, ny, nz)
end


function update_stresses!(fields::Fields3D, pml::CPML3D,
                           elastic::Elastic{3}, domain::Domain{3},
                           time::SimTime, params::SimParams{3},
                           backend, block_size)
    nx, ny, nz = domain.shape
    dx  = T_spacing(domain.coordinates[1])
    dy  = T_spacing(domain.coordinates[2])
    dz  = T_spacing(domain.coordinates[3])
    dt  = time.dt
    N   = params.N_fd
    c_fd  = params.c_fd
    c_dat = params.c_data
    clu   = elastic.c_lookup
    plu   = domain.pml_lookup
    nr    = (nx, ny, nz)

    _launch(_sxx_syy_szz_3d!, backend, block_size, nr,
            fields.sxx, fields.syy, fields.szz,
            fields.vx, fields.vy, fields.vz,
            c_dat, clu, plu,
            pml.vx_x, pml.vy_y, pml.vz_z,
            pml.x_odd, pml.y_evn, pml.z_evn,
            c_fd, dx, dy, dz, dt, N, nx, ny, nz)

    _launch(_sxy_3d!, backend, block_size, nr,
            fields.sxy, fields.vx, fields.vy,
            c_dat, clu, plu,
            pml.vy_x, pml.vx_y,
            pml.x_evn, pml.y_odd,
            c_fd, dx, dy, dt, N, nx, ny, nz)

    _launch(_sxz_3d!, backend, block_size, nr,
            fields.sxz, fields.vx, fields.vz,
            c_dat, clu, plu,
            pml.vz_x, pml.vx_z,
            pml.x_evn, pml.z_odd,
            c_fd, dx, dz, dt, N, nx, ny, nz)

    _launch(_syz_3d!, backend, block_size, nr,
            fields.syz, fields.vy, fields.vz,
            c_dat, clu, plu,
            pml.vz_y, pml.vy_z,
            pml.y_odd, pml.z_odd,
            c_fd, dy, dz, dt, N, nx, ny, nz)
end
