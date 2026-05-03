using ElasticFDSG
ElasticFDSG.devmode!(true)
using Test, LinearAlgebra, Statistics
using JLD2, YAML, Einsum, UnPack

println("Running DC source validation test...")

# velocity model
h = 10.0
xcoords = 0:h:1000
ycoords = 0:h:300
zcoords = 0:h:1000
nx, ny, nz = length(xcoords), length(ycoords), length(zcoords)

X = getindex.(Iterators.product(xcoords, ycoords, zcoords), 1)
Y = getindex.(Iterators.product(xcoords, ycoords, zcoords), 2)
Z = getindex.(Iterators.product(xcoords, ycoords, zcoords), 3)

velmod = zeros(Float32, 13, nx, ny, nz)
velmod[1,:,:,:] .= X
velmod[2,:,:,:] .= Y
velmod[3,:,:,:] .= Z
velmod[4,:,:,:] .= 3000
velmod[5,:,:,:] .= 2000
velmod[6,:,:,:] .= 2500

# receiver array
src_x, src_y, src_z = (50, 100, 500)
rcv = [Dict("x" => 950, "y" => 200, "z" => z) for z in 50:100:950]

# config
config = Dict(
    "settings" => Dict(
        "device"                   => "cpu",
        "precision"                => "Float32",
        "spatial_derivative_order" => 4,
        "verbose"                  => false,
        "output_file"              => nothing
    ),
    "time" => Dict(
        "start"    => 0.0,
        "end"      => 0.7,
        "timestep" => 0.001
    ),
    "source" => Dict(
        "dominant_frequency" => 20,
        "wavelet_type"       => "ricker",
        "wavelet_center"     => 0.0625,
        "location"           => Dict("x" => src_x, "y" => src_y, "z" => src_z),
        "seismic_moment"     => 1e15,
        "moment_tensor"      => Dict(
            "Mxx" => 0, "Myy" => 0, "Mzz" => 0,
            "Mxy" => 0, "Mxz" => 0, "Myz" => 0.0,
            "anisotropic" => false
        )
    ),
    "boundaries" => Dict(
        "xstart"    => "absorbing",
        "xend"      => "absorbing",
        "ystart"    => "absorbing",
        "yend"      => "absorbing",
        "zstart"    => "absorbing",
        "zend"      => "absorbing",
        "pml_layer" => 10
    ),
    "receivers" => Dict(
        "geophones" => rcv,
        "das"       => Dict("x_aligned" => [], "y_aligned" => [], "z_aligned" => []),
        "snapshots" => Dict(
            "plane_positions" => [Dict("x" => src_x, "y" => src_y, "z" => src_z)],
            "times"           => Float64[0.1],
            "fields"          => ["vx"]
        )
    )
)

# DC moment tensor
M_dc = normalize([1 0 0; 0 0 0; 0 0 -1])
mt_indices = Dict("Mxx"=>(1,1), "Myy"=>(2,2), "Mzz"=>(3,3),
                  "Mxy"=>(1,2), "Mxz"=>(1,3), "Myz"=>(2,3))
for (key, (i, j)) in mt_indices
    config["source"]["moment_tensor"][key] = M_dc[i, j]
end

# run simulation
sim_dc = ElasticFDSG.runsim(config, velmod)

# analytical solution helpers
STF(t, t0, f0)    = @. (1 - 2*(π*f0)^2*(t-t0)^2) * exp(-(π*f0)^2*(t-t0)^2)
STF_D1(t, t0, f0) = @. -2*(π*f0)^2*(t-t0)*(3 - 2*(π*f0*(t-t0))^2) * exp(-(π*f0*(t-t0))^2)

function STF_CONV(t, t0, f0, a, b)
    Δt = t[2] - t[1]
    nt = length(t)
    res = zeros(nt)
    for i in 1:nt, τ in a:Δt:b
        res[i] += (1 - 2*(π*f0)^2*(t[i]-t0-τ)^2) * exp(-(π*f0)^2*(t[i]-t0-τ)^2) * Δt * τ
    end
    return res
end

function analytical_displacement(pos, t, t0, fdom, M0, ξ, MT, vp, vs, rho)
    δ(i, j) = i == j ? 1 : 0
    nt = length(t)

    r_vec = pos - ξ
    r = norm(r_vec)
    r < eps() && return zeros(3, nt)
    γ = r_vec / r

    FNE = M0 / (4π * rho * r^4)
    FIP = M0 / (4π * rho * r^2 * vp^2)
    FIS = M0 / (4π * rho * r^2 * vs^2)
    FFP = M0 / (4π * rho * r   * vp^3)
    FFS = M0 / (4π * rho * r   * vs^3)

    RNE = zeros(3); RIP = zeros(3); RIS = zeros(3); RFP = zeros(3); RFS = zeros(3)
    @einsum RNE[n] = (15*γ[n]*γ[p]*γ[q] - 3*γ[n]*δ(p,q) - 3*γ[p]*δ(n,q) - 3*γ[q]*δ(n,p)) * MT[p,q]
    @einsum RIP[n] = (6*γ[n]*γ[p]*γ[q] - γ[n]*δ(p,q) - γ[p]*δ(n,q) - γ[q]*δ(n,p)) * MT[p,q]
    @einsum RIS[n] = -(6*γ[n]*γ[p]*γ[q] - γ[n]*δ(p,q) - γ[p]*δ(n,q) - 2*γ[q]*δ(n,p)) * MT[p,q]
    @einsum RFP[n] = (γ[n]*γ[p]*γ[q]) * MT[p,q]
    @einsum RFS[n] = -(γ[n]*γ[p]*γ[q] - δ(n,p)*γ[q]) * MT[p,q]

    stf_ne = STF_CONV(t, t0, fdom, r/vp, r/vs)
    stf_ip = STF(t, t0 + r/vp, fdom)
    stf_is = STF(t, t0 + r/vs, fdom)
    stf_fp = STF_D1(t, t0 + r/vp, fdom)
    stf_fs = STF_D1(t, t0 + r/vs, fdom)

    u = zeros(3, nt)
    for k in 1:3
        u[k, :] .= FNE*RNE[k]*stf_ne .+ FIP*RIP[k]*stf_ip .+
                   FIS*RIS[k]*stf_is .+ FFP*RFP[k]*stf_fp .+ FFS*RFS[k]*stf_fs
    end
    return u
end

function compute_analytical_solution(fdsg, M0, t0, vp, vs, rho)
    @unpack t, dt, nt = fdsg.time
    @unpack fdom = fdsg.source
    @unpack Mxx, Mxy, Mxz, Myy, Myz, Mzz = fdsg.source

    dx = fdsg.domain.coordinates[1][2] - fdsg.domain.coordinates[1][1]
    dy = fdsg.domain.coordinates[2][2] - fdsg.domain.coordinates[2][1]
    dz = fdsg.domain.coordinates[3][2] - fdsg.domain.coordinates[3][1]

    ξ  = Float64[fdsg.source.x, fdsg.source.y, fdsg.source.z]
    MT = Float64[Mxx Mxy Mxz; Mxy Myy Myz; Mxz Myz Mzz]

    u = zeros(Float64, size(fdsg.geophones.data))
    v = zeros(Float64, size(fdsg.geophones.data))
    n_rec = size(u, 1)

    # staggered grid position of receivers
    rcv_pos = [
        hcat(fdsg.geophones.coords[1, :] .+ dx/2, fdsg.geophones.coords[2, :],        fdsg.geophones.coords[3, :]),
        hcat(fdsg.geophones.coords[1, :],          fdsg.geophones.coords[2, :] .+ dy/2, fdsg.geophones.coords[3, :]),
        hcat(fdsg.geophones.coords[1, :],          fdsg.geophones.coords[2, :],        fdsg.geophones.coords[3, :] .+ dz/2),
    ]

    for i in 1:n_rec
        for comp in 1:3
            pos   = Float64[rcv_pos[comp][i, 1], rcv_pos[comp][i, 2], rcv_pos[comp][i, 3]]
            u_all = analytical_displacement(pos, t, t0, fdom, M0, ξ, MT, vp, vs, rho)
            u[i, comp, :] .= u_all[comp, :]
        end
        for k in 1:3, ti in 5:nt-4
            v[i, k, ti] = ( 1/280*u[i,k,ti-4] - 4/105*u[i,k,ti-3] + 1/5*u[i,k,ti-2] - 4/5*u[i,k,ti-1]
                           + 4/5*u[i,k,ti+1] - 1/5*u[i,k,ti+2] + 4/105*u[i,k,ti+3] - 1/280*u[i,k,ti+4]) / dt
        end
    end
    return u, v
end

# compute analytical solution
M0  = config["source"]["seismic_moment"]
t0  = config["source"]["wavelet_center"]
vp  = Float64(velmod[4, 1, 1, 1])
vs  = Float64(velmod[5, 1, 1, 1])
rho = Float64(velmod[6, 1, 1, 1])

println("  Computing analytical solution...")
_, v_dc = compute_analytical_solution(sim_dc, M0, t0, vp, vs, rho)

# metrics
function cc_and_shift(a, b, dt)
    cc = [sum(a .* circshift(b, s)) for s in -(length(a)÷4):(length(a)÷4)]
    cc ./= (norm(a) * norm(b))
    idx = argmax(cc)
    shift = (idx - length(a)÷4 - 1) * dt
    return maximum(cc), shift
end
amp_ratio_dB(a, b) = 20 * log10(maximum(abs, b) / maximum(abs, a))

t_vec = collect(sim_dc.time.t)
tlims = (0.25, 0.7)
tids  = findall(t_vec .>= tlims[1] .&& t_vec .<= tlims[2])
Δt    = sim_dc.time.dt

ana = v_dc
num = sim_dc.geophones.data

val_cc = Float64[]
val_ar = Float64[]
for comp in 1:3
    cc, _ = cc_and_shift(vec(ana[:, comp, tids]), vec(num[:, comp, tids]), Δt)
    push!(val_cc, cc)
    push!(val_ar, amp_ratio_dB(ana[:, comp, tids], num[:, comp, tids]))
end

println("  CC values vx, vy, vz:  ", round.(val_cc, digits=4))
println("  AR values vx, vy, vz:  ", round.(val_ar, digits=4), " dB")

# tests (thresholds for h = 10)
@testset "DC source validation (h=10m)" begin
    @test all(val_cc .> 0.98)
    @test all(val_ar .> -0.15)
end