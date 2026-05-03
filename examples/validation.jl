#= 
This script validates the numerical solution of ElasticFDSG against the analytical solution for a point moment tensor source in a homogeneous medium. It runs three simulations with different moment tensor types (ISO, DC, CLVD) and compares the seismograms at the receiver locations to the analytical solution using cross-correlation and amplitude ratio metrics.
=#


# NOTE: due to the project structure, i must switch back and forth to import dependencies 
using Pkg; Pkg.activate(joinpath(@__DIR__, ".."))
using ElasticFDSG
ElasticFDSG.devmode!(true)
using Test, LinearAlgebra, Statistics

using Pkg; Pkg.activate(joinpath(@__DIR__, "..", "test"))
using JLD2, YAML, Einsum, UnPack

using Pkg; Pkg.activate(joinpath(@__DIR__, "..", ".dev"))
using GLMakie, LaTeXStrings, Metal

# velocity model
h = 5.0;
xcoords = 0:h:1000;
ycoords = 0:h:300;
zcoords = 0:h:1000;
nx, ny, nz = length(xcoords), length(ycoords), length(zcoords);

X = getindex.(Iterators.product(xcoords, ycoords, zcoords), 1);
Y = getindex.(Iterators.product(xcoords, ycoords, zcoords), 2);
Z = getindex.(Iterators.product(xcoords, ycoords, zcoords), 3);

velmod = zeros(Float32, 13, nx, ny ,nz);
velmod[1,:,:,:] .= X;
velmod[2,:,:,:] .= Y;
velmod[3,:,:,:] .= Z;
velmod[4,:,:,:] .= 3000;
velmod[5,:,:,:] .= 2000;
velmod[6,:,:,:] .= 2500;

# receiver array
src_x, src_y, src_z = (50, 100, 500)
rcv = [Dict("x" => 950, "y"=>200, "z"=>z) for z in 50:100:950]
rcv_x, rcv_y, rcv_z = ([x["x"] for x in rcv], [y["y"] for y in rcv], [z["z"] for z in rcv])

# config
config = Dict(
    "settings" => Dict(
        "device" => "metal",
        "precision" => "Float32",            
        "spatial_derivative_order" => 4,      
        "verbose" => false,                   
        "output_file" => nothing
    ),
    "time" => Dict(
        "start" => 0.0,
        "end" => 0.7,
        "timestep" => 0.001
    ),
    "source" => Dict(
        "dominant_frequency" => 20,
        "wavelet_type" => "ricker",          
        "wavelet_center" => 0.0625,
        "location" => Dict(
            "x" => src_x,
            "y" => src_y,
            "z" => src_z
        ),
        "seismic_moment" => 1e15,
        "moment_tensor" => Dict(
            "Mxx" => 0,
            "Myy" => 0,
            "Mzz" => 0,
            "Mxy" => 0,
            "Mxz" => 0,
            "Myz" => 0.,
            "anisotropic" => false
        )
    ),
    "boundaries" => Dict(
        "xstart" => "absorbing",
        "xend" => "absorbing",
        "ystart" => "absorbing",
        "yend" => "absorbing",
        "zstart" => "absorbing",
        "zend" => "absorbing",
        "pml_layer" => 10
    ),
    "receivers" => Dict(
        "geophones" => rcv,
        "das" => Dict(
            "x_aligned" => [],  
            "y_aligned" => [],  
            "z_aligned" => []  
        ),
        "snapshots" => Dict(
            "plane_positions" => [
                Dict("x" => src_x, "y" => src_y, "z" => src_z)
            ],
            "times" => collect(0:0.05:0.75),           
            "fields" => ["vx", "vy", "vz"]  
        )
    )
)

# moment tensors
M_iso  = normalize([1 0 0; 0 1 0; 0 0 1])
M_dc   = normalize([1 0 0; 0 0 0; 0 0 -1])
M_clvd = normalize([1 0 0; 0 1 0; 0 0 -2])

# fill config for each source type
config_iso  = deepcopy(config)
config_dc   = deepcopy(config)
config_clvd = deepcopy(config)

mt_indices = Dict("Mxx"=>(1,1), "Myy"=>(2,2), "Mzz"=>(3,3),
                  "Mxy"=>(1,2), "Mxz"=>(1,3), "Myz"=>(2,3))
for (key, (i,j)) in mt_indices
    config_iso["source"]["moment_tensor"][key]  = M_iso[i,j]
    config_dc["source"]["moment_tensor"][key]   = M_dc[i,j]
    config_clvd["source"]["moment_tensor"][key] = M_clvd[i,j]
end
configs = (config_iso, config_dc, config_clvd);

# run simulations
mtensor_types = ["ISO", "DC", "CLVD"]
Simulations = []
for i in 1:3
    println("Running simulation for $(mtensor_types[i]) source...")
    sim = ElasticFDSG.runsim(configs[i], velmod)
    push!(Simulations, sim)
end
sim_iso, sim_dc, sim_clvd = Simulations;

# analytical solution
STF(t, t0, f0)    = @. (1 - 2*(π*f0)^2*(t-t0)^2) * exp(-(π*f0)^2*(t-t0)^2);
STF_D1(t, t0, f0) = @. -2*(π*f0)^2*(t-t0)*(3 - 2*(π*f0*(t-t0))^2) * exp(-(π*f0*(t-t0))^2);

function STF_CONV(t, t0, f0, a, b)
    Δt = t[2] - t[1]
    nt = length(t)
    res = zeros(nt)
    for i in 1:nt, τ in a:Δt:b
        res[i] += (1 - 2*(π*f0)^2*(t[i]-t0-τ)^2) * exp(-(π*f0)^2*(t[i]-t0-τ)^2) * Δt * τ
    end
    return res
end;

function analytical_displacement(pos, t, t0, fdom, M0, ξ, MT, vp, vs, rho)
    δ(i,j) = i == j ? 1 : 0
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
    @einsum RIS[n] = -(6*γ[n]*γ[p]*γ[q] -  γ[n]*δ(p,q) - γ[p]*δ(n,q) - 2*γ[q]*δ(n,p)) * MT[p,q]
    @einsum RFP[n] = (γ[n]*γ[p]*γ[q]) * MT[p,q]
    @einsum RFS[n] = -(γ[n]*γ[p]*γ[q] - δ(n,p)*γ[q]) * MT[p,q]

    stf_ne = STF_CONV(t, t0, fdom, r/vp, r/vs)
    stf_ip = STF(t, t0 + r/vp, fdom)
    stf_is = STF(t, t0 + r/vs, fdom)
    stf_fp = STF_D1(t, t0 + r/vp, fdom)
    stf_fs = STF_D1(t, t0 + r/vs, fdom)

    u = zeros(3, nt)
    for k in 1:3
        u[k,:] .= FNE*RNE[k]*stf_ne .+ FIP*RIP[k]*stf_ip .+
                  FIS*RIS[k]*stf_is .+ FFP*RFP[k]*stf_fp .+ FFS*RFS[k]*stf_fs
    end
    return u
end;

function compute_analytical_solution(fdsg, M0, t0, vp, vs, rho)

    @unpack t, dt, nt = fdsg.time
    @unpack fdom = fdsg.source
    @unpack Mxx, Mxy, Mxz, Myy, Myz, Mzz = fdsg.source

    dx = fdsg.domain.coordinates[1][2] - fdsg.domain.coordinates[1][1]
    dy = fdsg.domain.coordinates[2][2] - fdsg.domain.coordinates[2][1]
    dz = fdsg.domain.coordinates[3][2] - fdsg.domain.coordinates[3][1]

    # Source is injected into stress divergence at the main grid node (normal stresses)
    # → no staggered offset for ξ
    ξ  = Float64[fdsg.source.x, fdsg.source.y, fdsg.source.z]
    MT = Float64[Mxx Mxy Mxz; Mxy Myy Myz; Mxz Myz Mzz]

    u = zeros(Float64, size(fdsg.geophones.data))
    v = zeros(Float64, size(fdsg.geophones.data))

    n_rec = size(u, 1)

    # Each velocity component lives at a different staggered sub-node.
    # Evaluate the analytical solution at the correct sub-node per component.
    rcv_pos = [  # (n_rec × 3) position matrix per component
        hcat(fdsg.geophones.coords[1,:] .+ dx/2, fdsg.geophones.coords[2,:],        fdsg.geophones.coords[3,:]),          # vx
        hcat(fdsg.geophones.coords[1,:],          fdsg.geophones.coords[2,:] .+ dy/2, fdsg.geophones.coords[3,:]),        # vy
        hcat(fdsg.geophones.coords[1,:],          fdsg.geophones.coords[2,:],        fdsg.geophones.coords[3,:] .+ dz/2), # vz
    ]

    for i in 1:n_rec
        for comp in 1:3
            pos   = Float64[rcv_pos[comp][i,1], rcv_pos[comp][i,2], rcv_pos[comp][i,3]]
            u_all = analytical_displacement(pos, t, t0, fdom, M0, ξ, MT, vp, vs, rho)
            u[i, comp, :] .= u_all[comp, :]   # take only the matching component at this position
        end
        # 8th order central finite difference for velocity
        for k in 1:3, t in 5:nt-4
            v[i,k,t] = ( 1/280*u[i,k,t-4] - 4/105*u[i,k,t-3] + 1/5*u[i,k,t-2] - 4/5*u[i,k,t-1]
                        + 4/5*u[i,k,t+1] - 1/5*u[i,k,t+2] + 4/105*u[i,k,t+3] - 1/280*u[i,k,t+4]) / dt
        end
    end
    return u,v
end;


M0  = config["source"]["seismic_moment"]
t0  = config["source"]["wavelet_center"]
vp  = Float64(velmod[4,1,1,1])
vs  = Float64(velmod[5,1,1,1])
rho = Float64(velmod[6,1,1,1])

u_iso, v_iso   = compute_analytical_solution(sim_iso,  M0, t0, vp, vs, rho);
u_dc, v_dc     = compute_analytical_solution(sim_dc,   M0, t0, vp, vs, rho);
u_clvd, v_clvd = compute_analytical_solution(sim_clvd, M0, t0, vp, vs, rho);

analytics = [v_iso, v_dc, v_clvd];
numerics  = [sim_iso.geophones.data, sim_dc.geophones.data, sim_clvd.geophones.data];

# metrics
begin
    function cc_and_shift(a, b, dt)
        cc = [sum(a .* circshift(b, s)) for s in -(length(a)÷4):(length(a)÷4)]
        cc ./= (norm(a) * norm(b))
        idx = argmax(cc)
        shift = (idx - length(a)÷4 - 1) * dt
        return maximum(cc), shift
    end
    amp_ratio_dB(a, b) = 20 * log10(maximum(abs, b) / maximum(abs, a))

    t = collect(sim_iso.time.t)
    tlims = (0.25, 0.7) # only compare non-zero region of the seismograms
    tids = findall(t .>= tlims[1] .&& t .<= tlims[2])
    Δt = sim_iso.time.dt

    val_cc = []
    val_ar = [] 
    for i in 1:3
        ana = analytics[i]  
        num = numerics[i]

        for comp in 1:3
            cc, shift = cc_and_shift(vec(ana[:, comp, tids]), vec(num[:, comp, tids]), Δt)
            push!(val_cc, cc)        
            amp_ratio = amp_ratio_dB(ana[:, comp, tids], num[:, comp, tids])
            push!(val_ar, amp_ratio)
        end
    end
end

# tests
if h == 5
    @test all(val_cc .> 0.99)
    @test all(val_ar .> -0.2)  
elseif h == 10
    @test all(val_cc .> 0.98)
    @test all(val_ar .> -0.38)  
end;


# figure
norm11(a) = a ./ maximum(abs, a);

fig = Figure(size=(1000, 1000), fontsize=18, font = :bold)
titles = ["ISO", "DC", "CLVD"]
components = ["V_x", "V_y", "V_z"]

first_num = nothing
first_ana = nothing
counter = 1
digits = 4
for comp in 1:3
    for i in 1:3
        ana = analytics[i]  
        num = numerics[i]
    
        ylab = i == 1 ? "Receiver Index" : ""
        xlab = comp == 3 ? "Time [sec]" : ""

        ax = Axis(fig[comp, i],
                title = latexstring("\$$(components[comp])\$ — $(titles[i]) \\\\ CC: $(round(val_cc[counter],digits=digits)),  AR: $(round(val_ar[counter],digits=digits)) dB"),
                yreversed = true,
                ylabel = ylab, xlabel=xlab)

        n_rcv = size(ana,1)
        for n in 1:n_rcv
            local_scale = maximum(abs, num[n, comp, :]) 
            offset = n * 2.0 
            h_num = lines!(ax, t, (num[n, comp, :]) ./ local_scale .+ offset, color = :blue)
            h_ana = lines!(ax, t .+ Δt/2, (ana[n, comp, :]) ./ local_scale .+ offset, color = :red, linestyle = :dash)

            if isnothing(first_num)
                first_num = h_num
                first_ana = h_ana
            end
        end
        xlims!(ax, tlims...)
        counter += 1
    end
end
Legend(fig[2, 4], [first_num, first_ana], ["Numerical", "Analytical"])
display(fig)
GLMakie.save(joinpath(@__DIR__, "validation.png"), fig)