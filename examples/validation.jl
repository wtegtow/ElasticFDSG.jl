"""
This script computes the analytical solutions for the response of homogeneous elastic solid subjected to:
    - an isotropic (explosive) moment tensor source (ISO),
    - a double-couple moment tensor source (DC),
    - a compensated linear vector dipole (CLVD),
for a given source-receiver geometry. The analytical solutions are then compared with forward-modeled seismograms.
""";

using ElasticFDSG
using JLD2, YAML 
using GLMakie
Makie.inline!(true)
using LinearAlgebra, Einsum, UnPack, Statistics

# create velocity model 
xcoords = 0:10:1000 
ycoords = 0:10:300
zcoords = 0:10:1000 
nx, ny, nz = length(xcoords), length(ycoords), length(zcoords) 

X = getindex.(Iterators.product(xcoords, ycoords, zcoords), 1)
Y = getindex.(Iterators.product(xcoords, ycoords, zcoords), 2)
Z = getindex.(Iterators.product(xcoords, ycoords, zcoords), 3)

velmod = zeros(Float32, 13, nx, ny ,nz);
velmod[1,:,:,:] .= X;
velmod[2,:,:,:] .= Y;
velmod[3,:,:,:] .= Z;
velmod[4,:,:,:] .= 3000
velmod[5,:,:,:] .= 2000
velmod[6,:,:,:] .= 2500

VELMODPATH = joinpath(@__DIR__, "velmod.jld2")
jldsave(VELMODPATH; velmod);

# define source-receiver locations 
"""
IMPORTANT:
    Due to the staggered grid, velocity components are not defined exactly on fault planes.
    For this reason, receivers should not sit on extended fault planes.
""";

src_x, src_y, src_z = (50, 100, 500)
rcv = [Dict("x" => 950, "y"=>200, "z"=>z) for z in 50:100:950]
rcv_x, rcv_y, rcv_z = ([x["x"] for x in rcv], [y["y"] for y in rcv], [z["z"] for z in rcv])

# create configuration templates for 3 different moment tensor sources 
config = Dict(
    "settings" => Dict(
        "device" => "metal",                 
        "precision" => "Float32",            
        "spatial_derivative_order" => 4,      
        "show_progress_in_console" => true,   
        "output_file" => @__DIR__
    ),
    "time" => Dict(
        "start" => 0.0,
        "end" => 1,
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
                Dict("x" => 400, "y" => 400, "z" => 400)
            ],
            "times" => collect(0:0.05:0.75),           
            "fields" => ["vx", "vy", "vz"]  
        )
    )
)
# define moment tensors
M_iso  = normalize([1 0 0; 0 1 0; 0 0 1])
M_dc   = normalize([0 1 0; 1 0 0; 0 0 0])
M_clvd = normalize([1 0 0; 0 1 0; 0 0 -2])

# fill moment tensors to configs 
config_iso = deepcopy(config)
config_dc = deepcopy(config)
config_clvd = deepcopy(config)
indices = Dict(
    "Mxx" => (1,1),
    "Myy" => (2,2),
    "Mzz" => (3,3),
    "Mxy" => (1,2),
    "Mxz" => (1,3),
    "Myz" => (2,3)
)

for (key, (i,j)) in indices
    config_iso["source"]["moment_tensor"][key] = M_iso[i,j]   
    config_dc["source"]["moment_tensor"][key] = M_dc[i,j]   
    config_clvd["source"]["moment_tensor"][key] = M_clvd[i,j]   
end

# save all configs as .yaml files
configs= (
    (joinpath(@__DIR__, "cf_iso.yaml"), config_iso),
    (joinpath(@__DIR__, "cf_dc.yaml"), config_dc),
    (joinpath(@__DIR__, "cf_clvd.yaml"), config_clvd)
);
for (name, cf) in configs 
    open(name, "w") do file
        YAML.write(file, cf)
    end
end;

# run 3 forward simulations for each config file 
Simulations = []
for i in 1:3 
    sim = ElasticFDSG.dim3.runsim(configs[i][1], VELMODPATH; return_ = true)
    push!(Simulations, sim)
end; 
sim_iso, sim_dc, sim_clvd = Simulations;

# compute analytical solutions
function compute_analytical_solution(fdsg, M0, t0, vp, vs, rho)
    # Helper
    function STF(t, t0, f0)
    return @. (1-(2*π^2*f0^2)*(t-t0)^2)*(exp((-π^2*f0^2)*(t-t0)^2)) 
    end;
    function STF_D1(t, t0, f0)
        return @. -2*(π*f0)^2*(t-t0)*(3-2*(π*f0*(t-t0))^2)*exp(-(π*f0*(t-t0))^2)
    end;
    function STF_CONV(t, t0, f0, a, b)
        Δt = t[2]-t[1]
        nt = size(t,1)
        res = zeros(nt)
        for i in 1:nt 
            for τ in a:Δt:b 
                res[i] += (1-(2*π^2*f0^2)*(t[i]-t0-τ)^2)*(exp((-π^2*f0^2)*(t[i]-t0-τ)^2))*Δt*τ
            end 
        end
        return res
    end
    δ(i,j) = i == j ? 1 : 0

    # Main
    # unpack variables
    @unpack dx = fdsg.domain
    @unpack t, dt, nt = fdsg.time
    @unpack stf, fdom = fdsg.source
    @unpack Mxx, Mxy, Mxz, Myy, Myz, Mzz, x, y, z = fdsg.source 
    src_x = x; src_y = y; src_z = z
    M = [Mxx Mxy Mxz; Mxy Myy Myz; Mxz Myz Mzz]
    @show M

    # analytical displacement and velocities
    u = similar(fdsg.geophones.data)[:,:,:] .= 0;
    v = similar(fdsg.geophones.data)[:,:,:] .= 0;

    # receiver 
    n_rec = size(u, 1)
    rcv_x = fdsg.geophones.coords[1,:]
    rcv_y = fdsg.geophones.coords[2,:]
    rcv_z = fdsg.geophones.coords[3,:]

    # compute analytical solution for every receiver
    for i in 1:n_rec

        ξ = Float64[src_x + dx/2, src_y, src_z]
        x = Float64[rcv_x[i], rcv_y[i], rcv_z[i]]
        r = norm(x - ξ)
        γ = (x - ξ) / r

        FNE = M0 / (4π * rho * r^4)
        FIP = M0 / (4π * rho * r^2 * vp^2)
        FIS = M0 / (4π * rho * r^2 * vs^2)
        FFP = M0 / (4π * rho * r * vp^3)
        FFS = M0 / (4π * rho * r * vs^3)

        RNE = zeros(3)
        @einsum RNE[n] = (15*γ[n]*γ[p]*γ[q] - 
                            3*γ[n]*δ(p,q)  - 
                            3*γ[p]*δ(n,q)  -
                            3*γ[q]*δ(n,p)) * M[p,q]

        RIP = zeros(3)
        @einsum RIP[n] = (6*γ[n]*γ[p]*γ[q] - 
                            γ[n]*δ(p,q)  - 
                            γ[p]*δ(n,q)  -
                            γ[q]*δ(n,p)) * M[p,q]

        RIS = zeros(3)
        @einsum RIS[n] = -(6*γ[n]*γ[p]*γ[q]  - 
                            γ[n]*δ(p,q)    - 
                            γ[p]*δ(n,q)    -
                            2*γ[q]*δ(n,p)) * M[p,q]

        RFP = zeros(3)
        @einsum RFP[n] = (γ[n]*γ[p]*γ[q]) * M[p,q]

        RFS = zeros(3)
        @einsum RFS[n] =  -(γ[n]*γ[p]*γ[q] - δ(n,p) * γ[q]) * M[p,q]

        STF_NE = STF_CONV(t, t0, fdom, r/vp, r/vs)
        STF_IP = STF(t, t0 + r/vp, fdom)
        STF_IS = STF(t, t0 + r/vs, fdom)
        STF_FP = STF_D1(t, t0 + r/vp, fdom)
        STF_FS = STF_D1(t, t0 + r/vs, fdom)

        # displacements
        for k in 1:3 
            u[i,k,:] .= FNE * RNE[k] * STF_NE +
                        FIP * RIP[k] * STF_IP +
                        FIS * RIS[k] * STF_IS + 
                        FFP * RFP[k] * STF_FP +
                        FFS * RFS[k] * STF_FS
        end
        # velocities
        for k in 1:3, t in 2:nt-1 
            v[i,k,t] = (u[i,k,t+1] - u[i,k,t-1]) / (2*dt)
        end
    end;
    return u, v
end


M0 = config["source"]["seismic_moment"]
t0 = config["source"]["wavelet_center"]
vp = velmod[4,1,1,1]
vs = velmod[5,1,1,1]
rho = velmod[6,1,1,1]
u_iso, v_iso   = compute_analytical_solution(sim_iso,  M0, t0, vp, vs, rho);
u_dc, v_dc     = compute_analytical_solution(sim_dc,   M0, t0, vp, vs, rho);
u_clvd, v_clvd = compute_analytical_solution(sim_clvd, M0, t0, vp, vs, rho);

t = collect(sim_iso.time.t)
analytics = [v_iso, v_dc, v_clvd]
numerics = [sim_iso.geophones.data, sim_dc.geophones.data, sim_clvd.geophones.data]

# compare, visualize and compute errors
rmse(a, b) = sqrt(mean((a .- b).^2))
nrmse(a, b) = rmse(a, b) / (maximum(a) - minimum(a))
mae(a, b) = mean(abs.(a .- b))

begin
fig = Figure(size=(1000, 1000))
titles = ["Iso", "DC", "CLVD"]
components = ["Vx", "Vy", "Vz"]
first_num = nothing
first_ana = nothing

for comp in 1:3
    for i in 1:3
        ana = analytics[i]  
        num = numerics[i]
        
        rmerr  = nrmse(ana[:, comp, :], num[:, comp, :])
        maeerr = mae(ana[:, comp, :], num[:, comp, :])
        digits = 4

        ylab = i == 1 ? "Receiver Index" : ""
        xlab = comp == 3 ? "Time [sec]" : ""

        ax = Axis(fig[comp, i],
                  title = "$(components[comp]) — $(titles[i])\nRMSE: $(round(rmerr,digits=digits)) | MAE: $(round(maeerr,digits=digits))",
                  yreversed = true,
                  ylabel = ylab, xlabel=xlab)

        n_rcv = size(ana,1)
        for n in 1:n_rcv
            local_scale = maximum(abs, num[n, comp, :]) 
            offset = n * 2.0 
            h_num = lines!(ax, t, (num[n, comp, :]) ./ local_scale .+ offset, color = :blue)
            h_ana = lines!(ax, t, (ana[n, comp, :]) ./ local_scale .+ offset, color = :red, linestyle = :dash)

            if first_num === nothing
                first_num = h_num
                first_ana = h_ana
            end
        end
        xlims!(ax, 0.2, 0.7)
    end
end
Legend(fig[2, 4], [first_num, first_ana], ["numerical", "analytic"])
display(fig)
save(joinpath(@__DIR__, "validation.png"), fig)
end

# delete config files and velocity model to keep the repository clean
rm(VELMODPATH)
for (name, cf) in configs 
    rm(name)
end;