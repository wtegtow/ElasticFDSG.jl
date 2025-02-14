include(joinpath(@__DIR__, "../../../src/shared/wavelets.jl"))
include(joinpath(@__DIR__,"create_velmod.jl"))

using CairoMakie, YAML, HDF5, JLD2 
CairoMakie.activate!(type = "png")
# utils 
# -----------------------------------------------

# extract content from result file 
function extract_hdf5_content(file_path::String)
    file = h5open(file_path, "r")
    result = Dict()
    
    function process_group(group, prefix="")
        for name in keys(group)
            path = joinpath(prefix, name)
            obj = group[name]
            if obj isa HDF5.Group
                process_group(obj, path)
            elseif obj isa HDF5.Dataset
                result[path] = read(obj)
            else
                ;
            end
        end
    end
    
    process_group(file)
    close(file)
    return result
end;

# euclidian 2d distance
dist_3d(x1,x2,y1,y2,z1,z2) = sqrt((x1-x2)^2 + (y1-y2)^2 + (z1-z2)^2)

# dummy convolution
function convolve(a, b)
    n = length(a)
    m = length(b)
    
    result_len = n + m - 1
    result = zeros(eltype(a), result_len)
    
    for i in 1:result_len
        for j in max(1, i+1-m):min(i, n)
            result[i] += a[j] * b[i - j + 1]
        end
    end
    
    return result 
end;


# main 
# ----------------------------

components = ["on_vx", "on_vy", "on_vz"]

# select a random component for the tests
comp = rand(components)
comp = "on_vz"

# create the corresponding velmod
vp0 = 5000
vs0 = 2500 
rho0 = 2800
create_3dtest_velmod(comp, vp0, vs0, rho0)

# run sim
path_to_velmod = joinpath(@__DIR__, "velmod.jld2")
path_to_config = joinpath(@__DIR__, "$(comp).yaml")
ElasticFDSG.dim3.runsim(path_to_config, path_to_velmod)

# extract informations 
# -----------------------------------
config = YAML.load_file(path_to_config)
content = extract_hdf5_content(joinpath(@__DIR__, "results.h5"));

# simulated time
time = content["time"]
nt = length(time)
dt = time[2] - time[1]

# source locations
source_x = config["source"]["location"]["x"]
source_y = config["source"]["location"]["y"]
source_z = config["source"]["location"]["z"]

# reconstruct source time function
wavelet_type =  config["source"]["wavelet_type"]
wavelet_center = config["source"]["wavelet_center"]
fdom = config["source"]["dominant_frequency"]
STF = wavelet_type == "ricker" ? ricker(time, wavelet_center, fdom) : gauss1d(time, wavelet_center, fdom)
STF .*= config["source"]["amplitude"]


images1 = []

for geo_pair in 1:3
    # geophone positions (velocity measurements)
    geo_coords = content["geophones_coords"][:, geo_pair]
    geo_x, geo_y, geo_z = geo_coords

    # geophone data 
    geo_vx = content["geophones_data"][geo_pair,1,:] + content["geophones_data"][geo_pair+3,1,:]
    geo_vy = content["geophones_data"][geo_pair,2,:] + content["geophones_data"][geo_pair+3,2,:]  
    geo_vz = content["geophones_data"][geo_pair,3,:] + content["geophones_data"][geo_pair+3,3,:]  
    geo_v = nothing
    if comp == "on_vx"
        geo_v = geo_vx
    elseif comp == "on_vy"
        geo_v = geo_vy
    elseif comp == "on_vz"
        geo_v = geo_vz
    end

    # Greens Function 
    # -----------------------------------
    G = zeros(nt) 
    v = zeros(nt) # analytical velocity
    r = dist_3d(geo_x, source_x, geo_y, source_y, geo_z, source_z)

    # p wave 
    λp = 1 / ( 4 * π * rho0 * vp0^2 * r)
    t_p = r / vp0 
    p_arr_id = ceil(Int, t_p / dt)
    G[p_arr_id]= 1/dt * λp

    # s wave 
    λs = 1 / ( 4 * π * rho0 * vs0^2 * r)
    t_s = r / vs0
    s_arr_id = round(Int, t_s / dt)
    G[s_arr_id]= λs * 1/dt

    # convolve for displacement
    u = convolve(G, STF .* dt) 

    # derivative for velocity
    for i in 2:nt-1 
        v[i] = (u[i+1] - u[i-1]) / (2*dt)
    end

    begin 
        fig = Figure()
        ax = Axis(fig[1,1])
        lines!(ax, time, v, color="red", linestyle=:dash)
        lines!(ax, content["time"], geo_v, color=:blue)
        display(fig)
    end

    image_infos = [time, v, geo_v, r]

    push!(images1, image_infos)

end




begin 
    h = 700 
    w = 300
    fig = Figure(size=(h,w))

    ax00 = Axis(fig[1,1:2], title="Station:1 Dist $(images1[1][4])m \n N = 1", ylabel="vz")
    ax10 = Axis(fig[2,1:2], title="N = 10",ylabel="vz")

    ax01 = Axis(fig[1,3:4], title="Station:2 Dist $(images1[2][4])m \n N = 1")
    ax11 = Axis(fig[2,3:4], title="N = 10")

    ax02 = Axis(fig[1,5:6], title="Station:3 Dist $(images1[3][4])m \n N = 1")
    ax12 = Axis(fig[2,5:6], title="N = 10")

    # ---- 1
    lines!(ax00, images1[1][1], images1[1][3], color=:blue)
    lines!(ax00, images1[1][1], images1[1][2], color="red", linestyle=:dash)
    xlims!(ax00, 0, 0.3)

    lines!(ax10, images10[1][1], images10[1][3], color=:blue)
    lines!(ax10, images10[1][1], images10[1][2], color="red", linestyle=:dash)
    xlims!(ax10, 0, 0.3)


    # ----- 2
    lines!(ax01, images1[2][1], images1[2][3], color=:blue)
    lines!(ax01, images1[2][1], images1[2][2], color="red", linestyle=:dash)
    xlims!(ax01, 0.1, 0.5)

    lines!(ax11, images10[2][1], images10[2][3], color=:blue)
    lines!(ax11, images10[2][1], images10[2][2], color="red", linestyle=:dash)
    xlims!(ax11, 0.1, 0.5)

    # ----- 3
    lines!(ax02, images1[3][1], images1[3][3], color=:blue)
    lines!(ax02, images1[3][1], images1[3][2], color="red", linestyle=:dash)
    xlims!(ax02, 0.15, 0.6)

    lines!(ax12, images10[3][1], images10[3][3], color=:blue)
    lines!(ax12, images10[3][1], images10[3][2], color="red", linestyle=:dash)
    xlims!(ax12, 0.15, 0.6)


    display(fig)

    save(joinpath(@__DIR__,"ana3d$(h)$(w).png"), fig)
end 