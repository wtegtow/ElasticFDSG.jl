include(joinpath(@__DIR__, "../../../src/shared/wavelets.jl"))
include(joinpath(@__DIR__,"create_velmod.jl"))

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

# euclidian 3d distance
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

function relative_error(a,b)
    error_rel = mean(isapprox.(a, b, atol=std(a)/1) .== false)
    return error_rel 
end 

# main 
# ----------------------------

function test_3d(abs_error_tol = 0.05)

    """
    This function compares analytically derived seismograms with numerically derived seismograms 
    for point sources acting only on vx, vy, or vz. Because 3D models are computationally too expensive 
    for tests, the function randomly selects one of the three scenarios.
    
    The test is considered passed if the relative error between the numerical 
    seismograms and the analytical solutions for all geophones is below the specified 
    absolute error tolerance (default: 5%).

    """;

    error_tolerance_not_exceeded = true
    components = ["on_vx", "on_vy", "on_vz"]

    # select a random component for the tests
    comp = rand(components)

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

    images = []

    for geo_pair in 1:3
        # geophone positions (velocity measurements)
        geo_coords = content["geophones_coords"][:, geo_pair]
        geo_x, geo_y, geo_z = geo_coords

        # geophone data 
        """
        Note: The analytical solution does not take into account the radiation pattern of the source.
            For this reason, the maximum energy of each wave type is assumed to be radiated to the receiver.
            To account for this, two seismograms are added here: one aligned along the force direction towards 
            the source and receiving the maximum P-wave energy, and another at the same distance but
            perpendicular to the source along the force direction to receive the maximum S-wave energy.

            See source-receiver geometries in the .yaml files.
        """;
        
        geo_vx = content["geophones_data"][geo_pair,1,:] + content["geophones_data"][geo_pair+3,1,:]
        geo_vy = content["geophones_data"][geo_pair,2,:] + content["geophones_data"][geo_pair+3,2,:]  
        geo_vz = content["geophones_data"][geo_pair,3,:] + content["geophones_data"][geo_pair+3,3,:]  
        # select corresponding seismogram
        geo_v = nothing
        if comp == "on_vx"
            geo_v = geo_vx
        elseif comp == "on_vy"
            geo_v = geo_vy
        elseif comp == "on_vz"
            geo_v = geo_vz
        end

        # analytical solution 
        # -----------------------------------
        G = zeros(nt) # Green's function
        v = zeros(nt) # analytical velocity
        r = dist_3d(geo_x, source_x, geo_y, source_y, geo_z, source_z)

        # p wave 
        λp = 1 / ( 4 * π * rho0 * vp0^2 * r)
        t_p = r / vp0 
        p_arr_id = round(Int, t_p / dt)
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

        # calculate relative error 
        re = relative_error(v, geo_v)

        if re > abs_error_tol
            error_tolerance_not_exceeded = false
        end

        @test error_tolerance_not_exceeded

        #begin 
        #    fig = Figure()
        #    ax = Axis(fig[1,1])
        #    lines!(ax, time, v, color="red", linestyle=:dash)
        #    lines!(ax, content["time"], geo_v, color=:blue)
        #    display(fig)
        #    image_infos = [time, v, geo_v, r]
        #    push!(images, image_infos)
        #end

    end

end

test_3d()