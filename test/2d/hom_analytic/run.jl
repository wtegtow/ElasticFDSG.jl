#using ElasticFDSG, CairoMakie, HDF5, YAML, Test, Statistics

# Utils 
# ----------------------------------------

# load wavelets
include(joinpath(@__DIR__, "../../../src/shared/wavelets.jl"))

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
dist_2d(x1,x2,y1,y2) = sqrt((x1-x2)^2 + (y1-y2)^2)

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
    
    return result ./ 2
end;

function relative_error(a,b)
    error_rel = mean(isapprox.(a, b, atol=std(a)/1) .== false)
    return error_rel 
end 

# analytical solutions
function analytical_solution_2d(content, config_path, vp, vs, rho, geo_num, mode)

    # extract informations 
    # -----------------------------------

    config = YAML.load_file(config_path)

    # simulated time
    time = content["time"]
    nt = length(time)
    dt = time[2] - time[1]

    # source locations
    source_x = config["source"]["location"]["x"]
    source_y = config["source"]["location"]["y"]

    # reconstruct source time function
    wavelet_type =  config["source"]["wavelet_type"]
    wavelet_center = config["source"]["wavelet_center"]
    fdom = config["source"]["dominant_frequency"]
    STF = wavelet_type == "ricker" ? ricker(time, wavelet_center, fdom) : gauss1d(time, wavelet_center, fdom)
    STF .*= config["source"]["amplitude"]

    # geophone positions (velocity measurements)
    geo_coords = content["geophones_coords"][:, geo_num]
    geo_x, geo_y = geo_coords

    # geophone data (only one component needed)
    geo_vx = content["geophones_data"][geo_num,1,:]
    geo_vy = content["geophones_data"][geo_num,2,:]

    geo_v = nothing
    if mode == "on_x"
        geo_v = geo_vx
    elseif mode == "on_y"
        geo_v = geo_vy
    end

    # Greens Function 
    # -----------------------------------
    G = zeros(nt) 
    v = zeros(nt) # analytical velocity
    r = dist_2d(geo_x, source_x, geo_y, source_y)

    # p wave 
    for it in 1:nt 
        if (time[it] - r / vp) >= 0 
            G[it] = 1 / ( 2 * π * rho * vp^2) * (1 / sqrt(time[it]^2 - (r/vp)^2))
        end 
    end 

    # s wave 
    for it in 1:nt 
        if (time[it] - r / vs) >= 0 
            G[it] = 1 / ( 2 * π * rho * vs^2) * (1 / sqrt(time[it]^2 - (r/vs)^2))
        end 
    end 

    # convolve for displacement
    u = convolve(G, STF .* dt)
    # derivative for velocity
    for i in 2:nt-1 
        v[i] = (u[i+1] - u[i-1]) / (2*dt)
    end

    images = []
    push!(images, [time, geo_v, v])

    #begin     
    #    fig = Figure() 
    #    ax = Axis(fig[1,1], title="$geo_num " * mode, xlabel="time [sec]")
    #    num = lines!(ax, time, geo_v, color="blue", label="Numerical")
    #    ana = lines!(ax, time, v, color="red", linestyle=:dash, label="Analytical")
    #    axislegend(ax)
    #    display(fig)
    #end

    error_rel = relative_error(v, geo_v)
    
    return error_rel, images

end


# ----------------------------------------
# Test
# ----------------------------------------

function test_2Dhom(abs_error_tol = 0.05)

    """
    Runs two 2D simulations for a homogeneous isotropic velocity model to validate 
    the accuracy of numerical seismograms. The simulations apply forces to the 
    velocity components vx and vy separately.

    The configuration file places 8 geophones around the source location, and 
    their recorded responses are compared with analytical solutions for both components 
    (vx and vy).

    The test is considered passed if the relative error between the numerical 
    seismograms and the analytical solutions for all geophones is below the specified 
    absolute error tolerance (default: 5%).

    """

    error_tolerance_not_exceeded = true

    # preconfigured velocity model file and config files
    velmod_path = joinpath(@__DIR__, "velmod.jld2")
    cfpath_on_vx = joinpath(@__DIR__, "on_vx.yaml")
    cfpath_on_vy = joinpath(@__DIR__, "on_vy.yaml")
    
    # run both simulations
    ElasticFDSG.dim2.runsim(cfpath_on_vx, velmod_path)
    ElasticFDSG.dim2.runsim(cfpath_on_vy, velmod_path)

    # extract saved results
    on_vx = joinpath(@__DIR__, "on_vx.h5");
    on_vy = joinpath(@__DIR__, "on_vy.h5");

    content_on_vx = extract_hdf5_content(on_vx);
    content_on_vy = extract_hdf5_content(on_vy);

    # from velocity model
    vp = 5000
    vs = 2500 
    rho = 2800

    # compare with analytical solutions
    for geo_num in 1:9
        error_on_vx, _ = analytical_solution_2d(content_on_vx, cfpath_on_vx, vp, vs, rho, geo_num, "on_x")
        error_on_vy, _ = analytical_solution_2d(content_on_vy, cfpath_on_vy, vp, vs, rho, geo_num, "on_y")
        
        # if any error exceeds tolerance, trigger tests
        if error_on_vx > abs_error_tol 
           println("Test failed for vx-geophone $geo_num. Error of $error_on_vx exceeded tolerance of $abs_error_tol")
           error_tolerance_not_exceeded = false
        end

        if error_on_vy > abs_error_tol 
            println("Test failed for vy-geophone $geo_num. Error of $error_on_vy exceeded tolerance of $abs_error_tol")
            error_tolerance_not_exceeded = false
         end

        @test error_tolerance_not_exceeded

    end

end;

test_2Dhom()
