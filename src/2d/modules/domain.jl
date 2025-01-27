struct Domain{T, U, V}

    velmod::T

    x0::U 
    dx::U
    xend::U 
    nx::Int64
    xcoords::V

    y0::U 
    dy::U
    yend::U 
    ny::Int64
    ycoords::V 

    dim::Tuple{Int64, Int64} 

    inner_id::Tuple{Vector{Int64}, Vector{Int64}} 

    pml_id::Tuple{Vector{Int64}, Vector{Int64}} 
    pml_points::Vector{Tuple{Int64, Int64}}
    pml_points_hash::Matrix{Int64}

    memory::Float64
    number_of_nodes::Int64 

end

# utils 
function load_velmod(VELMODPATH)
    if endswith(VELMODPATH,"npz")
        velmodfile = npzread(joinpath(VELMODPATH))
        velmod = velmodfile["velmod"];
        elseif endswith(VELMODPATH,"npy")
        velmod = npzread(joinpath(VELMODPATH))
        elseif endswith(VELMODPATH,"jld2")
        velmodfile = jldopen(joinpath(VELMODPATH))
        velmod = velmodfile["velmod"];
        close(velmodfile)
        else
        error("no valid velocity model file found.")
    end;
    return velmod 
end 

function convert_to_array_hashmap(pml_hash_dict, ny, nx)
    hash_table = fill(-1, ny, nx)
    for ((y, x), idx) in pml_hash_dict
        hash_table[y, x] = idx
    end
    return hash_table
end;

# init domain 
function init_domain(VELMODPATH, settings::Settings)
    
    N = settings.N
    nlayer_pml = settings.config["pml"]["nlayer"]

    # load velocity model from file
    velmod = settings.float.(load_velmod(VELMODPATH))

    # process 
    x0 = velmod[1,begin,begin]
    xend = velmod[1,end,end]
    y0 = velmod[2,begin,begin]
    yend = velmod[2,end,end]

    file_ny, file_nx = size(velmod[1,:,:])

    dx = abs(velmod[1,begin,begin] - velmod[1,begin,begin+1])
    dy = abs(velmod[2,begin,begin] - velmod[2,begin+1,begin])

    xcoords_inner = LinRange(x0, xend, file_nx)
    @assert length(xcoords_inner) == file_nx "Lin  X range failed"

    ycoords_inner = LinRange(y0, yend, file_ny)
    @assert length(ycoords_inner) == file_ny "Lin Y range failed"

    RIGHT = settings.config["boundaries"]["right"] == "absorbing" ? true : false
    LEFT  = settings.config["boundaries"]["left"] == "absorbing" ? true : false
    TOP   = settings.config["boundaries"]["top"] == "absorbing" ? true : false
    BOTTOM = settings.config["boundaries"]["bottom"] == "absorbing" ? true : false

    xleft  = LEFT   == true ? (N+nlayer_pml)*dx  : N*dx 
    xright = RIGHT  == true ? (N+nlayer_pml)*dx  : N*dx 
    ytop   = TOP    == true ? (N+nlayer_pml)*dy  : N*dy
    ybot   = BOTTOM == true ? (N+nlayer_pml)*dy  : N*dy

    nleft = round(Int, xleft/dx)  
    nright  = round(Int,xright/dx)
    ntop  = round(Int,ytop/dy)
    nbot = round(Int,ybot/dy)

    # expand coordinates
    xcoords = LinRange(xcoords_inner[1] - xleft, xcoords_inner[end] + xright, file_nx + nleft + nright)
    ycoords = LinRange(ycoords_inner[1] - ytop, ycoords_inner[end] + ybot, file_ny + ntop + nbot)
    nx = length(xcoords)
    ny = length(ycoords)
    dim = (ny, nx)

    # find indices of inner_domain
    inner_domain = (ycoords_inner, xcoords_inner)
    inner_xid = [argmin(abs.(x .- xcoords)) for x in xcoords_inner]
    inner_yid = [argmin(abs.(y .- ycoords)) for y in ycoords_inner]
    inner_id = (inner_yid, inner_xid)

    @assert length(inner_id[1]) == file_ny "error in domain.jl"
    @assert length(inner_id[2]) == file_nx "error in domain.jl"

    # find indices of 1d pml profiles
    pml_xid = [inner_xid[begin] - nlayer_pml : inner_xid[1] - 1 ;
               inner_xid[end] + 1 : inner_xid[end] + nlayer_pml]

    pml_yid = [inner_yid[begin] - nlayer_pml : inner_yid[1] - 1 ;
               inner_yid[end] + 1 : inner_yid[end] + nlayer_pml]
    
    pml_id = (pml_yid, pml_xid)

    # calculate all pml points (indices)
    yline = collect(pml_yid[begin]:1:pml_yid[end]) 
    xline = collect(pml_xid[begin]:1:pml_xid[end])
    xpoints = [(n, m) for n in yline, m in pml_xid]
    ypoints = [(n, m) for n in pml_yid, m in xline]
    pml_points = union(xpoints, ypoints) 
    # create hashmap for pml points 
    pml_points_dict = Dict((y, x) => i for (i, (y, x)) in enumerate(pml_points))
    pml_points_hash = convert_to_array_hashmap(pml_points_dict, ny, nx)
 
    # memory approximation
    number_of_nodes = nx*ny
    element_size = Base.sizeof(settings.float)
    memory = ((number_of_nodes * element_size) / (1024^3)) * 17 # 16 fields + 1 hashmap array
    
    # types 
    T = typeof(velmod)
    U = typeof(x0)
    V = typeof(xcoords)

    domain = Domain{T, U, V}(
                    velmod,
                    x0,dx,xend,nx,xcoords,
                    y0,dy,yend,ny,ycoords,
                    dim,
                    inner_id, 
                    pml_id, pml_points, pml_points_hash,
                    memory, number_of_nodes)

    return domain 

end;