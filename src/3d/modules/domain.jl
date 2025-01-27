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

    z0::U 
    dz::U
    zend::U 
    nz::Int64
    zcoords::V 

    dim::Tuple{Int64, Int64, Int64} 

    inner_id::Tuple{Vector{Int64}, Vector{Int64}, Vector{Int64}}

    pml_id::Vector{Vector{Int64}}
    pml_points::Vector{Any}
    pml_points_hash::Array{Int64, 3}

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


# init domain 
function init_domain(VELMODPATH, settings::Settings)
    
    N = settings.N
    nlayer_pml = settings.config["pml"]["nlayer"]

    # load velocity model from file
    velmod = settings.float.(load_velmod(VELMODPATH))

    # read coordinates 
    x0 = velmod[1,begin,begin,begin]
    xend = velmod[1,end,end,end]
    y0 = velmod[2,begin,begin,begin]
    yend = velmod[2,end,end,end]
    z0 = velmod[3,begin,begin,begin]
    zend = velmod[3,end,end,end]

    # due to name collision, i rename variables above
    x_inner_start = x0
    x_inner_end   = xend 
    y_inner_start = y0
    y_inner_end = yend 
    z_inner_start = z0 
    z_inner_end = zend


    file_nx, file_ny, file_nz = size(velmod[1,:,:,:])

    dx = abs(velmod[1,begin,begin,begin] - velmod[1,begin+1,begin,begin])
    dy = abs(velmod[2,begin,begin,begin] - velmod[2,begin,begin+1,begin])
    dz = abs(velmod[3,begin,begin,begin] - velmod[3,begin,begin,begin+1])

    dx = settings.float(dx)
    dy = settings.float(dz)
    dz = settings.float(dz)

    # inner (physical) domain
    xcoords_inner = LinRange(x0, xend, file_nx)
    @assert length(xcoords_inner) == file_nx "Ups. Something went wrong with reading the xcoords"

    ycoords_inner = LinRange(y0, yend, file_ny)
    @assert length(ycoords_inner) == file_ny "Ups. Something went wrong with reading the ycoords"

    zcoords_inner = LinRange(z0, zend, file_nz)
    @assert length(zcoords_inner) == file_nz "Ups. Something went wrong with reading the zcoords"

    # boundaries 
    xstart = settings.config["boundaries"]["xstart"] == "absorbing" ? true : false
    xend   = settings.config["boundaries"]["xend"]   == "absorbing" ? true : false
    ystart = settings.config["boundaries"]["ystart"] == "absorbing" ? true : false
    yend   = settings.config["boundaries"]["yend"]   == "absorbing" ? true : false
    zstart = settings.config["boundaries"]["zstart"] == "absorbing" ? true : false
    zend   = settings.config["boundaries"]["zend"]   == "absorbing" ? true : false

    # pml points + ghost nodes (derivative points)
    nx_start = xstart == true ? nlayer_pml + N : N
    nx_end   = xend   == true ? nlayer_pml + N : N
    ny_start = ystart == true ? nlayer_pml + N : N
    ny_end   = yend   == true ? nlayer_pml + N : N
    nz_start = zstart == true ? nlayer_pml + N : N
    nz_end   = zend   == true ? nlayer_pml + N : N

    # expand coordinates
    xcoords = LinRange(xcoords_inner[1]   - nx_start * dx, 
                       xcoords_inner[end] + nx_end   * dx,
                       file_nx + nx_start + nx_end)
    
    ycoords = LinRange(ycoords_inner[1]   - ny_start * dy, 
                       ycoords_inner[end] + ny_end   * dy,
                       file_ny + ny_start + ny_end)
    
    zcoords = LinRange(zcoords_inner[1]   - nz_start * dz,
                       zcoords_inner[end] + nz_end   * dz,
                       file_nz + nz_start + nz_end)

    nx = length(xcoords)
    ny = length(ycoords)
    nz = length(zcoords)
    dim = (nz, ny, nx)

    # find indices of inner domain
    inner_xid = [argmin(abs.(x .- xcoords)) for x in xcoords_inner]
    inner_yid = [argmin(abs.(y .- ycoords)) for y in ycoords_inner]
    inner_zid = [argmin(abs.(z .- zcoords)) for z in zcoords_inner]

    inner_id = (inner_zid, inner_yid, inner_xid) # here order is changed to z,y,x 
    @assert length(inner_id[1]) == file_nz "error in domain.jl"
    @assert length(inner_id[2]) == file_ny "error in domain.jl"
    @assert length(inner_id[3]) == file_nx "error in domain.jl"

    # calculate indices of 1d pml profiles
    pml_xid = [xstart == true ? (N+1:N+nlayer_pml) : nothing ;
               xend   == true ? (inner_xid[end]+1:inner_xid[end]+nlayer_pml) : nothing]

    pml_yid = [ystart == true ? (N+1:N+nlayer_pml) : nothing ;
               yend   == true ? (inner_yid[end]+1:inner_yid[end]+nlayer_pml) : nothing]
        
    pml_zid = [zstart == true ? (N+1:N+nlayer_pml) : nothing ;
               zend   == true ? (inner_zid[end]+1:inner_zid[end]+nlayer_pml) : nothing]

    pml_id = [pml_zid, pml_yid, pml_xid]

    # calculate pml points (indices) and create array hashmap
    pml_points = []
    pml_points_hash = zeros(Int, dim);
    pml_points_hash .= -1

    pml_index = 1
    for z in 1:nz, y in 1:ny, x in 1:nx 
        
        if z in pml_zid ||
           y in pml_yid || 
           x in pml_xid 

           push!(pml_points, (z,y,x))
           pml_points_hash[z,y,x] = pml_index
           pml_index += 1

        end

    end
    
    pml_points = union(pml_points)
    @assert length(pml_points) == pml_index - 1 # because the index started with 1
 
    # memory approximation
    number_of_nodes = prod(dim)
    element_size = Base.sizeof(settings.float)
    memory = ((number_of_nodes * element_size) / (1024^3)) * 20 
    
    # types 
    T = typeof(velmod)
    U = typeof(x0)
    V = typeof(xcoords)

    domain = Domain{T, U, V}(
                    velmod,
                    x_inner_start,dx,x_inner_end,nx,xcoords,
                    y_inner_start,dy,y_inner_end,ny,ycoords,
                    z_inner_start,dz,z_inner_end,nz,zcoords,
                    dim,
                    inner_id, 
                    pml_id, pml_points, pml_points_hash,
                    memory, number_of_nodes)

    return domain 

end;