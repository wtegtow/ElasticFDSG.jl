struct Domain{T1,T2,T3,T4}

    x0::T1
    dx::T2
    xend::T1
    nx::Int64
    xcoords::T3

    z0::T1
    dz::T2
    zend::T1
    nz::Int64
    zcoords::T3

    dim::Tuple{Int64, Int64} 
    inner_id::Tuple{Vector{Int64}, Vector{Int64}}

    pml_id::Vector{Vector{Int64}}
    pml_lookup::T4

end


function init_domain(settings::Settings, velmod)
    
    N = settings.N
    nlayer_pml = settings.config["boundaries"]["pml_layer"]

    # read coordinates 
    x0 = velmod[1,begin,begin,begin]
    xend = velmod[1,end,end,end]
    z0 = velmod[2,begin,begin,begin]
    zend = velmod[2,end,end,end]

    # due to name collision, i rename variables above
    x_inner_start = x0
    x_inner_end   = xend 
    z_inner_start = z0 
    z_inner_end = zend

    file_nx, file_nz = size(velmod[1,:,:])

    dx = abs(velmod[1,begin,begin] - velmod[1,begin+1,begin])
    dz = abs(velmod[2,begin,begin] - velmod[2,begin,begin+1])

    dx = settings.float(dx)
    dz = settings.float(dz)

    # inner (physical) domain
    xcoords_inner = LinRange(x0, xend, file_nx)
    @assert length(xcoords_inner) == file_nx "Ups. Something went wrong with reading the xcoords"

    zcoords_inner = LinRange(z0, zend, file_nz)
    @assert length(zcoords_inner) == file_nz "Ups. Something went wrong with reading the zcoords"

    # boundaries 
    xstart = settings.config["boundaries"]["xstart"] == "absorbing" ? true : false
    xend   = settings.config["boundaries"]["xend"]   == "absorbing" ? true : false
    zstart = settings.config["boundaries"]["zstart"] == "absorbing" ? true : false
    zend   = settings.config["boundaries"]["zend"]   == "absorbing" ? true : false

    # pml points + ghost nodes (derivative points)
    nx_start = xstart == true ? nlayer_pml + N : N
    nx_end   = xend   == true ? nlayer_pml + N : N
    nz_start = zstart == true ? nlayer_pml + N : N
    nz_end   = zend   == true ? nlayer_pml + N : N

    # expand coordinates
    xcoords = LinRange(xcoords_inner[1]   - nx_start * dx, 
                       xcoords_inner[end] + nx_end   * dx,
                       file_nx + nx_start + nx_end)
    
    zcoords = LinRange(zcoords_inner[1]   - nz_start * dz,
                       zcoords_inner[end] + nz_end   * dz,
                       file_nz + nz_start + nz_end)

    nx = length(xcoords)
    nz = length(zcoords)
    dim = (nz, nx)

    # find indices of inner domain
    inner_xid = [argmin(abs.(x .- xcoords)) for x in xcoords_inner]
    inner_zid = [argmin(abs.(z .- zcoords)) for z in zcoords_inner]

    inner_id = (inner_zid, inner_xid) # here order is changed to z,y,x 
    @assert length(inner_id[1]) == file_nz "error in domain.jl"
    @assert length(inner_id[2]) == file_nx "error in domain.jl"

    # calculate indices of 1d pml profiles
    pml_xid = [xstart == true ? (N+1:N+nlayer_pml) : nothing ;
               xend   == true ? (inner_xid[end]+1:inner_xid[end]+nlayer_pml) : nothing]
        
    pml_zid = [zstart == true ? (N+1:N+nlayer_pml) : nothing ;
               zend   == true ? (inner_zid[end]+1:inner_zid[end]+nlayer_pml) : nothing]

    pml_id = [pml_zid, pml_xid]

    # calculate pml points (indices) and create array hashmap
    pml_points = []
    pml_lookup = zeros(Int, dim);
    pml_lookup .= -1

    pml_index = 1
    for z in 1:nz, x in 1:nx 
        
        if z in pml_zid ||
           x in pml_xid 

           push!(pml_points, (z,x))
           pml_lookup[z,x] = pml_index
           pml_index += 1

        end
    end
    
    pml_points = union(pml_points)
    @assert length(pml_points) == pml_index - 1 # because the index started with 1
 
    # types 
    T1 = typeof(x0)
    T2 = typeof(dx)
    T3 = typeof(xcoords)
    T4 = typeof(pml_lookup)

    domain = Domain{T1,T2,T3,T4}(
                    x_inner_start,dx,x_inner_end,nx,xcoords,
                    z_inner_start,dz,z_inner_end,nz,zcoords,
                    dim, inner_id, 
                    pml_id, pml_lookup)

    return domain 

end;