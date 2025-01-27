mutable struct Geophones{T, U}
    data::T
    n::Int64
    coords::U
    xids::Vector{Int64} 
    yids::Vector{Int64} 
    zids::Vector{Int64}
end


function init_geophones(settings::Settings, 
                        domain::Domain,
                        time::Time)

    geo_coords = nothing
    if !isnothing(settings.config["receivers"]["geophones"])
    
        ngeo = length(settings.config["receivers"]["geophones"])
        geo_coords = zeros(settings.float, 3, ngeo)
    
        for i in 1:ngeo
            geo_coords[1, i] = settings.config["receivers"]["geophones"][i]["x"]
            geo_coords[2, i] = settings.config["receivers"]["geophones"][i]["y"]
            geo_coords[3, i] = settings.config["receivers"]["geophones"][i]["z"]
        end 

    else
        geo_coords = settings.float.([domain.xcoords[1]; domain.ycoords[1]; domain.zcoords[1]]) # some dummy coordinates
    end;

    geo_x = geo_coords[1,:]
    geo_y = geo_coords[2,:]
    geo_z = geo_coords[3,:]

    x0 = domain.xcoords[begin]
    xend = domain.xcoords[end]
    y0 = domain.ycoords[begin]
    yend = domain.ycoords[end]
    z0 = domain.zcoords[begin]
    zend = domain.zcoords[end]

    # are receivers within the domain?
    if all(geo_x .>= x0)   && 
       all(geo_x .<= xend) && 
       all(geo_y .>= y0)   && 
       all(geo_y .<= yend) && 
       all(geo_z .>= z0)   &&
       all(geo_z .<= zend)
        ;
    else
        error("Some geophone receivers are defined outsite the domain.")
    end

    geo_x_ids = [argmin(abs.(domain.xcoords .- geo_xi))  for geo_xi in geo_x]
    geo_y_ids = [argmin(abs.(domain.ycoords .- geo_yi))  for geo_yi in geo_y]
    geo_z_ids = [argmin(abs.(domain.zcoords .- geo_zi))  for geo_zi in geo_z]

    ngeo = length(geo_x_ids)
    data = zeros(settings.float, ngeo, 3, time.nt)

    T = typeof(data)
    U = typeof(geo_coords)
    geophones = Geophones{T, U}(data, ngeo, geo_coords, geo_x_ids, geo_y_ids, geo_z_ids)
    
    return geophones

end;

function save_geophones!(geophones::Geophones, vx, vy, vz, ti)
    for n in 1:geophones.n  
        geophones.data[n,1,ti] = vx[geophones.zids[n],geophones.yids[n],geophones.xids[n]]
        geophones.data[n,2,ti] = vy[geophones.zids[n],geophones.yids[n],geophones.xids[n]]
        geophones.data[n,3,ti] = vz[geophones.zids[n],geophones.yids[n],geophones.xids[n]]
    end
end;
