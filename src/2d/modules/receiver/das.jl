mutable struct Fibers
    n::Int64
    orientation::Any
    data::Any
    coords::Any
    ids::Any
    cinv::Any
end

mutable struct DAS
    x_aligned::Fibers 
    y_aligned::Fibers 
end

function register_fibers(settings::Settings, domain::Domain, elastic::Elastic, time::Time, fiber_orientation)
 
    fiber_coords = []
    fiber_ids = []
    fiber_cinv = []
    fiber_data = []
    nfibers = 0

    if !isnothing(settings.config["receivers"]["das"][fiber_orientation])

        nfibers = length(settings.config["receivers"]["das"][fiber_orientation])

        for i in eachindex(settings.config["receivers"]["das"][fiber_orientation])

            fiber = settings.config["receivers"]["das"][fiber_orientation][i]

            if fiber_orientation == "x_aligned"
                start,step,ende = split(fiber["x_range"], ":")    
                fiber_x = parse(settings.float, start):parse(settings.float,step):parse(settings.float,ende)
                fiber_x = collect(fiber_x)
                fiber_y = fiber["y"]
                push!(fiber_coords, [fiber_x, fiber_y])

            elseif fiber_orientation == "y_aligned"
                fiber_x = fiber["x"]
                start,step,ende = split(fiber["y_range"], ":")    
                fiber_y = parse(settings.float, start):parse(settings.float,step):parse(settings.float,ende)
                fiber_y = collect(fiber_y)
                push!(fiber_coords, [fiber_x, fiber_y])
            end
        end

        for i in eachindex(fiber_coords)
            fiber_x, fiber_y = fiber_coords[i]
            if fiber_orientation == "x_aligned"
                fiber_x_ids = [argmin(abs.(domain.xcoords .- v)) for v in fiber_x]
                fiber_y_ids =  argmin(abs.(domain.ycoords .- fiber_y))
                push!(fiber_ids, [fiber_x_ids, fiber_y_ids])
            elseif fiber_orientation == "y_aligned"
                fiber_x_ids =  argmin(abs.(domain.xcoords .- fiber_x))
                fiber_y_ids = [argmin(abs.(domain.ycoords .- v)) for v in fiber_y]
                push!(fiber_ids, [fiber_x_ids, fiber_y_ids])
            end

        end

        for k in eachindex(fiber_ids)
            fiber_x_ids, fiber_y_ids = fiber_ids[k]
            fiber_dim = length(fiber_ids[k][1]) * length(fiber_ids[k][2])
            Cinv = zeros(settings.float, fiber_dim , 2, 2)
            c = 1
            for yid in fiber_y_ids, xid in fiber_x_ids 
                C = [elastic.c11[yid, xid] elastic.c13[yid, xid];
                     elastic.c13[yid, xid] elastic.c33[yid, xid]]
                Cinv[c,:,:] .= inv(C)
                c += 1
            end
            push!(fiber_cinv, Cinv)
        end

        for i in eachindex(fiber_ids)
            fiber_i_dim = length(fiber_ids[i][1]) * length(fiber_ids[i][2])
            data_i = zeros(settings.float, fiber_i_dim, time.nt)
            push!(fiber_data, data_i)
        end

    end
    
    return nfibers, fiber_coords, fiber_ids, fiber_cinv, fiber_data
end;


function init_das(settings::Settings, domain::Domain, elastic::Elastic, time::Time)

    nfibers, fiber_coords, fiber_ids, 
    fiber_cinv, fiber_data = register_fibers(settings, domain, elastic, time, "x_aligned")
    xfibers = Fibers(nfibers, "x_aligned", fiber_data, fiber_coords, fiber_ids, fiber_cinv)

    nfibers, fiber_coords, fiber_ids, 
    fiber_cinv, fiber_data = register_fibers(settings, domain, elastic, time, "y_aligned")
    yfibers = Fibers(nfibers, "y_aligned", fiber_data, fiber_coords, fiber_ids, fiber_cinv)

    das = DAS(xfibers, yfibers)

    return das 
end;


function save_das!(das::DAS, sxx, syy, ti, float = Float64)

    # x strains 
    for i in 1:das.x_aligned.n

        sxx_along_fiber = Array{float}(sxx[das.x_aligned.ids[i][2], das.x_aligned.ids[i][1]])
        syy_along_fiber = Array{float}(syy[das.x_aligned.ids[i][2], das.x_aligned.ids[i][1]])

        for channel in eachindex(das.x_aligned.ids[i][1])

            sigma = [sxx_along_fiber[channel] syy_along_fiber[channel]]
            das.x_aligned.data[i][channel, ti] = (sigma * das.x_aligned.cinv[i][channel,1,:])[1]
        
        end 
    end 

    # y strains 
    for i in 1:das.y_aligned.n

        sxx_along_fiber = Array{float}(sxx[das.y_aligned.ids[i][2], das.y_aligned.ids[i][1]])
        syy_along_fiber = Array{float}(syy[das.y_aligned.ids[i][2], das.y_aligned.ids[i][1]])

        for channel in eachindex(das.y_aligned.ids[i][2])

            sigma = [sxx_along_fiber[channel] syy_along_fiber[channel]]
            das.y_aligned.data[i][channel, ti] = (sigma * das.y_aligned.cinv[i][channel,2,:])[1]
        
        end 
    end 
end