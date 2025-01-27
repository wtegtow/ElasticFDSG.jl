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
    z_aligned::Fibers
end


function register_fibers(settings::Settings, domain::Domain, elastic::Elastic, time::Time, fiber_orientation)
 
    fiber_coords = []
    fiber_ids = []
    fiber_cinv = []
    fiber_data = []
    nfibers = 0

    if !isnothing(settings.config["receivers"]["das"][fiber_orientation])

        nfibers = length(settings.config["receivers"]["das"][fiber_orientation])
       
        # coordinates of fibers
        for i in 1:nfibers

            fiber = settings.config["receivers"]["das"][fiber_orientation][i]

            if fiber_orientation == "x_aligned"
                start,step,ende = split(fiber["x_range"], ":")    
                fiber_x = parse(settings.float, start):parse(settings.float,step):parse(settings.float,ende)
                fiber_x = collect(fiber_x)
                fiber_y = fiber["y"]
                fiber_z = fiber["z"]
                push!(fiber_coords, [fiber_x, fiber_y, fiber_z])

            elseif fiber_orientation == "y_aligned"
                fiber_x = fiber["x"]
                start,step,ende = split(fiber["y_range"], ":")    
                fiber_y = parse(settings.float, start):parse(settings.float,step):parse(settings.float,ende)
                fiber_y = collect(fiber_y)
                fiber_z = fiber["z"]
                push!(fiber_coords, [fiber_x, fiber_y, fiber_z])
            
            elseif fiber_orientation == "z_aligned"
                fiber_x = fiber["x"]
                fiber_y = fiber["y"]
                start,step,ende = split(fiber["z_range"], ":")    
                fiber_z = parse(settings.float, start):parse(settings.float,step):parse(settings.float,ende)
                fiber_z = collect(fiber_z)
                push!(fiber_coords, [fiber_x, fiber_y, fiber_z])
            end
        end

        # indices of fibers
        for i in 1:nfibers
            fiber_x, fiber_y, fiber_z = fiber_coords[i]
            if fiber_orientation == "x_aligned"
                fiber_x_ids = [argmin(abs.(domain.xcoords .- v)) for v in fiber_x]
                fiber_y_ids =  argmin(abs.(domain.ycoords .- fiber_y))
                fiber_z_ids =  argmin(abs.(domain.zcoords .- fiber_z))
                push!(fiber_ids, [fiber_x_ids, fiber_y_ids, fiber_z_ids])

            elseif fiber_orientation == "y_aligned"
                fiber_x_ids =  argmin(abs.(domain.xcoords .- fiber_x))
                fiber_y_ids = [argmin(abs.(domain.ycoords .- v)) for v in fiber_y]
                fiber_z_ids =  argmin(abs.(domain.zcoords .- fiber_z))
                push!(fiber_ids, [fiber_x_ids, fiber_y_ids, fiber_z_ids])
            
            elseif fiber_orientation == "z_aligned"
                fiber_x_ids =  argmin(abs.(domain.xcoords .- fiber_x))
                fiber_y_ids =  argmin(abs.(domain.ycoords .- fiber_y))
                fiber_z_ids = [argmin(abs.(domain.zcoords .- v)) for v in fiber_z]
                push!(fiber_ids, [fiber_x_ids, fiber_y_ids, fiber_z_ids])
            end
        end

        # inverted C-Tensors along fiber coordinates
        for k in 1:nfibers
            fiber_x_ids, fiber_y_ids, fiber_z_ids = fiber_ids[k]
            fiber_dim = length(fiber_x_ids) * length(fiber_y_ids) * length(fiber_z_ids)
            Cinv = zeros(settings.float, fiber_dim , 3, 3)
            c = 1
            for zid in fiber_z_ids, yid in fiber_y_ids, xid in fiber_x_ids 

                COORD = (zid, yid, xid)
                C = [elastic.c11[COORD...] elastic.c12[COORD...] elastic.c13[COORD...];
                     elastic.c12[COORD...] elastic.c22[COORD...] elastic.c23[COORD...];
                     elastic.c13[COORD...] elastic.c23[COORD...] elastic.c33[COORD...]]

                Cinv[c,:,:] .= inv(C)
                c += 1
            end
            push!(fiber_cinv, Cinv)
        end

        # data arrays 
        for i in 1:nfibers
            fiber_i_dim = length(fiber_ids[i][1]) * length(fiber_ids[i][2]) * length(fiber_ids[i][3])
            data_i = zeros(settings.float, fiber_i_dim, time.nt)
            push!(fiber_data, data_i)
        end

    end
    
    return nfibers, fiber_coords, fiber_ids, fiber_cinv, fiber_data
end;


function init_das(settings::Settings, domain::Domain, elastic::Elastic, time::Time)

    # x aligned fibers
    nfibers, fiber_coords, fiber_ids, 
    fiber_cinv, fiber_data = register_fibers(settings, domain, elastic, time, "x_aligned")

    xfibers = Fibers(nfibers, "x_aligned", fiber_data, fiber_coords, fiber_ids, fiber_cinv)

    # y aligned fibers
    nfibers, fiber_coords, fiber_ids, 
    fiber_cinv, fiber_data = register_fibers(settings, domain, elastic, time, "y_aligned")

    yfibers = Fibers(nfibers, "y_aligned", fiber_data, fiber_coords, fiber_ids, fiber_cinv)

    # z aligned fibers
    nfibers, fiber_coords, fiber_ids, 
    fiber_cinv, fiber_data = register_fibers(settings, domain, elastic, time, "z_aligned")

    zfibers = Fibers(nfibers, "z_aligned", fiber_data, fiber_coords, fiber_ids, fiber_cinv)

    # gather all fibers in DAS struct
    das = DAS(xfibers, yfibers, zfibers)

    return das 
end;


function save_das!(das::DAS, sxx, syy, szz, ti, FLOAT = Float64)

    # x strains 
    for i in 1:das.x_aligned.n

        COORD = (das.x_aligned.ids[i][3], das.x_aligned.ids[i][2], das.x_aligned.ids[i][1])
        sxx_along_fiber = Array{FLOAT}(sxx[COORD...])
        syy_along_fiber = Array{FLOAT}(syy[COORD...])
        szz_along_fiber = Array{FLOAT}(szz[COORD...])

        for channel in eachindex(das.x_aligned.ids[i][1])

            sigma = [sxx_along_fiber[channel] syy_along_fiber[channel] szz_along_fiber[channel]]
            das.x_aligned.data[i][channel, ti] = (sigma * das.x_aligned.cinv[i][channel,1,:])[1]
        
        end 
    end 

    # y strains 
    for i in 1:das.y_aligned.n

        COORD = (das.y_aligned.ids[i][3], das.y_aligned.ids[i][2], das.y_aligned.ids[i][1])
        sxx_along_fiber = Array{FLOAT}(sxx[COORD...])
        syy_along_fiber = Array{FLOAT}(syy[COORD...])
        szz_along_fiber = Array{FLOAT}(szz[COORD...])

        for channel in eachindex(das.y_aligned.ids[i][2])

            sigma = [sxx_along_fiber[channel] syy_along_fiber[channel] szz_along_fiber[channel]]
            das.y_aligned.data[i][channel, ti] = (sigma * das.y_aligned.cinv[i][channel,2,:])[1]
        
        end 
    end 

     # z strains 
     for i in 1:das.z_aligned.n

        COORD = (das.z_aligned.ids[i][3], das.z_aligned.ids[i][2], das.z_aligned.ids[i][1])
        sxx_along_fiber = Array{FLOAT}(sxx[COORD...])
        syy_along_fiber = Array{FLOAT}(syy[COORD...])
        szz_along_fiber = Array{FLOAT}(szz[COORD...])

        for channel in eachindex(das.z_aligned.ids[i][3])

            sigma = [sxx_along_fiber[channel] syy_along_fiber[channel] szz_along_fiber[channel]]
            das.z_aligned.data[i][channel, ti] = (sigma * das.z_aligned.cinv[i][channel,3,:])[1]
        
        end 
    end 
end