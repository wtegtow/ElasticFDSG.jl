mutable struct Snapshots{T}

    XY::T 
    XZ::T
    YZ::T

    fieldnames::Vector{String}

    t::Vector{Float64}
    tid::Vector{Int64}
    tid_map::Dict{Int64, Int64}

    xidx::Int64
    yidx::Int64
    zidx::Int64
end 

function init_snapshots(settings::Settings, domain::Domain, time::Time)

    tsnap = settings.config["receivers"]["snapshots"]["times"]
    fields = settings.config["receivers"]["snapshots"]["fields"]
    origin = settings.config["receivers"]["snapshots"]["origins"]
    
    if isnothing(tsnap) # some dummy snapshots
        tsnap = [time.t[1]]
        fields = ["vx"]
        origin = [abs(domain.xcoords[begin] - domain.xcoords[end])/2 
                  abs(domain.ycoords[begin] - domain.ycoords[end])/2 
                  abs(domain.zcoords[begin] - domain.zcoords[end])/2]
    end 


    # time indices of snapshots
    tsnap_ids = [argmin(abs.(time.t .- tsnap_i)) for tsnap_i in tsnap]
    tid_map = Dict(tsnap_ids[i] => i for i in eachindex(tsnap_ids))

    # spatial indices for snapshot planes 
    xidx = argmin(abs.(domain.xcoords .- origin["x"])) 
    yidx = argmin(abs.(domain.ycoords .- origin["y"]))
    zidx = argmin(abs.(domain.zcoords .- origin["z"]))

    nfields = length(fields)
    nsnaps = length(tsnap)

    XY = zeros(settings.float, nsnaps, nfields, domain.ny, domain.nx);
    XZ = zeros(settings.float, nsnaps, nfields, domain.nz, domain.nx);
    YZ = zeros(settings.float, nsnaps, nfields, domain.nz, domain.ny);

    T = typeof(XY)
    snapshots = Snapshots{T}(XY,XZ,YZ,
                            fields, 
                            tsnap, tsnap_ids, tid_map,
                            xidx, yidx, zidx);

    return snapshots
end;


function save_snapshots!(snapshots::Snapshots, 
                         vx, vy, vz, sxx, sxy, sxz, syy, syz, szz, 
                         ti, FLOAT = Float64)

    if haskey(snapshots.tid_map, ti)

        tid = snapshots.tid_map[ti]
        fieldnum = 1

        for field in snapshots.fieldnames

            if field == "vx"
                field_ = vx
            elseif field == "vy"
                field_ = vy
            elseif field == "vz"
                field_ = vz 

            elseif field == "sxx"
                field_ = sxx
            elseif field == "sxy"
                field_ = sxy
            elseif field == "sxz"
                field_ = sxz
            elseif field == "syy"
                field_ = syy
            elseif field == "syz"
                field_ = syz
            elseif field == "szz"
                field_ = szz
            end

            snapshots.XY[tid, fieldnum, :, :] .= Array{FLOAT}(field_[snapshots.zidx,:,:])
            snapshots.XZ[tid, fieldnum, :, :] .= Array{FLOAT}(field_[:,snapshots.yidx,:])
            snapshots.YZ[tid, fieldnum, :, :] .= Array{FLOAT}(field_[:,:,snapshots.xidx])

            fieldnum += 1
        end
    end
end;