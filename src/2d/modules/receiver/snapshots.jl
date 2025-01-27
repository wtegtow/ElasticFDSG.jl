mutable struct Snapshots{T}
    XY::T 
    fieldnames::Vector{String}
    t::Vector{Float64}
    tid::Vector{Int64}
    tid_map::Dict{Int64, Int64}
end 

function init_snapshots(settings::Settings, domain::Domain, time::Time)

    tsnap = settings.config["receivers"]["snapshots"]["times"]
    fields = settings.config["receivers"]["snapshots"]["fields"]

    if isnothing(tsnap) # some dummy snapshots
        tsnap = [time.t[1]]
        fields = ["vx"]
    end 

    tsnap_ids = [argmin(abs.(time.t .- tsnap_i)) for tsnap_i in tsnap]
    tid_map = Dict(tsnap_ids[i] => i for i in eachindex(tsnap_ids))

    nfields = length(fields)
    nsnaps = length(tsnap)
    XY = zeros(settings.float, nsnaps, nfields, domain.ny, domain.nx);

    T = typeof(XY)
    snapshots = Snapshots{T}(XY, fields, tsnap, tsnap_ids, tid_map);

    return snapshots
end;



function save_snapshots!(snapshots::Snapshots, vx, vy, sxx, sxy, syy, ti, float = Float64)
    if haskey(snapshots.tid_map, ti)
        tid = snapshots.tid_map[ti]
        fieldnum = 1
        for field in snapshots.fieldnames
            if field == "vx"
                snapshots.XY[tid,fieldnum,:,:] .= Array{float}(vx)
            elseif field == "vy"
                snapshots.XY[tid,fieldnum,:,:] .= Array{float}(vy)
            elseif field == "sxx"
                snapshots.XY[tid,fieldnum,:,:] .= Array{float}(sxx)
            elseif field == "sxy"
                snapshots.XY[tid,fieldnum,:,:] .= Array{float}(sxy)
            elseif field == "syy"
                snapshots.XY[tid,fieldnum,:,:] .= Array{float}(syy)
            end
            fieldnum += 1
        end
    end
end;