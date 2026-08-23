include(joinpath(@__DIR__, "forward.jl"))
include(joinpath(@__DIR__, "forces.jl"))

# Cast all kernel-relevant arrays to device (GPU)
function to_device!(fdsg::FDSG)
    AT = fdsg.device.array     # the array constructor for the target device
    AT == Array && return fdsg # already on CPU, no need to move
    N  = fdsg.config.dim

    d = fdsg.domain
    fdsg.domain  = Domain{N}(d.shape, d.coordinates, d.inner_ids, AT(d.pml_lookup))

    e = fdsg.elastic
    fdsg.elastic = Elastic{N}(e.c_tensors, AT(e.c_lookup), e.vmax, e.vmin)

    # parametric mutable structs: reconstruct so the type parameter updates
    FT = typeof(fdsg.fields).name.wrapper
    fdsg.fields = FT((AT(getfield(fdsg.fields, fn)) for fn in fieldnames(typeof(fdsg.fields)))...)

    PT = typeof(fdsg.pml).name.wrapper
    fdsg.pml = PT((AT(getfield(fdsg.pml, fn)) for fn in fieldnames(typeof(fdsg.pml)))...)

    return fdsg
end

function to_host!(fdsg::FDSG)
    fdsg.device.array == Array && return fdsg
    N = fdsg.config.dim

    d = fdsg.domain
    fdsg.domain  = Domain{N}(d.shape, d.coordinates, d.inner_ids, Array(d.pml_lookup))

    e = fdsg.elastic
    fdsg.elastic = Elastic{N}(e.c_tensors, Array(e.c_lookup), e.vmax, e.vmin)

    FT = typeof(fdsg.fields).name.wrapper
    fdsg.fields = FT((Array(getfield(fdsg.fields, fn)) for fn in fieldnames(typeof(fdsg.fields)))...)

    PT = typeof(fdsg.pml).name.wrapper
    fdsg.pml = PT((Array(getfield(fdsg.pml, fn)) for fn in fieldnames(typeof(fdsg.pml)))...)

    return fdsg
end

function solve!(fdsg::FDSG, block_size=nothing)

    fdsg = to_device!(fdsg)
    params = init_simparams(fdsg)
    backend = fdsg.device.backend
    
    if isnothing(block_size)
        block_size = fdsg.config.dim == 2 ? (32, 32) : (8, 8, 8)
    end

    fields    = fdsg.fields
    pml       = fdsg.pml
    elastic   = fdsg.elastic
    domain    = fdsg.domain
    time      = fdsg.time
    source    = fdsg.source
    geophones = fdsg.geophones
    das       = fdsg.das
    snapshots = fdsg.snapshots

    showinfo  = get(fdsg.config.dict["settings"], "verbose", true)
    if showinfo
        prog = Progress(time.nt; showspeed=true, desc="Solving...")
    end

    for ti in 1:time.nt

        update_velocities!(fields, pml, elastic, domain, time, params, backend, block_size)
        KernelAbstractions.synchronize(backend)

        update_stresses!(fields, pml, elastic, domain, time, params, backend, block_size)
        KernelAbstractions.synchronize(backend)

        GPUArrays.@allowscalar _stress_glut!(fields, source, domain, time, ti)

        GPUArrays.@allowscalar save_geophones!(geophones, fields, ti)
        save_das!(das, fields, domain, params.N_fd, ti)
        save_snapshots!(snapshots, fields, ti)

        showinfo && next!(prog)
    end

    fdsg = to_host!(fdsg)

    return fdsg
end
