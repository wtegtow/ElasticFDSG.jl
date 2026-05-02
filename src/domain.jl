struct Domain{N}
    shape::NTuple{N, Int}
    coordinates::NTuple{N, AbstractVector}
    inner_ids::NTuple{N, Vector{Int}}
    pml_lookup::AbstractArray{Int, N}
end 

# dim-specific helpers
_coord_vecs(vm::VelocityModel2D, fp) = (fp.(vm.X[:, 1]),   fp.(vm.Z[1, :]))
_coord_vecs(vm::VelocityModel3D, fp) = (fp.(vm.X[:, 1, 1]), fp.(vm.Y[1, :, 1]), fp.(vm.Z[1, 1, :]))
_axis_names(::VelocityModel2D) = ("x", "z")
_axis_names(::VelocityModel3D) = ("x", "y", "z")

_expand_coords(inner, pad_s, pad_e, d) =
    LinRange(inner[begin] - pad_s * d, inner[end] + pad_e * d, length(inner) + pad_s + pad_e)

_inner_ids(inner, expanded) = [argmin(abs.(v .- expanded)) for v in inner]

function _pml_ranges(abs_s, abs_e, N, nlayer, id_start, id_end)
    r = UnitRange{Int}[]
    abs_s && push!(r, (N+1):(N+nlayer))
    abs_e && push!(r, (id_end+1):(id_end+nlayer))
    return r
end

function init_domain(config::Config, velmod::VelocityModel)::Domain
    
    N          = config.dim
    boundaries = config.dict["boundaries"]
    nlayer_pml = Int(boundaries["pml_layer"])
    order_N    = Int(config.dict["settings"]["spatial_derivative_order"])
    fp         = eval(Symbol(config.dict["settings"]["precision"]))

    inner_vecs = _coord_vecs(velmod.fields, fp)
    axis_names = _axis_names(velmod.fields)

    spacings = map(v -> abs(v[2] - v[1]), inner_vecs)
    abs_s    = map(k -> get(boundaries, k * "start", "none") == "absorbing", axis_names)
    abs_e    = map(k -> get(boundaries, k * "end",   "none") == "absorbing", axis_names)
    pad_s    = map(a -> a ? nlayer_pml + order_N : order_N, abs_s)
    pad_e    = map(a -> a ? nlayer_pml + order_N : order_N, abs_e)

    coords    = map(_expand_coords, inner_vecs, pad_s, pad_e, spacings)
    shape     = map(length, coords)
    inner_ids = map(_inner_ids, inner_vecs, coords)

    pml_ranges = map((as, ae, iid) -> _pml_ranges(as, ae, order_N, nlayer_pml, iid[begin], iid[end]),
                     abs_s, abs_e, inner_ids)

    pml_sets   = map(r -> isempty(r) ? Set{Int}() : Set(reduce(vcat, collect.(r))), pml_ranges)
    pml_lookup = fill(-1, Tuple(shape))
    pml_index  = 1
    for idx in CartesianIndices(Tuple(shape))
        if any(idx[k] in pml_sets[k] for k in 1:N)
            pml_lookup[idx] = pml_index
            pml_index += 1
        end
    end

    return Domain{N}(Tuple(shape), coords, inner_ids, pml_lookup)
end;