include(joinpath(@__DIR__, "io.jl"))
include(joinpath(@__DIR__, "settings.jl"))
include(joinpath(@__DIR__, "domain.jl"))
include(joinpath(@__DIR__, "elastic.jl"))
include(joinpath(@__DIR__, "time.jl"))
include(joinpath(@__DIR__, "source.jl"))
include(joinpath(@__DIR__, "fields.jl"))
include(joinpath(@__DIR__, "pml.jl"))
include(joinpath(@__DIR__, "receiver.jl"))
include(joinpath(@__DIR__,"forces.jl"));
include(joinpath(@__DIR__,"forward.jl"));
include(joinpath(@__DIR__,"solver.jl"));
include(joinpath(@__DIR__, "iwindow.jl"))