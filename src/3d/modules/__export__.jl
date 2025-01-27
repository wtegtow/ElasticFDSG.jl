# io
include(joinpath(@__DIR__, "io/read.jl"))
include(joinpath(@__DIR__, "io/write.jl"))
# settings
include(joinpath(@__DIR__, "settings.jl"))
include(joinpath(@__DIR__, "device.jl"))
# setup
include(joinpath(@__DIR__, "domain.jl"))
include(joinpath(@__DIR__, "elastic.jl"))
include(joinpath(@__DIR__, "time.jl"))
include(joinpath(@__DIR__, "source.jl"))
include(joinpath(@__DIR__, "fields.jl"))
include(joinpath(@__DIR__, "pml.jl"))
# receiver
include(joinpath(@__DIR__, "receiver/das.jl"))
include(joinpath(@__DIR__, "receiver/geophones.jl"))
include(joinpath(@__DIR__, "receiver/snapshots.jl"))
# solver 
include(joinpath(@__DIR__, "solver/__export__.jl"))
# iwindow 
include(joinpath(@__DIR__, "iwindow.jl"))