include(joinpath(@__DIR__,"forces.jl")) 
include(joinpath(@__DIR__,"cpu/solver.jl")) 
include(joinpath(@__DIR__,"metal/solver.jl"))
include(joinpath(@__DIR__,"cuda/solver.jl"))

function init_solver(fdsg2d::FDSG2D)
    
    device = fdsg2d.settings.device 
    
    # cpu
    if device == "cpu"
        solver = solver_cpu! 
        return solver 

    # metal 
    elseif device == "gpu-metal" 
        solver = solver_metal
        return solver
    
    # cuda 
    elseif device == "gpu-cuda" 
        solver = solver_cuda
        return solver
    end
end

