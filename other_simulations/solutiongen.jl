using BayesianShockEstimation
using CSV
using DataFrames
using Serialization

# Read in parameters
parameters = DataFrame(CSV.File("data/params.csv"))
σ_w = parameters[!,"σ_w"][1]
σ_v = parameters[!,"σ_v"][1]
dx = parameters[!,"stepsize_x"][1]
dt = parameters[!,"stepsize_t"][1]
stepcount = parameters[!,"stepcount_f"][1]


max_solutions = 1000
thread_count = Threads.nthreads()
Threads.@threads for id in 1:thread_count
    # Read in data
    data = DataFrame(CSV.File("data/boundary_samples$(id).csv"))

    # Ensure data is properly parsed
    typelist = eltype.(eachcol(data))
    for (i, typecheck) in enumerate(typelist)
        if typecheck <: AbstractString 
            data[!,i] = tryparse.(Float64, data[!,i])
        end
    end

    # Run simulations
    directory = "data/thread$(id)"
    rm(directory, recursive=true, force=true)
    mkdir(directory)
    for (i, initialcondition) in enumerate(eachcol(data))
        if i <= max_solutions
         solution = godunov_burgers_1D(initialcondition, dx, dt, stepcount)

         filename = directory * "/$(i).jls"
         serialize(filename, solution)
        end
    end
end

