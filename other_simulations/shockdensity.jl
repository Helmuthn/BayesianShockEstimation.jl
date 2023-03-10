using BayesianShockEstimation
using Serialization
using CSV
using DataFrames

# Read in parameters
parameters = DataFrame(CSV.File("data/params.csv"))
dx = parameters[!,"stepsize_x"][1]
dt = parameters[!,"stepsize_t"][1]
stepcount = parameters[!,"stepcount_f"][1]
threshold = parameters[!,"threshold"][1]

# Set the radius range
radii = (1:.2:3) .* max(dx, dt) .+ min(dx,dt)/2


max_solutions = 1000
thread_count = Threads.nthreads()
Threads.@threads for id in 1:thread_count
    # Read in data
    directory = "data/thread$(id)"

    data = deserialize(directory * "/1.jls")
    shockcount = zeros(size(data))

    for (j, radius) in enumerate(radii)
        offsets = getoffsets(dx, dt, radius)
        shockcount .= 0
        
        for i in 1:max_solutions
            filename = directory * "/$(i).jls"
            data = deserialize(filename)

            shockcount .+= rangefilter(data, offsets) .> threshold
        end
        shockcount ./= (max_solutions * 2 * radius)
        filename = directory * "/result_r$(j).jls"
        serialize(filename, shockcount)
    end
end

