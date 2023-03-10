using BayesianShockEstimation
using CSV
using DataFrames

# Read in parameters
parameters = DataFrame(CSV.File("data/params.csv"))
σ_w = parameters[!,"σ_w"][1]
σ_v = parameters[!,"σ_v"][1]

# Read in data
data = DataFrame(CSV.File("data/boundary.csv"))
observations = data[!,"observations"]
true_boundary = data[!,"true_boundary"]

# Get directory name
directory = "data/node1"
#directory = "data/node$(ARGS[1])"



# Generate the initial boundary conditions
rts_sample_count = 10000
thread_count = Threads.nthreads()

estimates, kalman_variances, rts_variances = RTS_Smooth(observations, σ_w, σ_v) 
boundary_samples = zeros(thread_count, length(true_boundary),rts_sample_count) 
Threads.@threads for id in 1:thread_count
    for i in 1:rts_sample_count 
        boundary_samples[id, :, i] .= RTS_sample(estimates, kalman_variances, rts_variances, σ_w) 
    end
end

for id in 1:thread_count
    boundary_df = DataFrame(boundary_samples[id, :, :], :auto)
    CSV.write("data/boundary_samples$(id).csv", boundary_df)
end

