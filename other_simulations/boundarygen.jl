using BayesianShockEstimation
using CSV
using DataFrames

# Set Boundary Settings 
σ_w = 0.2 
σ_v = 1 
spatial_steps = 400 
boundary = BoundaryParams(σ_w, σ_v)

# Set PDE Settings
stepsize_x = 0.05 
stepsize_t = 0.005 
stepcount_f = 2000 
ballsize = 0.051 
threshold = 1

# Generate the boundary and observations
true_boundary = boundary.σ_w * randn(spatial_steps) 
for i in 2:length(true_boundary) 
    true_boundary[i] += true_boundary[i-1] 
end 
true_boundary .-= sum(true_boundary)/spatial_steps 
observations = true_boundary + boundary.σ_v * randn(spatial_steps)

# Save the boundary and observations
df = DataFrame(true_boundary = true_boundary, observations = observations)
CSV.write("data/boundary.csv", df)

df = DataFrame( σ_w = σ_w, 
                σ_v = σ_v, 
                spatial_steps = spatial_steps,
                stepsize_x = stepsize_x,
                stepsize_t = stepsize_t,
                stepcount_f = stepcount_f,
                ballsize = ballsize,
                threshold = threshold)

CSV.write("data/params.csv", df)

