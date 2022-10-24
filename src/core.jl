using SpecialFunctions: gamma
using Interpolations: linear_interpolation
# This file contains the core contributions of the work.

"""
    ShockParams

Struct containing all need information for deterministic shock simulations

# Parameters
 - `dynamics`: Representation of the dynamics
 - `stepsize`: Grid size
 - `ballsize`: Integration ball size
 - `threshold`: Minimum shock threshold
 - `dimension`: System dimensionality

"""
struct ShockParams
    dynamics
    stepsize
    ballsize
    threshold
    dimension::Int
end
export ShockParams

"""
    BoundaryParams

Struct containing all needed information for boundary sampling

# Parameters
 - `σ_w` : System noise standard deviation
 - `σ_v` : Measurement noise standard deviation
 - `t` : Sample interval

"""
struct BoundaryParams
    σ_w
    σ_v
    t
end
export BoundaryParams

"""
    ShockDensity(systemparams, samples)

Generates an approximation of the shock density function.

# Arguments
 - `systemparams` : Struct of parameters for system simulation
 - `samples`: Number of monte carlo samples

# Returns
A list of approximate densities
"""
function ShockDensity(systemparams::ShockParams, boundaryparams::BoundaryParams, samples)
    shockcounts = zeros(length(systemparams.positions))

    # Compute RTS Smoothing Details
    σ_w = boundaryparams.σ_w
    σ_v = boundaryparams.σ_v
    estimates, variances, RTSvariances = RTS_Smooth(samples, σ_w, σ_v)
    
    # Monte Carlo Integration Based On Stochastic Boundaries
    for i in range(samples)
        boundary = RTS_sample(estimates, variances, RTSvariances, σ_w)
        IncrementShockCounts!(shockcounts, systemparams, boundary)
    end

    # Normalize counts to the size of the ϵ-ball
    volume = BallVolume(systemparams)
    normalized_counts = shockcounts ./ volume

    # Convert to a density by normalizing to the number of samples
    shockdensity = normalized_counts ./ samples
    
    return shockdensity 
end
export ShockDensity


"""
    BallVolume(systemparams::ShockParams)::Float64

Helper function that computes the volume of the integration ball
"""
function BallVolume(systemparams::ShockParams)::Float64
    dim = systemparams.dimension
    R = systemparams.ballsize
    return π^(dim/2) * R^dim / gamma(dim/2 + 1) 
end


"""
    IncrementShockCounts!(shockcounts, systemparams, boundary)

Computes shocks for a given boundary.

# Arguments
 - `shockcounts`: Array of current shock counts
 - `systemparams` : Struct of parameters for system simulation
 - `T`: Sampling Interval
 - `boundary_values`: Boundary samples

# Returns
A list of approximate densities

# Notes
Uses the method of characteristics to find shocks.
"""
function IncrementShockCounts!(shockcounts, systemparams::ShockParams, T, boundary_values)

    stepsize = 0.1
    # Initialize a collection of ODEs evolving over space
    t_lst = Array(T * (1:0.1:length(boundary_values)))
    x_lst = zeros(length(t_lst))
    boundary = linear_interpolation(1:length(boundary_values), boundary_values)
    u_lst = boundary.(1:0.1:length(boundary_values))

    for space_step in 1:size(shockcounts)[2]    
        # Update all characteristic trajectories
        x_lst .+= stepsize
        t_lst .+= stepsize ./ u_lst
        
        removed = []
        introduced = []

        # Search for any out of order trajectories
        for i in 1:length(t_lst)-1
            # If they've permuted, then it's a shockwave
            if t_lst[i+1] < t_lst[i]
                push!(removed, (i, i+1)) 
            end
            # If the gap is larger than twice the stepsize, 
            # Add another characteristic curve
            if t_lst[i+1] - t_lst[i] > 2 * stepsize
                
            end
        end

        for shock in removed
            # Increment nearby shockcounts
            
        end

        for trajectory in introduced
            # Introduce new ODE trajectories
            
        end
    end
    

#    # Find all potential interesection positions
#    # Store as in an array of tuples (loc, time, t_1, t_2)
#    intersections = []
#    ind = 0
#    for i in 1:length(boundary_values)
#        for j in i+1:length(boundary_values)
#            intersect = ComputeIntersection(T * i, boundary_values[i], T * j, boundary_values[j]) 
#            if !isnothing(intersect)
#                ind += 1
#                intersections[ind] = (intersect[2], intersect[1], i, j)
#            end
#        end
#    end
#
#    # Sort the array by position of shock
#    sort!(intersections, by = x -> x[1])
#
#    # Prune the array to include at most one shock per boundary point
#    used_boundaries = []
#    ind = 0
#    preserved_indices = []
#    for ind in 1:length(intersections)
#        loc, time, t1, t2 = intersections[ind]
#        if !(t1 in used_boundaries) && !(t2 in used_boundaries)
#            push!(used_boundaries, t1, t2)
#            push!(preserved_indices, ind)
#        end
#    end
#    intersections = intersections[preserved_indices]
#
#    # Sort remaining shock points into curves.
#    # Note that shock curves must be piecewise quadratic
#
#    
#    # Increment shockcounts along curves
#    # For now, just do points
#    radius = systemparams.ballsize
#    gridsize = systemparams.stepsize
#    γ = systemparams.threshold
#
#    for shock in intersections
#
#        # If the shock is sufficiently large
#        magnitude = abs(boundary_values[shock[3]] - boundary_values[shock[4]])
#        if magnitude > γ
#
#            # Then increment nearby grid points
#            x_c, t_c = shock[1:2]/gridsize
#            min_vals = Int(ceil((shock[1:2] .- radius)/gridsize))
#            max_vals = Int(floor((shock[1:2] .+ radius)/gridsize))
#            for x in min_vals[1]:max_vals[1]
#                for t in min_vals[2]:max_vals[2]
#                    if (x - x_c)^2 + (t - t_c)^2 < radius^2
#                        shockcounts[x,t] += 1
#                    end
#                end
#            end
#        end
#    end
#    
end
export IncrementShockCounts!


"""
    ComputeIntersection(t1, u1, t2, u2)

Helper function to compute the intersection between
two characteristic curves in Burgers equation
starting at x=0.
"""
function ComputeIntersection(t1, u1, t2, u2)
    if u1 == u2
        return nothing 
    elseif (t2-t1)/(u2-u1) * u1 * u2 < 0
        return nothing
    else
        return ((u2*t2 - u1*t2)/(u2-u1), (t2-t1)/(u2-u1) * u1 * u2)
    end
end

"""
Sample from RTS smoother backwards sampler for random walk
"""
function RTS_sample(estimates, variances, RTSvariances, σ_w)
    # Start by generating the noise
    out = sqrt.(RTSvariances) .* randn(length(estimates))
    
    # Add in the means
    for i in length(out)-1:-1:1
        weight = σ_w^2/(σ_w^2 + variances[i])
        out[i] += weight * estimates[i] + (1 - weight) * out[i+1]
    end

    return out
end

"""
Helper Function computing variance for RTS sampler

Returns forward estimates with RTS variances
"""
function RTS_Smooth(observations, σ_w, σ_v)
    estimates, variances = KalmanFilterWalk(observations, σ_w, σ_v)
    RTSvariances = BackwardsVarianceWalk(variances, σ_w)
    return estimates, variances, RTSvariances
end
export RTS_Smooth

"""
Helper function computing Kalman filter for the given model
"""
function KalmanFilterWalk(observations, σ_w, σ_v)
    estimates = zeros(length(observations))
    variances = zeros(length(observations))

    # Initialize the Kalman filter
    estimates[1] = observations[0]
    variances[1] = σ_v^2

    for i in 2:length(observations)
        weight = σ_v^2/(variances[i-1] + σ_v^2 + σ_w^2)
        estimates[i] = weight * estimates[i-1] + (1 - weight) * observations[i]
        variances[i] = weight * (variances[i-1] + σ_w^2)
    end

    return estimates, variances
end

"""
Helper Function computing variance from RTS smoother
"""
function BackwardsVarianceWalk(variances, σ_w)
    RTSvariances = zeros(length(variances))
    for i in length(variances)-1:-1:1
        weight = σ_w^2/(variances[i] + σ_w^2)
        RTSvariances[i] = weight * variances[i] + (1-weight)^2 * RTSvariances[i+1]
    end
    return RTSvariances
end
