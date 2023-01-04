using SpecialFunctions: gamma
using Interpolations: linear_interpolation
using LinearAlgebra: dot, norm
# This file contains the core contributions of the work.

"""
    ShockParams

Struct containing all need information for deterministic shock simulations

# Parameters
 - `stepsize_x`: Spatial grid size
 - `stepsize_t`: Temporal grid size
 - `stepcount` : Number of temporal steps
 - `ballsize`  : Integration ball size
 - `threshold` : Minimum shock threshold

"""
struct ShockParams
    stepsize_x
    stepsize_t
    stepcount
    ballsize
    threshold
end
export ShockParams

"""
    BoundaryParams

Struct containing all needed information for boundary sampling

# Parameters
 - `σ_w` : System noise standard deviation
 - `σ_v` : Measurement noise standard deviation

"""
struct BoundaryParams
    σ_w
    σ_v
end
export BoundaryParams

"""
    ShockDensity(  systemparams::ShockParams, 
                   boundaryparams::BoundaryParams, 
                   observations, 
                   samplecount)

Generates an approximation of the shock density function.

# Arguments
 - `systemparams`   : Struct of parameters for system simulation
 - `boundaryparams` : Struct of parameters for boundary simulation
 - `observations`   : Boundary observations
 - `samplecount`    : Number of Monte Carlo samples

# Returns
A 2D array of approximate shockwave densities
"""
function ShockDensity(  systemparams::ShockParams, 
                        boundaryparams::BoundaryParams, 
                        observations, 
                        samplecount)

    shockcounts = zeros(length(observations)-1, length(systemparams.stepcount)-1)

    # Compute RTS Smoothing Details
    σ_w = boundaryparams.σ_w
    σ_v = boundaryparams.σ_v
    estimates, variances, RTSvariances = RTS_Smooth(observations, σ_w, σ_v)
    
    # Monte Carlo Integration Based On Stochastic Boundaries
    for _ in 1:samplecount
        boundary = RTS_sample(estimates, variances, RTSvariances, σ_w)
        IncrementShockCounts!(shockcounts, systemparams, boundary)
    end

    # Normalize counts to the size of the ϵ-ball
    volume = BallVolume(systemparams)
    normalized_counts = shockcounts ./ volume

    # Convert to a density by normalizing to the number of samples
    shockdensity = normalized_counts ./ samplecount
    
    return shockdensity 
end
export ShockDensity


"""
    BallVolume(systemparams::ShockParams)::Float64

Helper function that computes the volume of the integration ball
"""
function BallVolume(systemparams::ShockParams)::Float64
    dim = 2
    R = systemparams.ballsize
    return π^(dim/2) * R^dim / gamma(dim/2 + 1) 
end


"""
    IncrementShockCounts!(shockcounts, systemparams, boundary)

Computes shocks for a given boundary.

# Arguments
 - `shockcounts`     : Array of current shock counts
 - `systemparams`    : Struct of parameters for system simulation
 - `boundary_values` : Boundary samples

# Returns
A list of approximate densities

# Notes
Uses the method of characteristics to find shocks.
"""
function IncrementShockCounts!(shockcounts, systemparams::ShockParams, boundary_values)
    # Unpack the struct
    dx        = systemparams.stepsize_x
    dt        = systemparams.stepsize_t
    stepcount = systemparams.stepcount
    threshold = systemparams.threshold

    shocks = godunov_burgers_1D_shocks(boundary_values, dx, dt, stepcount, threshold)

    shockcounts .+= shocks
end
export IncrementShockCounts!

