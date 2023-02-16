using SpecialFunctions: gamma
using Interpolations: linear_interpolation
using LinearAlgebra: dot, norm
# This file contains the core contributions of the work.

"""
    ShockParams

Struct containing all needed information for deterministic shock simulations

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

    shockcounts = zeros(length(observations)-1, systemparams.stepcount-1)

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
    normalized_counts = shockcounts ./ systemparams.ballsize

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
    radius    = systemparams.ballsize

    shocks = godunov_burgers_1D_shocks(boundary_values, dx, dt, stepcount, threshold)

    shock_indicator = copy(shocks)
    offsets = getoffsets(dx, dt, radius)
    for offset in offsets
        CheckShiftedShock!(shock_indicator, shocks, offset)
    end
    

    shockcounts .+= shock_indicator
end
export IncrementShockCounts!

function getoffsets(dx, dt, radius)
    max_x = Int(floor(radius/dx))
    max_t = Int(floor(radius/dt))
    offsets = []
    for i in 0:max_x
        for j in 0:max_t
            if (i * dx)^2 + (j * dt)^2 < 1
                push!(offsets, (i,j))
                push!(offsets, (-i,j))
                push!(offsets, (i,-j))
                push!(offsets, (-i,-j))
            end
        end
    end
    return offsets
end

function CheckShiftedShock!(shock_indicator, shocks, offset)
    if offset[1] > 0
        x_indices = 1 + offset[1]:size(shocks)[1]
    else 
        x_indices = 1:size(shocks)[1] + offset[1]
    end
    if offset[2] > 0
        t_indices = 1+offset[2]:size(shocks)[2]
    else
        t_indices = 1:size(shocks)[2] + offset[2]
    end

    shifted_block  = @view shocks[x_indices, t_indices]

    if offset[1] > 0
        x_indices = 1:size(shocks)[1] - offset[1]
    else 
        x_indices = 1-offset[1]:size(shocks)[1]
    end
    if offset[2] > 0
        t_indices = 1:size(shocks)[2] - offset[2]
    else
        t_indices = 1-offset[2]:size(shocks)[2]
    end

    original_block = @view shock_indicator[x_indices, t_indices]

    original_block .|= shifted_block
end
