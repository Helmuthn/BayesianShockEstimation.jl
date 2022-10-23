using SpecialFunctions: gamma
# This file contains the core contributions of the work.

"""
    ShockParams

Struct containing all need information for deterministic shock simulations

# Parameters
 - `dynamics`: Representation of the dynamics
 - `positions`: List of sample positions
 - `ballsize`: Integration ball size
 - `threshold`: Minimum shock threshold
 - `dimension`: System dimensionality

"""
struct ShockParams
    dynamics
    positions
    ballsize
    threshold
    dimension::Int
end
export ShockParams

"""
    BoundaryParams

Struct containing all needed information for boundary sampling

# Parameters
 - `μ` : Mean at sample times
 - `σ` : Standard deviation at sample times
 - `t` : Sample times

"""
struct BoundaryParams
    μ::Array{Float64}
    σ::Array{Float64}
    t::Array{Float64}
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
    
    # Monte Carlo Integration Based On Stochastic Boundaries
    for i in range(samples)
        boundary = SampleBoundary(boundaryparams::BoundaryParams)
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
 - `boundary`: Boundary conditions being tested

# Returns
A list of approximate densities

# Notes
Uses the method of characteristics to find shocks.
"""
function IncrementShockCounts!(shockcounts, systemparams::ShockParams, boundary)
    # Find all pairwise intersection times
    intersect_locs = zeros(100,100)
    remaining_curves = 1:100

    shocks = []
    # Find closest interesctions, increment balls, then remove curves from list
    


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
    SampleBoundary(boundaryparams)

Generate a sample of a boundary.

# Arguments
  - `boundaryparams`: Parameters characterizing the distribution of the boundary

# Returns
A function representing a sample of boundary conditions.
"""
function SampleBoundary(boundaryparams::BoundaryParams)
    sample_count = length(boundaryparams.t)
    σ = boundaryparams.σ
    μ = boundaryparams.μ
    return σ .* randn(sample_count) + μ
end
export SampleBoundary

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
    RTSvariances = BackwardsCovarianceWalk(variances, σ_w)
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

function BackwardsCovarianceWalk(variances, σ_w)
    RTSvariances = zeros(length(variances))
    for i in length(variances)-1:-1:1
        weight = σ_w^2/(variances[i] + σ_w^2)
        RTSvariances[i] = weight * variances[i] + (1-weight)^2 * RTSvariances[i+1]
    end
    return RTSvariances
end

"""
Helper Function computing backwards smoothing according to RTS
"""
function BackwardsSample(estimates, variances, σ_w)
    out = zeros()
end
