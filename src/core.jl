module core
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
"""
struct BoundaryParams

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
    
end
export IncrementShockCounts!

"""
    StepCharacteristics
"""
function StepCharacteristics()

end

"""
    SampleBoundary(boundaryparams)

Generate a sample of a boundary.
"""
function SampleBoundary(boundaryparams::BoundaryParams)
    
end
export SampleBoundary
end
