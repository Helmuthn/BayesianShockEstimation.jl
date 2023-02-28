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

    shockcounts = zeros(length(observations), systemparams.stepcount)

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
    normalized_counts = shockcounts ./ systemparams.ballsize ./ 2

    # Convert to a density by normalizing to the number of samples
    shockdensity = normalized_counts ./ samplecount
    
    return shockdensity 
end
export ShockDensity



"""
    IncrementShockCounts!(shockcounts, systemparams, boundary)

Computes shocks for a given boundary.

# Arguments
 - `shockcounts`     : Array of current shock counts
 - `systemparams`    : Struct of parameters for system simulation
 - `boundary_values` : Boundary samples

# Returns
A list of approximate densities
"""
function IncrementShockCounts!(shockcounts, systemparams::ShockParams, boundary_values)
    # Unpack the struct
    dx        = systemparams.stepsize_x
    dt        = systemparams.stepsize_t
    stepcount = systemparams.stepcount
    threshold = systemparams.threshold
    radius    = systemparams.ballsize

    solution = godunov_burgers_1D(boundary_values, dx, dt, stepcount)

    offsets = getoffsets(dx, dt, radius)
    jumps = rangefilter(solution, offsets)

    shockcounts .+= jumps .> threshold
end
export IncrementShockCounts!

"""
    getoffsets(dx, dt, radius)

Gets the offset indices based on step sizes
and ball radius.

# Arguments
 - `dx`     : Spatial Stepsize
 - `dt`     : Temporal Stepsize
 - `radius` : Search Radius

# Returns
An array of tuples representing index offsets
"""
function getoffsets(dx, dt, radius)
    max_x = Int(floor(radius/dx))
    max_t = Int(floor(radius/dt))
    offsets = []
    for i in 0:max_x
        for j in 0:max_t
            if (i * dx)^2 + (j * dt)^2 < radius^2
                push!(offsets, (i,j))
                if i != 0
                    push!(offsets, (-i,j))
                end
                if j != 0
                    push!(offsets, (i,-j))
                end
                if i !=0 && j != 0
                    push!(offsets, (-i,-j))
                end
            end
        end
    end
    return offsets
end

"""
    rangefilter(image, offsets)

Computes the difference between the max filter
and the min filter applied to an image on a window
defined by offsets.

# Arguments
 - `image`   : The image to be filtered, 2D array
 - `offsets` : Array of tuples representing index offsets

# Returns
A 2D array representing the result of the range filter
"""
function rangefilter(image, offsets)
    M, N = size(image)
    out_max = copy(image)
    out_min = copy(image)
    for offset in offsets
        x_range = 1:M
        y_range = 1:N
        if offset[1] > 0
            x_range = (1 + offset[1]):M
        end
        if offset[1] < 0
            x_range = 1:(M+offset[1])
        end
        if offset[2] > 0
            y_range = (1 + offset[2]):N
        end
        if offset[2] < 0
            y_range = 1:(N+offset[2])
        end


        shifted_image = view(image, x_range .- offset[1],
                                    y_range .- offset[2])

        out_max_view = view(out_max, x_range, y_range)
        out_min_view = view(out_min, x_range, y_range)
        
        out_max_view .= max.(out_max_view, shifted_image)
        out_min_view .= min.(out_min_view, shifted_image)
    end
    return out_max .- out_min
end


