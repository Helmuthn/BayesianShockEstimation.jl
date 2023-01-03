using SpecialFunctions: gamma
using Interpolations: linear_interpolation
using LinearAlgebra: dot, norm
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
    radius = systemparams.ballsize
    threshold = systemparams.threshold
    
    # Initialize a collection of ODEs evolving over space
    t_lst = Array(T * (1:0.1:length(boundary_values)))
    x_lst = zeros(length(t_lst))
    boundary = linear_interpolation(1:length(boundary_values), boundary_values)
    u_lst = boundary.(1:0.1:length(boundary_values))

    shock_curves = []
    shock_curve_indices = zeros(length(u_lst))

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

        # For each shock found in the last search
        for shock in removed

            # Compute intersection (curve knot)
            # TODO
            intersect = (1,1)
            
            # Get the current shock curve
            curve_ind = max(shock_curve_indices[shock[1]], shock_curve_indices[shock[2]])

            # If the curve is new, create a new curve
            if curve_ind == 0
                append!(shock_curves, [])
                curve_ind = length(shock_curves)
            end

            # Append the knot to the curve, as well as the slope for Rankine-Hugoniot
            append!(shock_curves[curve_ind], (intersect, u_lst[shock[2]] - u_lst[shock[1]]))

            # Set neighboring labels
            if shock[1] > 1
                shock_curve_indices[shock[1]-1] = curve_ind
            end 
            if shock[2] < length(shock_curve_indices)
                shock_curve_indices[shock[2] + 1] = curve_ind
            end
        end

        for trajectory in introduced
            # Introduce new ODE trajectories in the center.
            
        end

    end

    # Delete the shocks
    # deleteat!()
    
    # Add the new trajectories
    # pushat!()

    # Identify all of the indices to increment
    increment_lst = FindPointsFromCurves(shock_curves, radius)

    # Increment the indices
    for loc in increment_lst
        shockcounts[loc[1], loc[2]] += 1
    end
end
export IncrementShockCounts!

"""
Helper function that converts curves to indices to increment
"""
function FindPointsFromCurves(shock_curves, radius)
    increment_lst = []
    for curve in shock_curves

        # For each pair of knots, get all indices in the associated interval
        for i in 1:length(curve)-1
            knot1 = curve[i]
            knot2 = curve[i+1]

            x_indices = []
            t_indices = []
            for x in x_indices
                for t in t_indices
                    if distanceToSegment(knot1, knot2, [x,t]) < radius
                        append!(increment_lst, [x,t])
                    end
                end 
            end
        end
    end

    return unique(increment_lst)
end

"""
Helper function computing distance from p3
to the line segment defined by the other points.
"""
function distanceToSegment(p1, p2, p3)
    segment = p2-p1
    dot_prod = dot(segment, p3)

    if dot_prod > eps()
        return norm(p2-p3)
    elseif dot_prod < eps()
        return norm(p1-p3)
    else
        num = abs((p2[1]-p1[1]) * (p1[2] - p3[2]) - (p1[1] - p3[1]) * (p2[2] - p1[1]))
        return num/norm(p2-p1) 
    end
end

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


#############################################
###### Burgers Equation Solver ##############
#############################################

"""
    flux_burgers(u)

Helper function returns flux in from the conservation form of Burgers' Equation.

# Arguments
 - `u` - The state value
"""
function flux_burgers(u)
    return 0.5 * u^2
end

"""
    riemann_solver_burgers(ul, ur)

Solves the Riemann problem for Burgers' equation and returns the average value.
Helper function for Godunov's scheme.

# Arguments
 - `ul` - The left state value
 - `ur  - The right state value

# Returns
The value of flux on the boundary between the regions.
"""
function riemann_solver_burgers(ul,ur)
    # Characteristic curves of ul intersect boundary
    if ul >= 0 && ur >= 0
        return flux_burgers(ul)
    
    # Characteristic curves of ur intersect boundary    
    elseif ul <= 0 && ur <= 0
        return flux_burgers(ur)
        
    # Rarefaction with zero speed curve
    elseif ul <= 0 && ur >= 0
        return 0.0
    
    # Shock direction needed
    elseif ul >= 0 && ur <= 0
        
        # Shock Travels Right
        if ul + ur > 0
            return flux_burgers(ul)
            
        # Shock Travels Left
        else
            return flux_burgers(ur)
        end
    end
end

"""
    riemann_solver_burgers!(boundary_flux_right, u)

Solve for all right boundary flux values.

# Arguments
 - `boundary_flux_right` - Output Array
 - `u` - State Values
"""
function riemann_solver_burgers!(boundary_flux_right, u)
    ul = @view u[1:end-1]
    ur = @view u[2:end]
    boundary_flux_right[1:end-1] .= riemann_solver_burgers.(ul, ur)
    boundary_flux_right[end] = riemann_solver_burgers(u[end],u[1])
end

"""
    flux_difference!(boundary_flux, boundary_flux_right)

Compute the difference between the right and left flux

# Arguments
 - `boundary_flux` - output array
 - `boundary_flux_right` - flux through right boundaries
"""
function flux_difference!(boundary_flux, boundary_flux_right)
    boundary_flux_right_t = @view boundary_flux_right[2:end]
    boundary_flux_left = @view boundary_flux_right[1:end-1]
    boundary_flux[2:end] .= boundary_flux_right_t .- boundary_flux_left
    boundary_flux[1] = boundary_flux_right[1] - boundary_flux_right[end]
end

"""
    godunov_burgers_1D_step!(u_next, boundary_flux, boundary_flux_right, u, dx, dt)

Step Burgers' equation according to Godunov's method

# Arguments
 - `u_next` - Array to store the value
 - `boundary_flux` - Array of boundary flux storage
 - `boundary_flux_right` - Array of right boundary flux storage
 - `u`  - Array of values for the previous timestep
 - `dx` - Spatial grid size
 - `dt` - Temporal grid size

# Returns
An array of values for the next timestep

# Notes
Assumes periodic boundary conditions
"""
function godunov_burgers_1D_step!(u_next, boundary_flux, boundary_flux_right, u, λ )
    # Get flux values
    riemann_solver_burgers!(boundary_flux_right, u)
    flux_difference!(boundary_flux, boundary_flux_right)
    boundary_flux .*= λ
    
    u_next .= u
    u_next .-= boundary_flux
end

"""
    godunov_burgers_1D(u0, dx, dt, stepcount)

Solve Burgers' equation numericallly with Godunov's method

# Arguments
 - `u0`        - Initial Boundary Condition
 - `dx`        - Spatial grid size
 - `dt`        - Temporal grid size
 - `stepcount` - Number of steps

# Returns
A 2D array representing the numerical solution of Burgers' Equation
"""
function godunov_burgers_1D(u0, dx, dt, stepcount)
    # Allocate flux memory
    boundary_flux = zeros(length(u0))
    boundary_flux_right = zeros(length(u0))
    
    # Allocate the output array
    out = zeros(length(u0), stepcount)
    out[:,1] .= u0
    
    λ = dt/dx

    for i in range(2,stepcount)
        u_next = @view out[:,i]
        u      = @view out[:,i-1]
        godunov_burgers_1D_step!(u_next, boundary_flux, boundary_flux_right, u, λ)
    end

    return out
end
export godunov_burgers_1D
