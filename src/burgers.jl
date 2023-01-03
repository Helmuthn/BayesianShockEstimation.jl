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

"""
    godunov_burgers_1D_shocks(u0, dx, dt, stepcount, threshold)

Finds spatial steps in the numericla solution of Burgers equation
with a slope greater than the threshold.
These points represent our approximations of shocks.

# Arguments
 - `u0`        - Initial Boundary Condition
 - `dx`        - Spatial grid size
 - `dt`        - Temporal grid size
 - `stepcount` - Number of steps
 - `threshold` - slope threshold
"""
function godunov_burgers_1D_shocks(u0, dx, dt, stepcount, threshold)
    solution = godunov_burgers_1D(u0, dx, dt, stepcount)
    shocks = zeros(Bool, (length(u0)-1,stepcount-1))
    
    threshold_x = threshold * dx
    threshold_t = threshold * dt

    # Check for spatial shocks
    solution_right = @view solution[2:end,:]
    solution_left  = @view solution[1:end-1,:]
    shocks .|= abs.(solution_right .- solution_left) .> threshold_x
    
    # Check for temporal shocks
    solution_up   = @view solution[:, 2:end]
    solution_down = @view solution[:, 1:end-1]
    shocks .|= abs.(solution_up .- solution_down) .> threshold_t

    return shocks
end
export godunov_burgers_1D_shocks
