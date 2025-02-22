using Plots

# Convert primitive variables (ρ, u, p) to conservative variables (ρ, ρ*u, E)
function primitive_to_conservative(ρ, u, p, γ)
    E = p / (γ - 1) + 0.5 * ρ * u^2
    return [ρ, ρ * u, E]
end

# Convert conservative variables (ρ, ρ*u, E) to primitive variables (ρ, u, p)
function conservative_to_primitive(U, γ)
    ρ = U[1]
    u = U[2] / ρ
    E = U[3]
    p = (γ - 1) * (E - 0.5 * ρ * u^2)
    return ρ, u, p
end

# Compute the flux F(U) for a given state U
function compute_flux(U, γ)
    ρ, u, p = conservative_to_primitive(U, γ)
    return [ρ * u, ρ * u^2 + p, u * (U[3] + p)]
end

# Compute the sound speed
function sos(p, ρ, γ)
    return sqrt(γ * p / ρ)
end

# Exact Riemann solver for the Euler equations
function exact_riemann_solver(UL, UR, γ)
    # Convert left/right states from conservative to primitive variables
    ρL, uL, pL = conservative_to_primitive(UL, γ)
    ρR, uR, pR = conservative_to_primitive(UR, γ)
    cL = sos(pL, ρL, γ)
    cR = sos(pR, ρR, γ)

    # Define helper functions f(p) and its derivative for the Newton iteration
    f(p, p_i, ρ_i, c_i) =
        p > p_i ? (p - p_i) * sqrt(2.0 / ((γ + 1) * ρ_i) / (p + (γ - 1) / (γ + 1) * p_i)) :
        (2 * c_i / (γ - 1)) * ((p / p_i)^((γ - 1) / (2γ)) - 1)
    df(p, p_i, ρ_i, c_i) =
        p > p_i ?
        sqrt(2.0 / ((γ + 1) * ρ_i) / (p + (γ - 1) / (γ + 1) * p_i)) *
        (1 - 0.5 * (p - p_i) / (p + (γ - 1) / (γ + 1) * p_i)) :
        (1 / (ρ_i * c_i)) * ((p / p_i)^(-((γ + 1) / (2γ))))

    # Initial guess for p_star (using the PVRS approximate solver)
    p_guess = 0.5 * (pL + pR) - 0.125 * (uR - uL) * (ρL + ρR) * (cL + cR)
    p_guess = max(1e-6, p_guess)

    # Newton iteration to compute p_star
    p_old = p_guess
    for _ in 1:100
        fL = f(p_old, pL, ρL, cL)
        fR = f(p_old, pR, ρR, cR)
        f_total = fL + fR + (uR - uL)
        dfL = df(p_old, pL, ρL, cL)
        dfR = df(p_old, pR, ρR, cR)
        dp = -f_total / (dfL + dfR)
        p_new = p_old + dp
        if abs(dp) < 1e-6 * p_old
            p_old = p_new
            break
        end
        p_old = p_new
    end
    p_star = p_old

    # Compute u_star.
    u_star = 0.5 * (uL + uR) + 0.5 * (f(p_star, pR, ρR, cR) - f(p_star, pL, ρL, cL))

    # Sample the state at x/t = 0:
    if u_star >= 0
        # Left star state
        if p_star > pL
            ρ_star =
                ρL *
                ((p_star / pL + (γ - 1) / (γ + 1)) / ((γ - 1) / (γ + 1) * p_star / pL + 1))
        else
            ρ_star = ρL * (p_star / pL)^(1 / γ)
        end
        U_star = primitive_to_conservative(ρ_star, u_star, p_star, γ)
    else
        # Right star state
        if p_star > pR
            ρ_star =
                ρR *
                ((p_star / pR + (γ - 1) / (γ + 1)) / ((γ - 1) / (γ + 1) * p_star / pR + 1))
        else
            ρ_star = ρR * (p_star / pR)^(1 / γ)
        end
        U_star = primitive_to_conservative(ρ_star, u_star, p_star, γ)
    end

    return compute_flux(U_star, γ)
end

# Godunov solver using the exact Riemann solver
function godunov_solver(num_cells, x_min, x_max, t_final, γ)
    dx = (x_max - x_min) / num_cells
    # Define cell centers
    x = [x_min + (i - 0.5) * dx for i in 1:num_cells]

    # Initialize solution using the Sod shock tube initial conditions
    U = Array{Float64}(undef, num_cells, 3)
    for i in 1:num_cells
        if x[i] < 0.5 * (x_min + x_max)
            ρ, u, p = 1.0, 0.0, 1.0
        else
            ρ, u, p = 0.125, 0.0, 0.1
        end
        U[i, :] = primitive_to_conservative(ρ, u, p, γ)
    end

    t = 0.0
    CFL = 0.9

    while t < t_final
        # Determine time step from CFL condition.
        max_speed = 0.0
        for i in 1:num_cells
            ρ, u, p = conservative_to_primitive(U[i, :], γ)
            c = sos(p, ρ, γ)
            max_speed = max(max_speed, abs(u) + c)
        end
        dt = CFL * dx / max_speed
        if t + dt > t_final
            dt = t_final - t
        end

        # Compute fluxes at cell interfaces using the exact Riemann solver.
        fluxes = Array{Float64}(undef, num_cells + 1, 3)
        for i in 2:num_cells
            fluxes[i, :] = exact_riemann_solver(U[i-1, :], U[i, :], γ)
        end
        # Boundary interfaces (transmissive BCs)
        fluxes[1, :] = compute_flux(U[1, :], γ)
        fluxes[num_cells+1, :] = compute_flux(U[num_cells, :], γ)

        # Update solution via the finite volume update.
        U_new = copy(U)
        for i in 1:num_cells
            U_new[i, :] .= U[i, :] .- dt / dx * (fluxes[i+1, :] .- fluxes[i, :])
        end

        U = U_new
        t += dt
    end
    return x, U
end

# Main execution
function main()
    γ = 1.4
    num_cells = 400
    x_min = 0.0
    x_max = 1.0
    t_final = 0.2  # final time

    x, U = godunov_solver(num_cells, x_min, x_max, t_final, γ)

    # Extract primitive variables for plotting.
    ρ_arr = zeros(length(x))
    u_arr = zeros(length(x))
    p_arr = zeros(length(x))
    for i in 1:length(x)
        ρ, u, p = conservative_to_primitive(U[i, :], γ)
        ρ_arr[i] = ρ
        u_arr[i] = u
        p_arr[i] = p
    end

    # Plot the results.
    plt1 = plot(x, ρ_arr; label="Density", xlabel="x", ylabel="Density", legend=:top)
    plt2 = plot(
        x,
        u_arr;
        label="Velocity",
        xlabel="x",
        ylabel="Velocity",
        legend=:top,
        color=:red,
    )
    plt3 = plot(
        x,
        p_arr;
        label="Pressure",
        xlabel="x",
        ylabel="Pressure",
        legend=:top,
        color=:green,
    )
    return plot(plt1, plt2, plt3; layout=(3, 1), size=(800, 600))
end

main()
