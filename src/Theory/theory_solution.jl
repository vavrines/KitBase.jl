"""
Solves the Riemann problem for the Euler equations.

## Arguments

  - `left_state`: Tuple (ρₗ, uₗ, pₗ)
  - `right_state`: Tuple (ρᵣ, uᵣ, pᵣ)
  - `x`: spatial coordinate at which to evaluate the solution
  - `t`: time at which to evaluate the solution (t ≥ 0)
  - `γ`: ratio of specific heats (default 1.4)

## Outputs

  - `(ρ, u, p)` at (x, t)
"""
function sample_riemann_solution(left_state, right_state, x, t, γ=1.4; tol=1e-6)
    ρₗ, uₗ, pₗ = left_state
    ρᵣ, uᵣ, pᵣ = right_state

    # At t = 0, return the initial condition
    t == 0 && return x < 0 ? left_state : right_state

    # Sound speeds in left and right states
    cₗ = sqrt(γ * pₗ / ρₗ)
    cᵣ = sqrt(γ * pᵣ / ρᵣ)

    # f-function for a given state (shock or rarefaction)
    function f(p, ρ, pᵢ, c)
        if p > pᵢ  # shock wave
            A = 2.0 / ((γ + 1) * ρ)
            B = (γ - 1) / (γ + 1) * pᵢ
            return (p - pᵢ) * sqrt(A / (p + B))
        else  # rarefaction wave
            return 2 * c / (γ - 1) * ((p / pᵢ)^((γ - 1) / (2γ)) - 1)
        end
    end

    # Derivative of f-function with respect to p
    function df(p, ρ, pᵢ, c)
        if p > pᵢ  # shock
            A = 2.0 / ((γ + 1) * ρ)
            B = (γ - 1) / (γ + 1) * pᵢ
            sqrt_term = sqrt(A / (p + B))
            return sqrt_term * (1 - 0.5 * (p - pᵢ) / (p + B))
        else  # rarefaction
            return (1 / (ρ * c)) * (p / pᵢ)^(-(γ + 1) / (2γ))
        end
    end

    # Find p_star using Newton's method
    p_old = 0.5 * (pₗ + pᵣ)  # initial guess
    for _ in 1:20
        fₗ = f(p_old, ρₗ, pₗ, cₗ)
        fᵣ = f(p_old, ρᵣ, pᵣ, cᵣ)
        func = fₗ + fᵣ + (uᵣ - uₗ)
        dfₗ = df(p_old, ρₗ, pₗ, cₗ)
        dfᵣ = df(p_old, ρᵣ, pᵣ, cᵣ)
        dp = -func / (dfₗ + dfᵣ)
        p_new = p_old + dp
        p_new = max(p_new, 1e-6)  # enforce positive pressure
        abs(dp) < tol && break
        p_old = p_new
    end
    p_star = p_old
    u_star = 0.5 * (uₗ + uᵣ) + 0.5 * (f(p_star, ρᵣ, pᵣ, cᵣ) - f(p_star, ρₗ, pₗ, cₗ))

    # Compute wave speeds
    # Left wave speeds
    if p_star > pₗ  # left shock
        Sₗ = uₗ - cₗ * sqrt((γ + 1) / (2γ) * (p_star / pₗ) + (γ - 1) / (2γ))
    else  # left rarefaction
        Sₕₗ = uₗ - cₗ               # head of rarefaction
        Sₜₗ = u_star - cₗ * (p_star / pₗ)^((γ - 1) / (2γ))  # tail of rarefaction
    end

    # Right wave speeds
    if p_star > pᵣ  # right shock
        Sᵣ = uᵣ + cᵣ * sqrt((γ + 1) / (2γ) * (p_star / pᵣ) + (γ - 1) / (2γ))
    else  # right rarefaction
        Sₕᵣ = uᵣ + cᵣ               # head of rarefaction
        Sₜᵣ = u_star + cᵣ * (p_star / pᵣ)^((γ - 1) / (2γ))  # tail of rarefaction
    end

    ξ = x / t  # self-similar variable

    # Determine the solution at (t, x) based on ξ
    if ξ < u_star  # solution from the left side
        if p_star > pₗ  # left shock
            if ξ < Sₗ
                return left_state
            else
                # state in the left star region (behind the shock)
                ρ_star_l =
                    ρₗ * ((p_star / pₗ) + (γ - 1) / (γ + 1)) /
                    (((γ - 1) / (γ + 1)) * (p_star / pₗ) + 1)
                return (ρ_star_l, u_star, p_star)
            end
        else  # left rarefaction
            if ξ < Sₕₗ
                return left_state
            elseif ξ > Sₜₗ
                # state in the left star region (after the fan)
                ρ_star_l = ρₗ * (p_star / pₗ)^(1 / γ)
                return (ρ_star_l, u_star, p_star)
            else
                # within the rarefaction fan
                u = (2 / (γ + 1)) * (cₗ + (γ - 1) / 2 * uₗ + ξ)
                c = (2 / (γ + 1)) * (cₗ + (γ - 1) / 2 * (uₗ - ξ))
                ρ = ρₗ * (c / cₗ)^(2 / (γ - 1))
                p = pₗ * (c / cₗ)^(2γ / (γ - 1))
                return (ρ, u, p)
            end
        end
    else  # solution from the right side
        if p_star > pᵣ  # right shock
            if ξ > Sᵣ
                return right_state
            else
                ρ_star_r =
                    ρᵣ * ((p_star / pᵣ) + (γ - 1) / (γ + 1)) /
                    (((γ - 1) / (γ + 1)) * (p_star / pᵣ) + 1)
                return (ρ_star_r, u_star, p_star)
            end
        else  # right rarefaction
            if ξ > Sₕᵣ
                return right_state
            elseif ξ < Sₜᵣ
                # state in the right star region (before the fan ends)
                ρ_star_r = ρᵣ * (p_star / pᵣ)^(1 / γ)
                return (ρ_star_r, u_star, p_star)
            else
                # within the rarefaction fan on the right
                u = (2 / (γ + 1)) * (-cᵣ + (γ - 1) / 2 * uᵣ + ξ)
                c = (2 / (γ + 1)) * (cᵣ - (γ - 1) / 2 * (ξ - uᵣ))
                ρ = ρᵣ * (c / cᵣ)^(2 / (γ - 1))
                p = pᵣ * (c / cᵣ)^(2γ / (γ - 1))
                return (ρ, u, p)
            end
        end
    end
end

function sample_riemann_solution(left_state, right_state, x::AV, t, γ=1.4; tol=1e-6)
    return [
        sample_riemann_solution(left_state, right_state, x[i], t, γ; tol=tol) for
        i in 1:length(x)
    ]
end
