"""
$(SIGNATURES)

Bose-Einstein integral `Li_{ν+1}(z)`
"""
be_integral(ν, z) = polylog(ν + 1, z) |> real


"""
$(SIGNATURES)

Fermi-Dirac integral `-Li_{ν+1}(-z)`
"""
fd_integral(ν, z) = -polylog(ν + 1, -z) |> real


"""
$(SIGNATURES)

Equilibrium distribution of Bosons
"""
function be_equilibrium(u, prim, β)
    A, U, λ = prim
    return @. β / sqrt(π) / (A^(-1) * exp(λ * (u - U)^2) - 1)
end

"""
$(SIGNATURES)
"""
function be_equilibrium(u, v, prim, β)
    A, U, V, λ = prim
    return @. β / π / (A^(-1) * exp(λ * ((u - U)^2 + (v - V)^2)) - 1)
end


"""
$(SIGNATURES)

Equilibrium distribution of Fermions
"""
function fd_equilibrium(u, prim, β)
    A, U, λ = prim
    return @. β / sqrt(π) / (A^(-1) * exp(λ * (u - U)^2) + 1)
end

"""
$(SIGNATURES)
"""
function fd_equilibrium(u, v, prim, β)
    A, U, V, λ = prim
    return @. β / π / (A^(-1) * exp(λ * ((u - U)^2 + (v - V)^2)) + 1)
end

# ------------------------------------------------------------
# Methods for computing primitive variables
# ------------------------------------------------------------

function Aeq_1d(u, p)
    ρ, ρe, β, fn = p
    return 4 * β^2 * fn(-0.5, u)^3 * ρe - ρ^3 * fn(0.5, u)
end
const Aprob0_1d =
    NonlinearSolve.NonlinearProblem{false}(Aeq_1d, 0.3, (0.556, 0.19, 2.0, fd_integral))

function Aeq_2d(u, p)
    ρ, ρe, β, fn = p
    return 2 * β * fn(0.0, u)^2 * ρe - ρ^2 * fn(1.0, u)
end
const Aprob0_2d =
    NonlinearSolve.NonlinearProblem{false}(Aeq_2d, 0.66, (1.77, 1.75, 2.0, fd_integral))

function Aprob(w, β, it = :fd)
    fn = eval(Symbol(string(it) * "_integral"))
    if length(w) == 3
        return NonlinearSolve.remake(
            Aprob0_1d,
            u0 = w[1],
            p = (w[1], w[end] - w[2]^2 / w[1] / 2, β, fn),
        )
    else
        return NonlinearSolve.remake(
            Aprob0_2d,
            u0 = w[1],
            p = (w[1], w[end] - (w[2]^2 + w[3]^2) / w[1] / 2, β, fn),
        )
    end
end


"""
$(SIGNATURES)

Transform conservative -> primitive variables
"""
function quantum_conserve_prim(w, β, it = :fd)
    prim = zero(w)
    fn = eval(Symbol(string(it) * "_integral"))

    prob = Aprob(w, β, it)
    sol = NonlinearSolve.solve(prob, NLSolveJL(), reltol = 1e-9)
    @assert sol.u[1] != 0
    prim[1] = sol.u[1]
    prim[2] = w[2] / w[1]

    if length(w) == 3 # 1D
        prim[3] = β^2 * fn(-0.5, prim[1])^2 / w[1]^2
    elseif length(w) == 4 # 2D
        prim[3] = w[3] / w[1]
        prim[4] = β * fn(0.0, prim[1]) / w[1]
    end

    return prim
end


"""
$(SIGNATURES)

Transform primitive -> conservative variables
"""
function quantum_prim_conserve(prim, β, it = :fd)
    w = zero(prim)
    fn = eval(Symbol(string(it) * "_integral"))

    if length(prim) == 3 # 1D
        w[1] = fn(-0.5, prim[1]) * β / sqrt(prim[end])
        w[2] = w[1] * prim[2]
        w[3] = fn(0.5, prim[1]) * β / (prim[end])^(3 / 2) / 4 + 0.5 * w[1] * prim[2]^2
    elseif length(prim) == 4 # 2D
        w[1] = fn(0.0, prim[1]) * β / prim[end]
        w[2] = w[1] * prim[2]
        w[3] = w[1] * prim[2]
        w[4] =
            0.5 * fn(1.0, prim[1]) * β / (prim[end])^2 +
            0.5 * w[1] * (prim[2]^2 + prim[3]^2)
    end
    return w
end
