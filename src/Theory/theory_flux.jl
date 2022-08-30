"""
$(SIGNATURES)

Theoretical flux of linear advection equation
"""
advection_flux(u, a) = a * u


"""
$(SIGNATURES)

Theoretical flux of Burgers' equation
"""
burgers_flux(u) = 0.5 * u^2


"""
$(SIGNATURES)

Theoretical fluxes of Euler Equations
"""
function euler_flux(w::AV, γ; frame = :cartesian::Symbol)
    prim = conserve_prim(w, γ)
    p = 0.5 * prim[1] / prim[end]

    if eltype(w) <: Integer
        F = similar(w, Float64)
        G = similar(w, Float64)
        H = similar(w, Float64)
    else
        F = similar(w)
        G = similar(w)
        H = similar(w)
    end

    if length(w) == 3
        F[1] = w[2]
        F[2] = w[2]^2 / w[1] + p
        F[3] = (w[3] + p) * w[2] / w[1]

        return (F,)
    elseif length(w) == 4
        F[1] = w[2]
        F[2] = w[2]^2 / w[1] + p
        F[3] = w[2] * w[3] / w[1]
        F[4] = (w[end] + p) * w[2] / w[1]

        G[1] = w[3]
        G[2] = w[3] * w[2] / w[1]
        G[3] = w[3]^2 / w[1] + p
        G[4] = (w[end] + p) * w[3] / w[1]

        return F, G
    elseif length(w) == 5
        F[1] = w[2]
        F[2] = w[2]^2 / w[1] + p
        F[3] = w[2] * w[3] / w[1]
        F[4] = w[2] * w[4] / w[1]
        F[5] = (w[end] + p) * w[2] / w[1]

        G[1] = w[3]
        G[2] = w[3] * w[2] / w[1]
        G[3] = w[3]^2 / w[1] + p
        G[4] = w[3] * w[4] / w[1]
        G[5] = (w[end] + p) * w[3] / w[1]

        H[1] = w[4]
        H[2] = w[4] * w[2] / w[1]
        H[3] = w[4] * w[3] / w[1]
        H[4] = w[4]^2 / w[1] + p
        H[5] = (w[end] + p) * w[4] / w[1]

        return F, G, H
    end
end


"""
$(SIGNATURES)

Flux Jacobian of Euler Equations
"""
function euler_jacobi(w::AV, γ)
    if eltype(w) <: Integer
        A = similar(w, Float64, 3, 3)
    else
        A = similar(w, 3, 3)
    end

    A[1, 1] = 0.0
    A[1, 2] = 1.0
    A[1, 3] = 0.0

    A[2, 1] = 0.5 * (γ - 3.0) * (w[2] / w[1])^2
    A[2, 2] = (3.0 - γ) * w[2] / w[1]
    A[2, 3] = γ - 1.0

    A[3, 1] = -γ * w[3] * w[2] / w[1]^2 + (γ - 1.0) * (w[2] / w[1])^3
    A[3, 2] = γ * w[3] / w[1] - 1.5 * (γ - 1.0) * (w[2] / w[1])^2
    A[3, 3] = γ * w[2] / w[1]

    return A
end
