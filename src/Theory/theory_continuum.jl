# ============================================================
# Continuum Theory
# ============================================================

"""
$(SIGNATURES)

Transform primitive -> conservative variables
"""
function prim_conserve(prim::AV, γ)
    if eltype(prim) <: Integer
        W = similar(prim, Float64)
    else
        W = similar(prim)
    end

    if length(prim) == 3 # 1D
        W[1] = prim[1]
        W[2] = prim[1] * prim[2]
        W[3] = 0.5 * prim[1] / prim[3] / (γ - 1.0) + 0.5 * prim[1] * prim[2]^2
    elseif length(prim) == 4 # 2D
        W[1] = prim[1]
        W[2] = prim[1] * prim[2]
        W[3] = prim[1] * prim[3]
        W[4] = 0.5 * prim[1] / prim[4] / (γ - 1.0) + 0.5 * prim[1] * (prim[2]^2 + prim[3]^2)
    elseif length(prim) == 5 # 3D
        W[1] = prim[1]
        W[2] = prim[1] * prim[2]
        W[3] = prim[1] * prim[3]
        W[4] = prim[1] * prim[4]
        W[5] =
            0.5 * prim[1] / prim[5] / (γ - 1.0) +
            0.5 * prim[1] * (prim[2]^2 + prim[3]^2 + prim[4]^2)
    else
        throw("prim -> w : dimension error")
    end

    return W
end

"""
$(SIGNATURES)
"""
prim_conserve(ρ, U, λ, γ) = prim_conserve([ρ, U, λ], γ)

"""
$(SIGNATURES)
"""
prim_conserve(ρ, U, V, λ, γ) = prim_conserve([ρ, U, V, λ], γ)

"""
$(SIGNATURES)
"""
prim_conserve(ρ, U, V, W, λ, γ) = prim_conserve([ρ, U, V, W, λ], γ)

"""
$(SIGNATURES)

Rykov model
"""
function prim_conserve(prim::AV, γ, Kr)
    if eltype(prim) <: Integer
        W = similar(prim, Float64, length(prim) - 1)
    else
        W = similar(prim, length(prim) - 1)
    end

    if length(prim) == 5 # 1D
        W[1] = prim[1]
        W[2] = prim[1] * prim[2]
        W[3] = 0.5 * prim[1] / prim[3] / (γ - 1.0) + 0.5 * prim[1] * prim[2]^2
        W[4] = 0.25 * prim[1] * Kr / prim[5]
    elseif length(prim) == 6 # 2D
        W[1] = prim[1]
        W[2] = prim[1] * prim[2]
        W[3] = prim[1] * prim[3]
        W[4] = 0.5 * prim[1] / prim[4] / (γ - 1.0) + 0.5 * prim[1] * (prim[2]^2 + prim[3]^2)
        W[5] = 0.25 * prim[1] * Kr / prim[6]
    else
        throw("prim -> w : dimension error")
    end

    return W
end


"""
$(SIGNATURES)

Transform multi-component primitive -> conservative variables
"""
function mixture_prim_conserve(prim::AM, γ)
    if eltype(prim) <: Integer
        w = similar(prim, Float64)
    else
        w = similar(prim)
    end

    for j in axes(w, 2)
        w[:, j] .= prim_conserve(prim[:, j], γ)
    end

    return w
end


"""
$(SIGNATURES)

Transform conservative -> primitive variables

scalar: pseudo primitive vector for scalar conservation laws
"""
conserve_prim(u) = [u, 0.5 * u, 1.0]

"""
$(SIGNATURES)
"""
conserve_prim(u, a) = [u, a, 1.0]

"""
$(SIGNATURES)

vector: primitive vector for Euler, Navier-Stokes and extended equations
"""
function conserve_prim(W::AV, γ)
    if eltype(W) <: Integer
        prim = similar(W, Float64)
    else
        prim = similar(W)
    end

    if length(W) == 3 # 1D
        prim[1] = W[1]
        prim[2] = W[2] / W[1]
        prim[3] = 0.5 * W[1] / (γ - 1.0) / (W[3] - 0.5 * W[2]^2 / W[1])
    elseif length(W) == 4 # 2D
        prim[1] = W[1]
        prim[2] = W[2] / W[1]
        prim[3] = W[3] / W[1]
        prim[4] = 0.5 * W[1] / (γ - 1.0) / (W[4] - 0.5 * (W[2]^2 + W[3]^2) / W[1])
    elseif length(W) == 5 # 3D
        prim[1] = W[1]
        prim[2] = W[2] / W[1]
        prim[3] = W[3] / W[1]
        prim[4] = W[4] / W[1]
        prim[5] = 0.5 * W[1] / (γ - 1.0) / (W[5] - 0.5 * (W[2]^2 + W[3]^2 + W[4]^2) / W[1])
    else
        throw("w -> prim : dimension dismatch")
    end

    return prim
end

"""
$(SIGNATURES)
"""
conserve_prim(ρ, M, E, γ) = conserve_prim([ρ, M, E], γ)

"""
$(SIGNATURES)
"""
conserve_prim(ρ, MX, MY, E, γ) = conserve_prim([ρ, MX, MY, E], γ)

"""
$(SIGNATURES)
"""
conserve_prim(ρ, MX, MY, MZ, E, γ) = conserve_prim([ρ, MX, MY, MZ, E], γ)

"""
$(SIGNATURES)

Rykov model
"""
function conserve_prim(w::AV, K, Kr)
    if eltype(w) <: Integer
        prim = similar(w, Float64, length(w) + 1)
    else
        prim = similar(w, length(w) + 1)
    end

    if length(w) == 4 # 1D
        prim[1] = w[1]
        prim[2] = w[2] / w[1]
        prim[3] = 0.25 * w[1] * (K + Kr + 1.0) / (w[3] - 0.5 * w[2]^2 / w[1])
        prim[4] = 0.25 * w[1] * (K + 1.0) / (w[3] - w[4] - 0.5 * w[2]^2 / w[1])
        prim[5] = 0.25 * w[1] * Kr / w[4]
    elseif length(w) == 5 # 2D
        prim[1] = w[1]
        prim[2] = w[2] / w[1]
        prim[3] = w[3] / w[1]
        prim[4] = 0.25 * w[1] * (K + Kr + 2.0) / (w[4] - 0.5 * (w[2]^2 + w[3]^2) / w[1])
        prim[5] = 0.25 * w[1] * (K + 2.0) / (w[4] - w[5] - 0.5 * (w[2]^2 + w[3]^2) / w[1])
        prim[6] = 0.25 * w[1] * Kr / w[5]
    else
        throw("w -> prim : dimension dismatch")
    end

    return prim
end


"""
$(SIGNATURES)

Transform multi-component conservative -> primitive variables
"""
function mixture_conserve_prim(W::AM, γ)
    if eltype(W) <: Integer
        prim = similar(W, Float64)
    else
        prim = similar(W)
    end

    for j in axes(prim, 2)
        prim[:, j] .= conserve_prim(W[:, j], γ)
    end

    return prim
end


"""
$(SIGNATURES)

Calculate mixture primitive variables from AAP model
"""
function aap_hs_prim(prim::AM, tau::AV, mi, ni, me, ne, kn)

    mixprim = similar(prim)

    if size(prim, 1) == 3

        mixprim[1, :] = deepcopy(prim[1, :])
        mixprim[2, 1] =
            prim[2, 1] +
            tau[1] / kn *
            (
                4.0 * sqrt(2.0) / (3.0 * sqrt(π)) * prim[1, 2] / (ni + ne) / (mi + me) *
                sqrt(1.0 / prim[end, 1] + 1.0 / prim[end, 2])
            ) *
            (prim[2, 2] - prim[2, 1])
        mixprim[2, 2] =
            prim[2, 2] +
            tau[2] / kn *
            (
                4.0 * sqrt(2.0) / (3.0 * sqrt(π)) * prim[1, 1] / (ni + ne) / (mi + me) *
                sqrt(1.0 / prim[end, 1] + 1.0 / prim[end, 2])
            ) *
            (prim[2, 1] - prim[2, 2])
        mixprim[3, 1] =
            1.0 / (
                1.0 / prim[end, 1] - 2.0 / 3.0 * (mixprim[2, 1] - prim[2, 1])^2 +
                tau[1] / kn * 2.0 * mi / (mi + me) *
                (
                    4.0 * sqrt(2.0) / (3.0 * sqrt(π)) * prim[1, 2] / (ni + ne) / (mi + me) *
                    sqrt(1.0 / prim[end, 1] + 1.0 / prim[end, 2])
                ) *
                (
                    1.0 / prim[end, 2] * me / mi - 1.0 / prim[end, 1] +
                    2.0 / 3.0 * me / mi * (prim[2, 2] - prim[2, 1])^2
                )
            )
        mixprim[3, 2] =
            1.0 / (
                1.0 / prim[end, 2] - 2.0 / 3.0 * (mixprim[2, 2] - prim[2, 2])^2 +
                tau[2] / kn * 2.0 * me / (mi + me) *
                (
                    4.0 * sqrt(2.0) / (3.0 * sqrt(π)) * prim[1, 1] / (ni + ne) / (mi + me) *
                    sqrt(1.0 / prim[end, 1] + 1.0 / prim[end, 2])
                ) *
                (
                    1.0 / prim[end, 1] * mi / me - 1.0 / prim[end, 2] +
                    2.0 / 3.0 * mi / me * (prim[2, 1] - prim[2, 2])^2
                )
            )

    elseif size(prim, 1) == 4

        mixprim[1, :] = deepcopy(prim[1, :])
        mixprim[2, 1] =
            prim[2, 1] +
            tau[1] / kn *
            (
                4.0 * sqrt(2.0) / (3.0 * sqrt(π)) * prim[1, 2] / (ni + ne) / (mi + me) *
                sqrt(1.0 / prim[end, 1] + 1.0 / prim[end, 2])
            ) *
            (prim[2, 2] - prim[2, 1])
        mixprim[2, 2] =
            prim[2, 2] +
            tau[2] / kn *
            (
                4.0 * sqrt(2.0) / (3.0 * sqrt(π)) * prim[1, 1] / (ni + ne) / (mi + me) *
                sqrt(1.0 / prim[end, 1] + 1.0 / prim[end, 2])
            ) *
            (prim[2, 1] - prim[2, 2])
        mixprim[3, 1] =
            prim[3, 1] +
            tau[1] / kn *
            (
                4.0 * sqrt(2.0) / (3.0 * sqrt(π)) * prim[1, 2] / (ni + ne) / (mi + me) *
                sqrt(1.0 / prim[end, 1] + 1.0 / prim[end, 2])
            ) *
            (prim[3, 2] - prim[3, 1])
        mixprim[3, 2] =
            prim[3, 2] +
            tau[2] / kn *
            (
                4.0 * sqrt(2.0) / (3.0 * sqrt(π)) * prim[1, 1] / (ni + ne) / (mi + me) *
                sqrt(1.0 / prim[end, 1] + 1.0 / prim[end, 2])
            ) *
            (prim[3, 1] - prim[3, 2])
        mixprim[4, 1] =
            1.0 / (
                1.0 / prim[end, 1] - 2.0 / 3.0 * (mixprim[2, 1] - prim[2, 1])^2 -
                2.0 / 3.0 * (mixprim[3, 1] - prim[3, 1])^2 +
                tau[1] / kn * 2.0 * mi / (mi + me) *
                (
                    4.0 * sqrt(2.0) / (3.0 * sqrt(π)) * prim[1, 2] / (ni + ne) / (mi + me) *
                    sqrt(1.0 / prim[end, 1] + 1.0 / prim[end, 2])
                ) *
                (
                    1.0 / prim[end, 2] * me / mi - 1.0 / prim[end, 1] +
                    2.0 / 3.0 * me / mi * (prim[2, 2] - prim[2, 1])^2 +
                    2.0 / 3.0 * me / mi * (prim[3, 2] - prim[3, 1])^2
                )
            )
        mixprim[4, 2] =
            1.0 / (
                1.0 / prim[end, 2] - 2.0 / 3.0 * (mixprim[2, 2] - prim[2, 2])^2 -
                2.0 / 3.0 * (mixprim[3, 2] - prim[3, 2])^2 +
                tau[2] / kn * 2.0 * me / (mi + me) *
                (
                    4.0 * sqrt(2.0) / (3.0 * sqrt(π)) * prim[1, 1] / (ni + ne) / (mi + me) *
                    sqrt(1.0 / prim[end, 1] + 1.0 / prim[end, 2])
                ) *
                (
                    1.0 / prim[end, 1] * mi / me - 1.0 / prim[end, 2] +
                    2.0 / 3.0 * mi / me * (prim[2, 1] - prim[2, 2])^2 +
                    2.0 / 3.0 * mi / me * (prim[3, 1] - prim[3, 2])^2
                )
            )

    elseif size(prim, 1) == 5

        mixprim[1, :] = deepcopy(prim[1, :])
        mixprim[2, 1] =
            prim[2, 1] +
            tau[1] / kn *
            (
                4.0 * sqrt(2.0) / (3.0 * sqrt(π)) * prim[1, 2] / (ni + ne) / (mi + me) *
                sqrt(1.0 / prim[end, 1] + 1.0 / prim[end, 2])
            ) *
            (prim[2, 2] - prim[2, 1])
        mixprim[2, 2] =
            prim[2, 2] +
            tau[2] / kn *
            (
                4.0 * sqrt(2.0) / (3.0 * sqrt(π)) * prim[1, 1] / (ni + ne) / (mi + me) *
                sqrt(1.0 / prim[end, 1] + 1.0 / prim[end, 2])
            ) *
            (prim[2, 1] - prim[2, 2])
        mixprim[3, 1] =
            prim[3, 1] +
            tau[1] / kn *
            (
                4.0 * sqrt(2.0) / (3.0 * sqrt(π)) * prim[1, 2] / (ni + ne) / (mi + me) *
                sqrt(1.0 / prim[end, 1] + 1.0 / prim[end, 2])
            ) *
            (prim[3, 2] - prim[3, 1])
        mixprim[3, 2] =
            prim[3, 2] +
            tau[2] / kn *
            (
                4.0 * sqrt(2.0) / (3.0 * sqrt(π)) * prim[1, 1] / (ni + ne) / (mi + me) *
                sqrt(1.0 / prim[end, 1] + 1.0 / prim[end, 2])
            ) *
            (prim[3, 1] - prim[3, 2])
        mixprim[4, 1] =
            prim[4, 1] +
            tau[1] / kn *
            (
                4.0 * sqrt(2.0) / (3.0 * sqrt(π)) * prim[1, 2] / (ni + ne) / (mi + me) *
                sqrt(1.0 / prim[end, 1] + 1.0 / prim[end, 2])
            ) *
            (prim[4, 2] - prim[4, 1])
        mixprim[4, 2] =
            prim[4, 2] +
            tau[2] / kn *
            (
                4.0 * sqrt(2.0) / (3.0 * sqrt(π)) * prim[1, 1] / (ni + ne) / (mi + me) *
                sqrt(1.0 / prim[end, 1] + 1.0 / prim[end, 2])
            ) *
            (prim[4, 1] - prim[4, 2])
        mixprim[5, 1] =
            1.0 / (
                1.0 / prim[end, 1] - 2.0 / 3.0 * (mixprim[2, 1] - prim[2, 1])^2 -
                2.0 / 3.0 * (mixprim[3, 1] - prim[3, 1])^2 -
                2.0 / 3.0 * (mixprim[4, 1] - prim[4, 1])^2 +
                tau[1] / kn * 2.0 * mi / (mi + me) *
                (
                    4.0 * sqrt(2.0) / (3.0 * sqrt(π)) * prim[1, 2] / (ni + ne) / (mi + me) *
                    sqrt(1.0 / prim[end, 1] + 1.0 / prim[end, 2])
                ) *
                (
                    1.0 / prim[end, 2] * me / mi - 1.0 / prim[end, 1] +
                    2.0 / 3.0 * me / mi * (prim[2, 2] - prim[2, 1])^2 +
                    2.0 / 3.0 * me / mi * (prim[3, 2] - prim[3, 1])^2 +
                    2.0 / 3.0 * me / mi * (prim[4, 2] - prim[4, 1])^2
                )
            )
        mixprim[5, 2] =
            1.0 / (
                1.0 / prim[end, 2] - 2.0 / 3.0 * (mixprim[2, 2] - prim[2, 2])^2 -
                2.0 / 3.0 * (mixprim[3, 2] - prim[3, 2])^2 -
                2.0 / 3.0 * (mixprim[4, 2] - prim[4, 2])^2 +
                tau[2] / kn * 2.0 * me / (mi + me) *
                (
                    4.0 * sqrt(2.0) / (3.0 * sqrt(π)) * prim[1, 1] / (ni + ne) / (mi + me) *
                    sqrt(1.0 / prim[end, 1] + 1.0 / prim[end, 2])
                ) *
                (
                    1.0 / prim[end, 1] * mi / me - 1.0 / prim[end, 2] +
                    2.0 / 3.0 * mi / me * (prim[2, 1] - prim[2, 2])^2 +
                    2.0 / 3.0 * mi / me * (prim[3, 1] - prim[3, 2])^2 +
                    2.0 / 3.0 * mi / me * (prim[4, 1] - prim[4, 2])^2
                )
            )

    else

        throw("AAP mixture: dimension dismatch")

    end

    return mixprim

end


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
