# ============================================================
# Continuum Theory
# ============================================================

"""
$(SIGNATURES)

Transform primitive -> conservative variables
"""
function prim_conserve!(W::AV, prim::AV, γ)
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

    return nothing
end

"""
$(SIGNATURES)

Polyatomic gas

* w: ρ, ρV⃗, ρEₜ, ρEᵣ
* prim: ρ, V⃗, λ, λₜ, λᵣ
"""
function prim_conserve!(W::AV, prim::AV, γ, Kr)
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
    elseif length(prim) == 7 # 3D
        W[1] = prim[1]
        W[2] = prim[1] * prim[2]
        W[3] = prim[1] * prim[3]
        W[4] = prim[1] * prim[4]
        W[5] =
            0.5 * prim[1] / prim[5] / (γ - 1.0) +
            0.5 * prim[1] * (prim[2]^2 + prim[3]^2 + prim[4]^2)
        W[6] = 0.25 * prim[1] * Kr / prim[7]
    else
        throw("prim -> w : dimension error")
    end

    return W
end


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
    prim_conserve!(W, prim, γ)

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

Polyatomic gas
"""
function prim_conserve(prim::AV, γ, Kr)
    if eltype(prim) <: Integer
        W = similar(prim, Float64, length(prim) - 1)
    else
        W = similar(prim, length(prim) - 1)
    end
    prim_conserve!(W, prim, γ, Kr)

    return W
end


"""
$(SIGNATURES)

Transform multi-component primitive -> conservative variables
"""
function mixture_prim_conserve!(w, prim::AM, γ)
    for j in axes(w, 2)
        _w = @view w[:, j]
        _prim = @view prim[:, j]
        prim_conserve!(_w, _prim, γ)
    end

    return w
end

"""
$(SIGNATURES)

Transform multi-component polyatomic primitive -> conservative variables
"""
function mixture_prim_conserve!(w, prim::AM, γ, Kr)
    for j in axes(w, 2)
        _w = @view w[:, j]
        _prim = @view prim[:, j]
        prim_conserve!(_w, _prim, γ, Kr)
    end

    return w
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
    mixture_prim_conserve!(w, prim, γ)

    return w
end

"""
$(SIGNATURES)

Transform multi-component primitive -> conservative variables
"""
function mixture_prim_conserve(prim::AM, γ, Kr)
    if eltype(prim) <: Integer
        w = similar(prim, Float64, size(prim, 1) - 1, size(prim, 2))
    else
        w = similar(prim, size(prim, 1) - 1, size(prim, 2))
    end
    mixture_prim_conserve!(w, prim, γ, Kr)

    return w
end


"""
$(SIGNATURES)

vector: primitive vector for Euler, Navier-Stokes and extended equations
"""
function conserve_prim!(prim, W::AV, γ)
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

    return nothing
end

"""
$(SIGNATURES)

Polyatomic gas

* w: ρ, ρV⃗, ρEₜ, ρEᵣ
* prim: ρ, V⃗, λ, λₜ, λᵣ
"""
function conserve_prim!(prim, w::AV, K, Kr)
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
    elseif length(w) == 6 # 3D
        prim[1] = w[1]
        prim[2] = w[2] / w[1]
        prim[3] = w[3] / w[1]
        prim[4] = w[4] / w[1]
        prim[5] =
            0.25 * w[1] * (K + Kr + 3.0) / (w[5] - 0.5 * (w[2]^2 + w[3]^2 + w[4]^2) / w[1])
        prim[6] =
            0.25 * w[1] * (K + 3.0) /
            (w[5] - w[6] - 0.5 * (w[2]^2 + w[3]^2 + w[4]^2) / w[1])
        prim[7] = 0.25 * w[1] * Kr / w[6]
    else
        throw("w -> prim : dimension dismatch")
    end

    return nothing
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
    conserve_prim!(prim, W, γ)

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
    conserve_prim!(prim, w, K, Kr)

    return prim
end


"""
$(SIGNATURES)

Transform multi-component conservative -> primitive variables
"""
function mixture_conserve_prim!(prim, W::AM, γ)
    for j in axes(prim, 2)
        _prim = @view prim[:, j]
        _w = @view W[:, j]
        conserve_prim!(_prim, _w, γ)
    end

    return nothing
end

"""
$(SIGNATURES)

Transform multi-component polyatomic conservative -> primitive variables
"""
function mixture_conserve_prim!(prim, W::AM, K, Kr)
    for j in axes(prim, 2)
        _prim = @view prim[:, j]
        _w = @view W[:, j]
        conserve_prim!(_prim, _w, K, Kr)
    end

    return nothing
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
    mixture_conserve_prim!(prim, W, γ)

    return prim
end

"""
$(SIGNATURES)

Transform multi-component polyatomic conservative -> primitive variables
"""
function mixture_conserve_prim(W::AM, K, Kr)
    if eltype(W) <: Integer
        prim = similar(W, Float64, size(W, 1) + 1, size(W, 2))
    else
        prim = similar(W, size(W, 1) + 1, size(W, 2))
    end
    mixture_conserve_prim!(prim, W, K, Kr)

    return prim
end


"""
$(SIGNATURES)

Calculate heat capacity ratio (monatomic gas)
"""
function heat_capacity_ratio(K, D::Integer)
    if D == 1
        return (K + 3.0) / (K + 1.0)
    elseif D == 2
        return (K + 4.0) / (K + 2.0)
    elseif D == 3
        return (K + 5.0) / (K + 3.0)
    end
end

"""
$(SIGNATURES)

Calculate heat capacity ratio (diatomic gas)
"""
function heat_capacity_ratio(K, Nr, D::Integer)
    if D == 1
        return (K + 3.0 + Nr) / (K + 1.0 + Nr)
    elseif D == 2
        return (K + 4.0 + Nr) / (K + 2.0 + Nr)
    elseif D == 3
        return (K + 5.0 + Nr) / (K + 3.0 + Nr)
    end
end


"""
$(SIGNATURES)

Calculate internal degrees of freedom (monatomic gas)
"""
function internal_dof(γ, D::Integer)
    if D == 1
        return (3 - γ) / (γ - 1)
    elseif D == 2
        return (4 - 2 * γ) / (γ - 1)
    elseif D == 3
        return (5 - 3 * γ) / (γ - 1)
    end
end


"""
$(SIGNATURES)

Calculate speed of sound
"""
sound_speed(λ::Real, γ::Real) = (0.5 * γ / λ)^0.5

"""
$(SIGNATURES)
"""
sound_speed(prim::AV, γ) = sound_speed(prim[end], γ)

"""
$(SIGNATURES)

Calculate sound speed in mixture
"""
function sound_speed(prim::AM, γ)
    c = similar(prim, axes(prim, 2))
    for j in eachindex(c)
        c[j] = sound_speed(prim[end, j], γ)
    end

    return maximum(c)
end


"""
$(SIGNATURES)

Calculate reference viscosity with variable hard sphere (VHS) model
"""
ref_vhs_vis(Kn, alpha, omega) =
    5.0 * (alpha + 1.0) * (alpha + 2.0) * √π /
    (4.0 * alpha * (5.0 - 2.0 * omega) * (7.0 - 2.0 * omega)) * Kn

# ------------------------------------------------------------
# Entropy
# ------------------------------------------------------------

"""
$(SIGNATURES)

Maxwell Boltzmann entropy
"""
maxwell_boltzmann(f) = f * log(f) - f


"""
$(SIGNATURES)

Prim of Maxwell Boltzmann entropy
"""
maxwell_boltzmann_prime(x) = log(x)


"""
$(SIGNATURES)

Dual of Maxwell Boltzmann entropy
"""
maxwell_boltzmann_dual(f) = exp(f)


"""
$(SIGNATURES)

Dual prim of Maxwell Boltzmann entropy
"""
maxwell_boltzmann_dual_prime(f) = exp(f)


"""
$(SIGNATURES)

Reconstruct mathematical entropy from Legendre dual
"""
function kinetic_entropy(α::AA, m::AA, weights::AV)
    B = KitBase.maxwell_boltzmann_dual_prime.(α' * m)[:]
    return sum(maxwell_boltzmann.(B) .* weights)
end


"""
$(SIGNATURES)

Calculate entropy from fluid variables
"""
fluid_entropy(ρ, c, γ) = c^2 / γ / ρ^(γ - 1.0)

function fluid_entropy(prim::AV, γ)
    c = sound_speed(prim, γ)
    return fluid_entropy(prim[1], c, γ)
end
