"""
$(SIGNATURES)
    
RHS-ODE of BGK equation
"""
function bgk_ode!(df, f::AA{T}, p, t) where {T}
    M, τ = p
    df .= (M .- f) ./ τ
end


"""
$(SIGNATURES)
    
RHS-ODE of ES-BGK equation
"""
function esbgk_ode!(df, f::AA{T}, p, t) where {T}
    gs = p[1:end-1]
    τ = p[end]
    G = esbgk(f, gs...)
    df .= (G .- f) ./ τ
end


"""
$(SIGNATURES)

Calculate collision time with variable hard sphere (VHS) model
"""
vhs_collision_time(ρ, λ, μᵣ, ω) = μᵣ * 2.0 * λ^(1.0 - ω) / ρ

"""
$(SIGNATURES)

For Rykov model prim[end] should be λₜ
"""
vhs_collision_time(prim::AV, muRef, omega) =
    vhs_collision_time(prim[1], prim[end], muRef, omega)

# ------------------------------------------------------------
# ν model
# ------------------------------------------------------------

"""
$(SIGNATURES)

Regulator to construct ν-BGK model

_L. Mieussens and H. Struchtrup. Phys. Fluids, 16(8): 2797-2813, 2004._
"""
function νbgk_relaxation_time(τ::T, u, prim::AV, tp)::T where {T}
    c = @. abs(u - prim[2])
    return @. τ / mieussens_frequency(c, tp)
end

function νbgk_relaxation_time(τ::T, u, v, prim::AV, tp)::T where {T}
    c = @. sqrt((u - prim[2])^2 + (v - prim[3])^2)
    return @. τ / mieussens_frequency(c, tp)
end

function νbgk_relaxation_time(τ::T, u, v, w, prim::AV, tp)::T where {T}
    c = @. sqrt((u - prim[2])^2 + (v - prim[3])^2 + (w - prim[4])^2)
    return @. τ / mieussens_frequency(c, tp)
end

mieussens_frequency(η, ::Type{Class{1}}) = @. 0.431587 * η^1.791288

mieussens_frequency(η, ::Type{Class{2}}) = @. 0.0268351 * (1 + 14.2724 * η^2)

mieussens_frequency(η, ::Type{Class{3}}) = @. 0.0365643 * (1 + 10 * η^2.081754)

mieussens_frequency(η, ::Type{Class{4}}) = @. 0.1503991 * (1 + 0.92897 * η^4) # usually this isn't used


"""
$(SIGNATURES)

Regulator to construct ν-Shakhov model

_Yuan & Wu, J. Fluid Mech., 942, 2022._

``ν = A(p/μ)(νeq0(ξ)+2νeq0(0))``
"""
function νshakhov_relaxation_time(τ::T, u, prim::AV)::T where {T}
    c = @. abs(u - prim[2])
    return @. 1 / 0.0871 * τ / (νeq0(c) + 2 * νeq0(eps()))
end

νeq0(x) = 1.5 * (exp(-x^2) + sqrt(pi) / 2 * (1 / x + 2 * x) * erf(x))

# ------------------------------------------------------------
# AAP model
# ------------------------------------------------------------

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

Calculate mixture collision time from AAP model
"""
function aap_hs_collision_time(prim::AM, mi, ni, me, ne, kn)
    ν = similar(prim, 2)

    ν[1] =
        prim[1, 1] / (mi * (ni + ne)) * 4.0 * sqrt(π) / 3.0 *
        sqrt(1.0 / prim[end, 1] + 1.0 / prim[end, 1]) / (sqrt(2.0) * π * kn) +
        prim[1, 2] / (me * (ni + ne)) * 4.0 * sqrt(π) / 3.0 *
        sqrt(1.0 / prim[end, 1] + 1.0 / prim[end, 2]) / (sqrt(2.0) * π * kn)
    ν[2] =
        prim[1, 1] / (mi * (ni + ne)) * 4.0 * sqrt(π) / 3.0 *
        sqrt(1.0 / prim[end, 1] + 1.0 / prim[end, 2]) / (sqrt(2.0) * π * kn) +
        prim[1, 2] / (me * (ni + ne)) * 4.0 * sqrt(π) / 3.0 *
        sqrt(1.0 / prim[end, 2] + 1.0 / prim[end, 2]) / (sqrt(2.0) * π * kn)

    return 1.0 ./ ν
end


"""
$(SIGNATURES)

Source term of AAP model in DifferentialEquations.jl
"""
function aap_hs_diffeq!(du, u, p, t)
    if length(u) == 6
        I₁, I₂, I₃, E₁, E₂, E₃ = u
        w = [
            I₁ E₁
            I₂ E₂
            I₃ E₃
        ]
    elseif length(u) == 8
        I₁, I₂, I₃, I₄, E₁, E₂, E₃, E₄ = u
        w = [
            I₁ E₁
            I₂ E₂
            I₃ E₃
            I₄ E₄
        ]
    elseif length(u) == 10
        I₁, I₂, I₃, I₄, I₅, E₁, E₂, E₃, E₄, E₅ = u
        w = [
            I₁ E₁
            I₂ E₂
            I₃ E₃
            I₄ E₄
            I₅ E₅
        ]
    else
    end

    τᵢ, τₑ, mi, ni, me, ne, kn, γ = p
    τ = [τᵢ, τₑ]

    # modified variables
    prim = mixture_conserve_prim(w, γ)
    mixprim = aap_hs_prim(prim, τ, mi, ni, me, ne, kn)
    mixw = mixture_conserve_prim(mixprim, γ)

    if length(u) == 6
        du[1] = (mixw[1, 1] - I₁) / τᵢ
        du[2] = (mixw[2, 1] - I₂) / τᵢ
        du[3] = (mixw[3, 1] - I₃) / τᵢ
        du[4] = (mixw[1, 2] - E₁) / τₑ
        du[5] = (mixw[2, 2] - E₂) / τₑ
        du[6] = (mixw[3, 2] - E₃) / τₑ
    elseif length(u) == 8
        du[1] = (mixw[1, 1] - I₁) / τᵢ
        du[2] = (mixw[2, 1] - I₂) / τᵢ
        du[3] = (mixw[3, 1] - I₃) / τᵢ
        du[4] = (mixw[4, 1] - I₄) / τᵢ
        du[5] = (mixw[1, 2] - E₁) / τₑ
        du[6] = (mixw[2, 2] - E₂) / τₑ
        du[7] = (mixw[3, 2] - E₃) / τₑ
        du[8] = (mixw[4, 2] - E₄) / τₑ
    elseif length(u) == 10
        du[1] = (mixw[1, 1] - I₁) / τᵢ
        du[2] = (mixw[2, 1] - I₂) / τᵢ
        du[3] = (mixw[3, 1] - I₃) / τᵢ
        du[4] = (mixw[4, 1] - I₄) / τᵢ
        du[5] = (mixw[5, 1] - I₅) / τᵢ
        du[6] = (mixw[1, 2] - E₁) / τₑ
        du[7] = (mixw[2, 2] - E₂) / τₑ
        du[8] = (mixw[3, 2] - E₃) / τₑ
        du[9] = (mixw[4, 2] - E₄) / τₑ
        du[10] = (mixw[5, 2] - E₅) / τₑ
    end

    return nothing
end
