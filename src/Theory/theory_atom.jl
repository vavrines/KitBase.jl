"""
    ref_vhs_vis(Kn, alpha, omega)

Calculate reference viscosity with variable hard sphere (VHS) model

"""
ref_vhs_vis(Kn, alpha, omega) =
    5.0 * (alpha + 1.0) * (alpha + 2.0) * √π /
    (4.0 * alpha * (5.0 - 2.0 * omega) * (7.0 - 2.0 * omega)) * Kn


"""
    vhs_collision_time(ρ, λ, μᵣ, ω)
    vhs_collision_time(prim::AV{T}, muRef, omega) where {T<:Real}

Calculate collision time with variable hard sphere (VHS) model

"""
vhs_collision_time(ρ, λ, μᵣ, ω) = μᵣ * 2.0 * λ^(1.0 - ω) / ρ

vhs_collision_time(prim::AV{T}, muRef, omega) where {T<:Real} =
    muRef * 2.0 * prim[end]^(1.0 - omega) / prim[1] # for rykov model prim[end] should be λₜ


"""
    aap_hs_collision_time(
        prim::AA{<:Real,2},
        mi::Real,
        ni::Real,
        me::Real,
        ne::Real,
        kn::Real,
    )

Calculate mixture collision time from AAP model
"""
function aap_hs_collision_time(prim::AM{T}, mi, ni, me, ne, kn) where {T<:Real}

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
    aap_hs_diffeq!(du, u, p, t)

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


"""
    boltzmann_ode!(df, f, p, t)

RHS-ODE of Boltzmann equation
"""
function boltzmann_ode!(df, f::AA{T,3}, p, t) where {T<:Real}
    Kn, M, phi, psi, phipsi = p
    df .= boltzmann_fft(f, Kn, M, phi, psi, phipsi)
end


"""
    boltzmann_nuode!(df, f, p, t)

RHS-ODE of Boltzmann equation with non-uniform velocity
"""
function boltzmann_nuode!(df, f::AA{T,3}, p, t) where {T<:Real}
    Kn, M, phi, psi, phipsi, u, v, w, vnu, u1, v1, w1, vuni = p

    nu = length(u)
    nv = length(v)
    nw = length(w)
    nu1 = length(u1)
    nv1 = length(v1)
    nw1 = length(w1)

    curve = itp.RegularGridInterpolator((u, v, w), f)
    _f = reshape(curve(vuni), nu1, nv1, nw1)
    _df = boltzmann_fft(_f, Kn, M, phi, psi, phipsi)

    curve1 = itp.RegularGridInterpolator((u1, v1, w1), _df)
    df .= reshape(curve1(vnu), nu, nv, nw)
end


"""
    bgk_ode!(df, f, p, t)
    
RHS-ODE of BGK equation
"""
function bgk_ode!(df, f::AA{T}, p, t) where {T<:Real}
    g, τ = p
    df .= (g .- f) ./ τ
end
