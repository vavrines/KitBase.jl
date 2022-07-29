"""
$(SIGNATURES)

1D1F1V
"""
function step!(
    w::T3,
    prim::T3,
    f::T4,
    fwL::T1,
    ffL::T2,
    fwR::T1,
    ffR::T2,
    u::T5,
    weights::T5,
    γ,
    μᵣ,
    ω,
    Pr,
    dx,
    dt,
    RES,
    AVG,
    collision,
) where {T1<:AV,T2<:AV,T3<:AV,T4<:AV,T5<:AV}

    #--- store W^n and calculate H^n,\tau^n ---#
    w_old = deepcopy(w)

    if collision == :shakhov
        q = heat_flux(f, prim, u, weights)
        M_old = maxwellian(u, prim)
        S = shakhov(u, M_old, q, prim, Pr)
    else
        S = zeros(axes(f))
    end

    #--- update W^{n+1} ---#
    @. w += (fwL - fwR) / dx
    prim .= conserve_prim(w, γ)

    #--- record residuals ---#
    @. RES += (w - w_old)^2
    @. AVG += abs(w)

    #--- calculate M^{n+1} and tau^{n+1} ---#
    M = maxwellian(u, prim)
    M .+= S
    τ = vhs_collision_time(prim, μᵣ, ω)

    #--- update distribution function ---#
    for i in eachindex(u)
        f[i] = (f[i] + (ffL[i] - ffR[i]) / dx + dt / τ * M[i]) / (1.0 + dt / τ)
    end

end

"""
$(SIGNATURES)

1D1F3V
"""
function step!(
    w::T3,
    prim::T3,
    f::T4,
    fwL::T1,
    ffL::T2,
    fwR::T1,
    ffR::T2,
    uVelo::T5,
    vVelo::T5,
    wVelo::T5, # avoid conflict with w
    weights::T5,
    γ,
    μᵣ,
    ω,
    Pr,
    dx,
    dt,
    RES,
    AVG,
    collision,
) where {T1<:AA{<:FN,1},T2<:AA{<:FN,3},T3<:AA{<:FN,1},T4<:AA{<:FN,3},T5<:AA{<:FN,3}}

    #--- store W^n and calculate shakhov term ---#
    w_old = deepcopy(w)

    if collision == :shakhov
        q = heat_flux(f, prim, uVelo, vVelo, wVelo, weights)
        M_old = maxwellian(uVelo, vVelo, wVelo, prim)
        S = shakhov(uVelo, vVelo, wVelo, M_old, q, prim, Pr, K)
    else
        S = zeros(axes(f))
    end

    #--- update W^{n+1} ---#
    @. w += (fwL - fwR) / dx
    prim .= conserve_prim(w, γ)

    #--- record residuals ---#
    @. RES += (w - w_old)^2
    @. AVG += abs(w)

    #--- calculate M^{n+1} and tau^{n+1} ---#
    M = maxwellian(uVelo, vVelo, wVelo, prim)
    M .+= S
    τ = vhs_collision_time(prim, μᵣ, ω)

    #--- update distribution function ---#
    for k in axes(wVelo, 3), j in axes(vVelo, 2), i in axes(uVelo, 1)
        f[i, j, k] =
            (f[i, j, k] + (ffL[i, j, k] - ffR[i, j, k]) / dx + dt / τ * M[i, j, k]) /
            (1.0 + dt / τ)
    end

end

"""
$(SIGNATURES)

1D1F3V fast spectral method
"""
function step!(
    w::T3,
    prim::T3,
    f::T4,
    fwL::T1,
    ffL::T2,
    fwR::T1,
    ffR::T2,
    γ,
    Kn_bz,
    nm,
    phi,
    psi,
    phipsi,
    dx,
    dt,
    RES,
    AVG,
    collision,
) where {T1<:AV,T2<:AA{T5,3},T3<:AV,T4<:AA{T6,3}} where {T5,T6}

    @assert collision == :fsm

    w_old = deepcopy(w)
    @. w += (fwL - fwR) / dx
    prim .= conserve_prim(w, γ)

    @. RES += (w - w_old)^2
    @. AVG += abs(w)

    Q = zero(f[:, :, :])
    boltzmann_fft!(Q, f, Kn_bz, nm, phi, psi, phipsi)

    for k in axes(f, 3), j in axes(f, 2), i in axes(f, 1)
        f[i, j, k] += (ffL[i, j, k] - ffR[i, j, k]) / dx + dt * Q[i, j, k]
    end

end

"""
$(SIGNATURES)

2D1F2V
"""
function step!(
    w::T1,
    prim::T1,
    h::T2,
    fwL::T1,
    fhL::T2,
    fwR::T1,
    fhR::T2,
    fwD::T1,
    fhD::T2,
    fwU::T1,
    fhU::T2,
    u::T3,
    v::T3,
    weights::T3,
    γ,
    μᵣ,
    ω,
    Pr,
    Δs,
    dt,
    RES,
    AVG,
    collision,
) where {T1<:AV,T2<:Union{AM,AV},T3<:Union{AM,AV}}

    #--- store W^n and calculate shakhov term ---#
    w_old = deepcopy(w)

    if collision == :shakhov
        q = heat_flux(h, b, prim, u, v, weights)

        MH_old = maxwellian(u, v, prim)
        SH = shakhov(u, v, MH_old, q, prim, Pr)
    else
        SH = zero(h)
    end

    #--- update W^{n+1} ---#
    @. w += (fwL - fwR + fwD - fwU) / Δs
    prim .= conserve_prim(w, γ)

    #--- record residuals ---#
    @. RES += (w - w_old)^2
    @. AVG += abs(w)

    #--- calculate M^{n+1} and tau^{n+1} ---#
    MH = maxwellian(u, v, prim)
    MH .+= SH
    τ = vhs_collision_time(prim, μᵣ, ω)

    #--- update distribution function ---#
    for i in eachindex(u)
        h[i] =
            (h[i] + (fhL[i] - fhR[i] + fhD[i] - fhU[i]) / Δs + dt / τ * MH[i]) /
            (1.0 + dt / τ)
    end

end

"""
$(SIGNATURES)

2D1F2V @ triangle
"""
function step!(
    w::T1,
    prim::T1,
    f::T2,
    fw1::T1,
    ff1::T2,
    fw2::T1,
    ff2::T2,
    fw3::T1,
    ff3::T2,
    u::T3,
    v::T3,
    weights::T3,
    K,
    γ,
    μᵣ,
    ω,
    Pr,
    Δs,
    dirc::AV,
    dt,
    RES,
    AVG,
    collision,
) where {T1<:AV,T2<:Union{AM,AV},T3<:Union{AM,AV}}

    #--- store W^n and calculate shakhov term ---#
    w_old = deepcopy(w)

    if collision == :shakhov
        q = heat_flux(h, b, prim, u, v, weights)

        M_old = maxwellian(u, v, prim)
        S = shakhov(u, v, M_old, q, prim, Pr, K)
    else
        S = zero(f)
    end

    #--- update W^{n+1} ---#
    @. w -= (fw1 * dirc[1] + fw2 * dirc[2] + fw3 * dirc[3]) / Δs
    prim .= conserve_prim(w, γ)

    #--- record residuals ---#
    @. RES += (w - w_old)^2
    @. AVG += abs(w)

    #--- calculate M^{n+1} and tau^{n+1} ---#
    M = maxwellian(u, v, prim)
    M .+= S
    τ = vhs_collision_time(prim, μᵣ, ω)

    #--- update distribution function ---#
    for i in eachindex(u)
        f[i] =
            (
                f[i] - (ff1[i] * dirc[1] + ff2[i] * dirc[2] + ff3[i] * dirc[3]) / Δs +
                dt / τ * M[i]
            ) / (1.0 + dt / τ)
    end

end

"""
$(SIGNATURES)

2D1F3V fast spectral method
"""
function step!(
    w::T3,
    prim::T3,
    f::T4,
    fwL::T1,
    ffL::T2,
    fwR::T1,
    ffR::T2,
    fwD::T1,
    ffD::T2,
    fwU::T1,
    ffU::T2,
    γ,
    Kn_bz,
    nm,
    phi,
    psi,
    phipsi,
    Δs,
    dt,
    RES,
    AVG,
    collision,
) where {T1<:AA{<:FN,1},T2<:AA{<:FN,3},T3<:AA{<:FN,1},T4<:AA{<:FN,3}}

    @assert collision == :fsm

    w_old = deepcopy(w)
    @. w += (fwL - fwR + fwD - fwU) / Δs

    prim .= conserve_prim(w, γ)

    @. RES += (w - w_old)^2
    @. AVG += abs(w)

    Q = zero(f[:, :, :])
    boltzmann_fft!(Q, f, Kn_bz, nm, phi, psi, phipsi)

    for k in axes(f, 3), j in axes(f, 2), i in axes(f, 1)
        f[i, j, k] +=
            (ffL[i, j, k] - ffR[i, j, k] + ffD[i, j, k] - ffU[i, j, k]) / Δs +
            dt * Q[i, j, k]
    end

    return nothing

end
