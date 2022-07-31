"""
$(SIGNATURES)

1D2F1V
"""
function step!(
    w::T3,
    prim::T3,
    h::T4,
    b::T4,
    fwL::T1,
    fhL::T2,
    fbL::T2,
    fwR::T1,
    fhR::T2,
    fbR::T2,
    u::T5,
    weights::T5,
    K,
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

    #--- store W^n and calculate shakhov term ---#
    w_old = deepcopy(w)

    if collision == :shakhov
        q = heat_flux(h, b, prim, u, weights)
        MH_old = maxwellian(u, prim)
        MB_old = MH_old .* K ./ (2.0 * prim[end])
        SH, SB = shakhov(u, MH_old, MB_old, q, prim, Pr, K)
    else
        SH = zeros(axes(h))
        SB = zeros(axes(b))
    end

    #--- update W^{n+1} ---#
    @. w += (fwL - fwR) / dx
    prim .= conserve_prim(w, γ)

    #--- record residuals ---#
    @. RES += (w - w_old)^2
    @. AVG += abs(w)

    #--- calculate M^{n+1} and tau^{n+1} ---#
    MH = maxwellian(u, prim)
    MB = MH .* K ./ (2.0 * prim[end])
    MH .+= SH
    MB .+= SB
    τ = vhs_collision_time(prim, μᵣ, ω)

    #--- update distribution function ---#
    for i in eachindex(u)
        h[i] = (h[i] + (fhL[i] - fhR[i]) / dx + dt / τ * MH[i]) / (1.0 + dt / τ)
        b[i] = (b[i] + (fbL[i] - fbR[i]) / dx + dt / τ * MB[i]) / (1.0 + dt / τ)
    end

end

"""
$(SIGNATURES)

1D2F1V mixture
"""
function step!(
    w::T3,
    prim::T3,
    h::T4,
    b::T4,
    fwL::T1,
    fhL::T2,
    fbL::T2,
    fwR::T1,
    fhR::T2,
    fbR::T2,
    u::T5,
    weights::T5,
    inK,
    γ,
    mi,
    ni,
    me,
    ne,
    Kn,
    Pr,
    dx,
    dt,
    RES,
    AVG,
    collision,
) where {T1<:AM,T2<:AM,T3<:AM,T4<:AM,T5<:AM}

    #--- update conservative flow variables ---#
    # w^n
    w_old = deepcopy(w)
    prim_old = deepcopy(prim)

    # flux -> w^{n+1}
    @. w += (fwL - fwR) / dx
    prim .= mixture_conserve_prim(w, γ)

    # temperature protection
    if prim[end, 1] < 0
        @warn "negative temperature update of component 1"
        w .= w_old
        prim .= prim_old
    elseif prim[end, 2] < 0
        @warn "negative temperature update of component 2"
        w .= w_old
        prim .= prim_old
    end

    # source -> w^{n+1}
    #=
    # DifferentialEquations.jl
    tau = get_tau(prim, mi, ni, me, ne, Kn)
    for j in axes(w, 2)
        prob = ODEProblem(aap_hs_diffeq!,
            vcat(w[1:end,j,1], w[1:end,j,2]),
            dt,
            (tau[1], tau[2], mi, ni, me, ne, Kn, γ)
        )
        sol = solve(prob, Rosenbrock23())

        w[:,j,1] .= sol[end][1:end÷2]
        w[:,j,2] .= sol[end][end÷2+1:end]
    end
    prim .= mixture_conserve_prim(w, γ)
    =#
    # explicit
    tau = aap_hs_collision_time(prim, mi, ni, me, ne, Kn)
    mprim = aap_hs_prim(prim, tau, mi, ni, me, ne, Kn)
    mw = mixture_prim_conserve(mprim, γ)
    for k in axes(w, 2)
        @. w[:, k] += (mw[:, k] - w_old[:, k]) * dt / tau[k]
    end
    prim .= mixture_conserve_prim(w, γ)

    #--- update particle distribution function ---#
    # flux -> f^{n+1}
    @. h += (fhL - fhR) / dx
    @. b += (fbL - fbR) / dx

    # source -> f^{n+1}
    tau = aap_hs_collision_time(prim, mi, ni, me, ne, Kn)

    # interspecies interaction
    #mprim = deepcopy(prim)
    mprim = aap_hs_prim(prim, tau, mi, ni, me, ne, Kn)

    H = mixture_maxwellian(u, mprim)
    B = similar(H)
    for j in axes(B, 2)
        B[:, j] = H[:, j] * inK / (2.0 * mprim[end, j])
    end

    # BGK term
    for k in axes(h, 2)
        @. h[:, k] = (h[:, k] + dt / tau[k] * H[:, k]) / (1.0 + dt / tau[k])
        @. b[:, k] = (b[:, k] + dt / tau[k] * B[:, k]) / (1.0 + dt / tau[k])
    end

    #--- record residuals ---#
    @. RES += (w_old - w)^2
    @. AVG += abs(w)

end

"""
$(SIGNATURES)

2D2F2V
"""
function step!(
    w::T1,
    prim::T1,
    h::T2,
    b::T2,
    fwL::T1,
    fhL::T2,
    fbL::T2,
    fwR::T1,
    fhR::T2,
    fbR::T2,
    fwD::T1,
    fhD::T2,
    fbD::T2,
    fwU::T1,
    fhU::T2,
    fbU::T2,
    u::T3,
    v::T3,
    weights::T3,
    K,
    γ,
    μᵣ,
    ω,
    Pr,
    Δs,
    dt,
    RES,
    AVG,
    collision,
) where {T1<:AV,T2<:AVOM,T3<:AVOM}

    #--- store W^n and calculate shakhov term ---#
    w_old = deepcopy(w)

    if collision == :shakhov
        q = heat_flux(h, b, prim, u, v, weights)

        MH_old = maxwellian(u, v, prim)
        MB_old = MH_old .* K ./ (2.0 * prim[end])
        SH, SB = shakhov(u, v, MH_old, MB_old, q, prim, Pr, K)
    else
        SH = zero(h)
        SB = zero(b)
    end

    #--- update W^{n+1} ---#
    @. w += (fwL - fwR + fwD - fwU) / Δs
    prim .= conserve_prim(w, γ)

    #--- record residuals ---#
    @. RES += (w - w_old)^2
    @. AVG += abs(w)

    #--- calculate M^{n+1} and tau^{n+1} ---#
    MH = maxwellian(u, v, prim)
    MB = MH .* K ./ (2.0 * prim[end])
    MH .+= SH
    MB .+= SB
    τ = vhs_collision_time(prim, μᵣ, ω)

    #--- update distribution function ---#
    for i in eachindex(u)
        h[i] =
            (h[i] + (fhL[i] - fhR[i] + fhD[i] - fhU[i]) / Δs + dt / τ * MH[i]) /
            (1.0 + dt / τ)
        b[i] =
            (b[i] + (fbL[i] - fbR[i] + fbD[i] - fbU[i]) / Δs + dt / τ * MB[i]) /
            (1.0 + dt / τ)
    end

end

"""
$(SIGNATURES)

2D2F2V @ triangle
"""
function step!(
    w::T1,
    prim::T1,
    h::T2,
    b::T2,
    fw1::T1,
    fh1::T2,
    fb1::T2,
    fw2::T1,
    fh2::T2,
    fb2::T2,
    fw3::T1,
    fh3::T2,
    fb3::T2,
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
) where {T1<:AV,T2<:AVOM,T3<:AVOM}

    #--- store W^n and calculate shakhov term ---#
    w_old = deepcopy(w)

    if collision == :shakhov
        q = heat_flux(h, b, prim, u, v, weights)

        MH_old = maxwellian(u, v, prim)
        MB_old = MH_old .* K ./ (2.0 * prim[end])
        SH, SB = shakhov(u, v, MH_old, MB_old, q, prim, Pr, K)
    else
        SH = zero(h)
        SB = zero(b)
    end

    #--- update W^{n+1} ---#
    @. w -= (fw1 * dirc[1] + fw2 * dirc[2] + fw3 * dirc[3]) / Δs
    prim .= conserve_prim(w, γ)

    #--- record residuals ---#
    @. RES += (w - w_old)^2
    @. AVG += abs(w)

    #--- calculate M^{n+1} and tau^{n+1} ---#
    MH = maxwellian(u, v, prim)
    MB = MH .* K ./ (2.0 * prim[end])
    MH .+= SH
    MB .+= SB
    τ = vhs_collision_time(prim, μᵣ, ω)

    #--- update distribution function ---#
    for i in eachindex(u)
        h[i] =
            (
                h[i] - (fh1[i] * dirc[1] + fh2[i] * dirc[2] + fh3[i] * dirc[3]) / Δs +
                dt / τ * MH[i]
            ) / (1.0 + dt / τ)
        b[i] =
            (
                b[i] - (fb1[i] * dirc[1] + fb2[i] * dirc[2] + fb3[i] * dirc[3]) / Δs +
                dt / τ * MB[i]
            ) / (1.0 + dt / τ)
    end

    return nothing

end
