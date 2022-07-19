"""
$(SIGNATURES)

Update flow variables with finite volume method

1D
"""
step!(
    KS::AbstractSolverSet,
    cell::AbstractControlVolume,
    faceL::T,
    faceR::T,
    p,
    coll = :bgk,
) where {T<:AbstractInterface} = step!(KS, KS.vs, KS.gas, cell, faceL, faceR, p, coll)

"""
$(SIGNATURES)

2D
"""
step!(
    KS::AbstractSolverSet,
    cell::AbstractControlVolume,
    faceL::T,
    faceR::T,
    faceD::T,
    faceU::T,
    p,
    coll = :bgk,
) where {T<:AbstractInterface} =
    step!(KS, KS.vs, KS.gas, cell, faceL, faceR, faceD, faceU, p, coll)

"""
$(SIGNATURES)

1D scalar
"""
function step!(
    KS,
    vs,
    gas::Scalar,
    cell::TC,
    faceL,
    faceR,
    p,
    coll = :bgk,
) where {TC<:Union{ControlVolume,ControlVolume1D}}
    dt, dx, RES, AVG = p
    step!(cell.w, cell.prim, faceL.fw, faceR.fw, gas.a, dx, RES, AVG)
end

"""
$(SIGNATURES)

1D0F
"""
function step!(
    KS,
    vs,
    gas::Gas,
    cell::TC,
    faceL,
    faceR,
    p,
    coll = :bgk,
) where {TC<:Union{ControlVolume,ControlVolume1D}}
    dt, dx, RES, AVG = p
    step!(cell.w, cell.prim, faceL.fw, faceR.fw, gas.γ, dx, RES, AVG)
end

"""
$(SIGNATURES)

1D0F mixture
"""
function step!(
    KS,
    vs,
    gas::Mixture,
    cell::TC,
    faceL,
    faceR,
    p,
    coll = :bgk,
) where {TC<:Union{ControlVolume,ControlVolume1D}}
    dt, dx, RES, AVG = p
    step!(
        cell.w,
        cell.prim,
        faceL.fw,
        faceR.fw,
        gas.γ,
        gas.mi,
        gas.ni,
        gas.me,
        gas.ne,
        gas.Kn[1],
        dx,
        dt,
        RES,
        AVG,
    )
end

"""
$(SIGNATURES)

1D1F1V
"""
function step!(
    KS,
    vs::AbstractVelocitySpace1D,
    gas::Gas,
    cell::TC,
    faceL,
    faceR,
    p,
    coll = :bgk,
) where {TC<:Union{ControlVolume1F,ControlVolume1D1F}}
    dt, dx, RES, AVG = p
    step!(
        cell.w,
        cell.prim,
        cell.f,
        faceL.fw,
        faceL.ff,
        faceR.fw,
        faceR.ff,
        KS.vs.u,
        KS.vs.weights,
        gas.γ,
        gas.μᵣ,
        gas.ω,
        gas.Pr,
        dx,
        dt,
        RES,
        AVG,
        coll,
    )
end

"""
$(SIGNATURES)

1D1F3V
"""
function step!(
    KS,
    vs::AbstractVelocitySpace3D,
    gas::Gas,
    cell::TC,
    faceL,
    faceR,
    p,
    coll = :bgk,
) where {TC<:Union{ControlVolume1F,ControlVolume1D1F}}
    dt, dx, RES, AVG = p
    if coll == :fsm
        step!(
            cell.w,
            cell.prim,
            cell.f,
            faceL.fw,
            faceL.ff,
            faceR.fw,
            faceR.ff,
            gas.γ,
            gas.fsm.Kn,
            gas.fsm.nm,
            gas.fsm.ϕ,
            gas.fsm.ψ,
            gas.fsm.χ,
            dx,
            dt,
            RES,
            AVG,
            coll,
        )
    else
        step!(
            cell.w,
            cell.prim,
            cell.f,
            faceL.fw,
            faceL.ff,
            faceR.fw,
            faceR.ff,
            vs.u,
            vs.v,
            vs.w,
            vs.weights,
            gas.γ,
            gas.μᵣ,
            gas.ω,
            gas.Pr,
            dx,
            dt,
            RES,
            AVG,
            coll,
        )
    end
end

"""
$(SIGNATURES)

1D2F1V
"""
function step!(
    KS,
    vs::AbstractVelocitySpace1D,
    gas::Gas,
    cell::TC,
    faceL,
    faceR,
    p,
    coll = :bgk,
) where {TC<:Union{ControlVolume2F,ControlVolume1D2F}}
    dt, dx, RES, AVG = p
    step!(
        cell.w,
        cell.prim,
        cell.h,
        cell.b,
        faceL.fw,
        faceL.fh,
        faceL.fb,
        faceR.fw,
        faceR.fh,
        faceR.fb,
        KS.vs.u,
        KS.vs.weights,
        gas.K,
        gas.γ,
        gas.μᵣ,
        gas.ω,
        gas.Pr,
        dx,
        dt,
        RES,
        AVG,
        coll,
    )
end

"""
$(SIGNATURES)

1D2F1V mixture
"""
function step!(
    KS,
    vs::AbstractVelocitySpace1D,
    gas::Mixture,
    cell::TC,
    faceL,
    faceR,
    p,
    coll = :bgk,
) where {TC<:Union{ControlVolume2F,ControlVolume1D2F}}
    dt, dx, RES, AVG = p
    step!(
        cell.w,
        cell.prim,
        cell.h,
        cell.b,
        faceL.fw,
        faceL.fh,
        faceL.fb,
        faceR.fw,
        faceR.fh,
        faceR.fb,
        KS.vs.u,
        KS.vs.weights,
        gas.K,
        gas.γ,
        gas.mi,
        gas.ni,
        gas.me,
        gas.ne,
        gas.Kn[1],
        gas.Pr,
        dx,
        dt,
        RES,
        AVG,
        coll,
    )
end

"""
$(SIGNATURES)

2D0F
"""
function step!(
    KS,
    vs,
    gas::Gas,
    cell::TC,
    faceL,
    faceR,
    faceD,
    faceU,
    p,
    coll = :bgk,
) where {TC<:Union{ControlVolume,ControlVolume2D}}
    dt, Δs, RES, AVG = p
    step!(
        cell.w,
        cell.prim,
        faceL.fw,
        faceR.fw,
        faceD.fw,
        faceU.fw,
        gas.γ,
        Δs,
        RES,
        AVG,
        coll,
    )
end

"""
$(SIGNATURES)

2D1F2V
"""
function step!(
    KS,
    vs::AbstractVelocitySpace2D,
    gas::Gas,
    cell::TC,
    faceL,
    faceR,
    faceD,
    faceU,
    p,
    coll = :bgk,
) where {TC<:Union{ControlVolume1F,ControlVolume2D1F}}
    dt, Δs, RES, AVG = p
    step!(
        cell.w,
        cell.prim,
        cell.f,
        faceL.fw,
        faceL.ff,
        faceR.fw,
        faceR.ff,
        faceD.fw,
        faceD.ff,
        faceU.fw,
        faceU.ff,
        vs.u,
        vs.v,
        vs.weights,
        gas.γ,
        gas.μᵣ,
        gas.ω,
        gas.Pr,
        Δs,
        dt,
        RES,
        AVG,
        coll,
    )
end

"""
$(SIGNATURES)

2D2F2V
"""
function step!(
    KS,
    vs::AbstractVelocitySpace2D,
    gas::Gas,
    cell::TC,
    faceL,
    faceR,
    faceD,
    faceU,
    p,
    coll = :bgk,
) where {TC<:Union{ControlVolume2F,ControlVolume2D2F}}
    dt, Δs, RES, AVG = p
    step!(
        cell.w,
        cell.prim,
        cell.h,
        cell.b,
        faceL.fw,
        faceL.fh,
        faceL.fb,
        faceR.fw,
        faceR.fh,
        faceR.fb,
        faceD.fw,
        faceD.fh,
        faceD.fb,
        faceU.fw,
        faceU.fh,
        faceU.fb,
        vs.u,
        vs.v,
        vs.weights,
        gas.K,
        gas.γ,
        gas.μᵣ,
        gas.ω,
        gas.Pr,
        Δs,
        dt,
        RES,
        AVG,
        coll,
    )
end

# ------------------------------------------------------------
# Low-level backends
# ------------------------------------------------------------

"""
$(SIGNATURES)

Scalar
"""
function step!(w::X, prim::AV{X}, fwL::X, fwR::X, a, dx, RES, AVG) where {X<:FN}
    #--- store W^n and calculate H^n,\tau^n ---#
    w_old = deepcopy(w)

    #--- update W^{n+1} ---#
    w += (fwL - fwR) / dx
    prim .= conserve_prim(w, a)

    #--- record residuals ---#
    RES += (w - w_old)^2
    AVG += abs(w)

    return w, RES, AVG
end

"""
$(SIGNATURES)

1D0F
"""
function step!(
    w::Y,
    prim::Y,
    fwL::X,
    fwR::X,
    γ,
    dx,
    RES,
    AVG,
) where {X<:AA{<:FN,1},Y<:AA{<:FN,1}} # 1D0F

    #--- store W^n and calculate H^n,\tau^n ---#
    w_old = deepcopy(w)

    #--- update W^{n+1} ---#
    @. w += (fwL - fwR) / dx
    prim .= conserve_prim(w, γ)

    #--- record residuals ---#
    @. RES += (w - w_old)^2
    @. AVG += abs(w)

    return nothing

end

"""
$(SIGNATURES)

1D0F Mixture
"""
function step!(
    w::T2,
    prim::T2,
    fwL::T1,
    fwR::T1,
    γ,
    mi,
    ni,
    me,
    ne,
    Kn,
    dx,
    dt,
    RES,
    AVG,
) where {T1<:AA{<:FN,2},T2<:AA{<:FN,2}}

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
    # explicit
    tau = aap_hs_collision_time(prim, mi, ni, me, ne, Kn)
    mprim = aap_hs_prim(prim, tau, mi, ni, me, ne, Kn)
    mw = mixture_prim_conserve(mprim, γ)
    for k in axes(w, 2)
        @. w[:, k] += (mw[:, k] - w_old[:, k]) * dt / tau[k]
    end
    prim .= mixture_conserve_prim(w, γ)

    #--- record residuals ---#
    @. RES += (w_old - w)^2
    @. AVG += abs(w)

end

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
    collision = :bgk::Symbol,
) where {T1<:AA{<:FN,1},T2<:AA{<:FN,1},T3<:AA{<:FN,1},T4<:AA{<:FN,1},T5<:AA{<:FN,1}}

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
    collision = :bgk::Symbol,
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
    collision = :fsm::Symbol,
) where {T1<:AA{<:FN,1},T2<:AA{<:FN,3},T3<:AA{<:FN,1},T4<:AA{<:FN,3}}

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
    collision = :bgk::Symbol,
) where {T1<:AA{<:FN,1},T2<:AA{<:FN,1},T3<:AA{<:FN,1},T4<:AA{<:FN,1},T5<:AA{<:FN,1}}

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
    collision = :bgk::Symbol,
) where {T1<:AA{<:FN,2},T2<:AA{<:FN,2},T3<:AA{<:FN,2},T4<:AA{<:FN,2},T5<:AA{<:FN,2}}

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

1D3F1V Rykov
"""
function step!(
    w::T3,
    prim::T3,
    h::T4,
    b::T4,
    r::T4,
    fwL::T1,
    fhL::T2,
    fbL::T2,
    frL::T2,
    fwR::T1,
    fhR::T2,
    fbR::T2,
    frR::T2,
    u::T5,
    weights::T5,
    K,
    Kr,
    μᵣ,
    ω,
    Pr,
    T₀,
    Z₀,
    σ,
    ω0,
    ω1,
    dx,
    dt,
    RES,
    AVG,
    collision = :rykov::Symbol,
) where {T1<:AA{<:FN,1},T2<:AA{<:FN,1},T3<:AA{<:Real,1},T4<:AA{<:FN,1},T5<:AA{<:FN,1}}

    #--- store W^n and calculate shakhov term ---#
    w_old = deepcopy(w)

    if collision == :rykov
        q = heat_flux(h, b, r, prim, u, weights)
    else
        q = zeros(2)
    end

    #--- update W^{n+1} ---#
    @. w += (fwL - fwR) / dx

    MHT = similar(h)
    MBT = similar(b)
    MRT = similar(r)
    MHR = similar(h)
    MBR = similar(b)
    MRR = similar(r)
    maxwellian!(MHT, MBT, MRT, MHR, MBR, MRR, u, prim, K, Kr)
    τ_old = vhs_collision_time(prim[1:end-1], μᵣ, ω)
    Zr = rykov_zr(1.0 / prim[4], T₀, Z₀)
    Er0_old = 0.5 * sum(@. weights * ((1.0 / Zr) * MRR + (1.0 - 1.0 / Zr) * MRT))

    w[4] += dt * (Er0_old - w_old[4]) / τ_old
    prim .= conserve_prim(w, K, Kr)

    #--- record residuals ---#
    @. RES += (w - w_old)^2
    @. AVG += abs(w)

    #--- calculate M^{n+1} and tau^{n+1} ---#
    maxwellian!(MHT, MBT, MRT, MHR, MBR, MRR, u, prim, K, Kr)

    SHT = similar(h)
    SBT = similar(b)
    SRT = similar(r)
    SHR = similar(h)
    SBR = similar(b)
    SRR = similar(r)
    rykov!(
        SHT,
        SBT,
        SRT,
        SHR,
        SBR,
        SRR,
        u,
        MHT,
        MBT,
        MRT,
        MHR,
        MBR,
        MRR,
        q,
        prim,
        Pr,
        K,
        σ,
        ω0,
        ω1,
    )

    MH = (1.0 - 1.0 / Zr) * (MHT + SHT) + 1.0 / Zr * (MHR + SHR)
    MB = (1.0 - 1.0 / Zr) * (MBT + SBT) + 1.0 / Zr * (MBR + SBR)
    MR = (1.0 - 1.0 / Zr) * (MRT + SRT) + 1.0 / Zr * (MRR + SRR)

    τ = vhs_collision_time(prim[1:end-1], μᵣ, ω)

    #--- update distribution function ---#
    for i in eachindex(u)
        h[i] = (h[i] + (fhL[i] - fhR[i]) / dx + dt / τ * MH[i]) / (1.0 + dt / τ)
        b[i] = (b[i] + (fbL[i] - fbR[i]) / dx + dt / τ * MB[i]) / (1.0 + dt / τ)
        r[i] = (r[i] + (frL[i] - frR[i]) / dx + dt / τ * MR[i]) / (1.0 + dt / τ)
    end

end

"""
$(SIGNATURES)

1D4F1V
"""
function step!(
    KS::T,
    cell::ControlVolume1D4F,
    faceL::Interface1D4F,
    faceR::Interface1D4F,
    dx,
    dt,
    RES,
    AVG,
    collision = :bgk::Symbol,
    isMHD = true::Bool,
) where {T<:AbstractSolverSet}

    #--- update conservative flow variables: step 1 ---#
    # w^n
    w_old = deepcopy(cell.w)
    prim_old = deepcopy(cell.prim)

    # flux -> w^{n+1}
    @. cell.w += (faceL.fw - faceR.fw) / dx
    cell.prim .= mixture_conserve_prim(cell.w, KS.gas.γ)

    # temperature protection
    if cell.prim[5, 1] < 0
        @warn ("ion temperature update is negative")
        cell.w .= w_old
        cell.prim .= prim_old
    elseif cell.prim[5, 2] < 0
        @warn ("electron temperature update is negative")
        cell.w .= w_old
        cell.prim .= prim_old
    end

    # source -> w^{n+1}
    if isMHD == false
        #=
        # DifferentialEquations.jl
        tau = get_tau(cell.prim, KS.gas.mi, KS.gas.ni, KS.gas.me, KS.gas.ne, KS.gas.Kn[1])
        for j in axes(wRan, 2)
        prob = ODEProblem( mixture_source,
            vcat(cell.w[1:5,j,1], cell.w[1:5,j,2]),
            dt,
            (tau[1], tau[2], KS.gas.mi, KS.gas.ni, KS.gas.me, KS.gas.ne, KS.gas.Kn[1], KS.gas.γ) )
        sol = solve(prob, Rosenbrock23())

        cell.w[1:5,j,1] .= sol[end][1:5]
        cell.w[1:5,j,2] .= sol[end][6:10]
        for k=1:2
        cell.prim[:,j,k] .= conserve_prim(cell.w[:,j,k], KS.gas.γ)
        end
        end
        =#

        # explicit
        tau = aap_hs_collision_time(
            cell.prim,
            KS.gas.mi,
            KS.gas.ni,
            KS.gas.me,
            KS.gas.ne,
            KS.gas.Kn[1],
        )
        mprim = aap_hs_prim(
            cell.prim,
            tau,
            KS.gas.mi,
            KS.gas.ni,
            KS.gas.me,
            KS.gas.ne,
            KS.gas.Kn[1],
        )
        mw = mixture_prim_conserve(mprim, KS.gas.γ)
        for k = 1:2
            @. cell.w[:, k] += (mw[:, k] - w_old[:, k]) * dt / tau[k]
        end
        cell.prim .= mixture_conserve_prim(cell.w, KS.gas.γ)
    end

    #--- update electromagnetic variables ---#
    # flux -> E^{n+1} & B^{n+1}
    cell.E[1] -= dt * (faceL.femR[1] + faceR.femL[1]) / dx
    cell.E[2] -= dt * (faceL.femR[2] + faceR.femL[2]) / dx
    cell.E[3] -= dt * (faceL.femR[3] + faceR.femL[3]) / dx
    cell.B[1] -= dt * (faceL.femR[4] + faceR.femL[4]) / dx
    cell.B[2] -= dt * (faceL.femR[5] + faceR.femL[5]) / dx
    cell.B[3] -= dt * (faceL.femR[6] + faceR.femL[6]) / dx
    cell.ϕ -= dt * (faceL.femR[7] + faceR.femL[7]) / dx
    cell.ψ -= dt * (faceL.femR[8] + faceR.femL[8]) / dx

    for i = 1:3
        if 1 ∈ vcat(isnan.(cell.E), isnan.(cell.B))
            @warn "NaN electromagnetic update"
        end
    end

    # source -> ϕ
    #@. cell.ϕ += dt * (cell.w[1,:,1] / KS.gas.mi - cell.w[1,:,2] / KS.gas.me) / (KS.gas.lD^2 * KS.gas.rL)

    # source -> U^{n+1}, E^{n+1} and B^{n+1}
    mr = KS.gas.mi / KS.gas.me
    A, b = em_coefficients(cell.prim, cell.E, cell.B, mr, KS.gas.lD, KS.gas.rL, dt)
    x = A \ b

    #--- calculate lorenz force ---#
    cell.lorenz[1, 1] =
        0.5 * (
            x[1] + cell.E[1] + (cell.prim[3, 1] + x[5]) * cell.B[3] -
            (cell.prim[4, 1] + x[6]) * cell.B[2]
        ) / KS.gas.rL
    cell.lorenz[2, 1] =
        0.5 * (
            x[2] + cell.E[2] + (cell.prim[4, 1] + x[6]) * cell.B[1] -
            (cell.prim[2, 1] + x[4]) * cell.B[3]
        ) / KS.gas.rL
    cell.lorenz[3, 1] =
        0.5 * (
            x[3] + cell.E[3] + (cell.prim[2, 1] + x[4]) * cell.B[2] -
            (cell.prim[3, 1] + x[5]) * cell.B[1]
        ) / KS.gas.rL
    cell.lorenz[1, 2] =
        -0.5 *
        (
            x[1] + cell.E[1] + (cell.prim[3, 2] + x[8]) * cell.B[3] -
            (cell.prim[4, 2] + x[9]) * cell.B[2]
        ) *
        mr / KS.gas.rL
    cell.lorenz[2, 2] =
        -0.5 *
        (
            x[2] + cell.E[2] + (cell.prim[4, 2] + x[9]) * cell.B[1] -
            (cell.prim[2, 2] + x[7]) * cell.B[3]
        ) *
        mr / KS.gas.rL
    cell.lorenz[3, 2] =
        -0.5 *
        (
            x[3] + cell.E[3] + (cell.prim[2, 2] + x[7]) * cell.B[2] -
            (cell.prim[3, 2] + x[8]) * cell.B[1]
        ) *
        mr / KS.gas.rL

    cell.E[1] = x[1]
    cell.E[2] = x[2]
    cell.E[3] = x[3]

    #--- update conservative flow variables: step 2 ---#
    cell.prim[2, 1] = x[4]
    cell.prim[3, 1] = x[5]
    cell.prim[4, 1] = x[6]
    cell.prim[2, 2] = x[7]
    cell.prim[3, 2] = x[8]
    cell.prim[4, 2] = x[9]

    cell.w .= mixture_prim_conserve(cell.prim, KS.gas.γ)

    #--- update particle distribution function ---#
    # flux -> f^{n+1}
    @. cell.h0 += (faceL.fh0 - faceR.fh0) / dx
    @. cell.h1 += (faceL.fh1 - faceR.fh1) / dx
    @. cell.h2 += (faceL.fh2 - faceR.fh2) / dx
    @. cell.h3 += (faceL.fh3 - faceR.fh3) / dx

    # force -> f^{n+1} : step 1
    for j in axes(cell.h0, 2)
        _h0 = @view cell.h0[:, j]
        _h1 = @view cell.h1[:, j]
        _h2 = @view cell.h2[:, j]
        _h3 = @view cell.h3[:, j]

        shift_pdf!(_h0, cell.lorenz[1, j], KS.vs.du[1, j], dt)
        shift_pdf!(_h1, cell.lorenz[1, j], KS.vs.du[1, j], dt)
        shift_pdf!(_h2, cell.lorenz[1, j], KS.vs.du[1, j], dt)
        shift_pdf!(_h3, cell.lorenz[1, j], KS.vs.du[1, j], dt)
    end

    # force -> f^{n+1} : step 2
    for k in axes(cell.h1, 3)
        @. cell.h3[:, k] +=
            2.0 * dt * cell.lorenz[2, k] * cell.h1[:, k] +
            (dt * cell.lorenz[2, k])^2 * cell.h0[:, k] +
            2.0 * dt * cell.lorenz[3, k] * cell.h2[:, k] +
            (dt * cell.lorenz[3, k])^2 * cell.h0[:, k]
        @. cell.h2[:, k] += dt * cell.lorenz[3, k] * cell.h0[:, k]
        @. cell.h1[:, k] += dt * cell.lorenz[2, k] * cell.h0[:, k]
    end

    # source -> f^{n+1}
    tau = aap_hs_collision_time(
        cell.prim,
        KS.gas.mi,
        KS.gas.ni,
        KS.gas.me,
        KS.gas.ne,
        KS.gas.Kn[1],
    )

    # interspecies interaction
    if isMHD == true
        prim = deepcopy(cell.prim)
    else
        prim = aap_hs_prim(
            cell.prim,
            tau,
            KS.gas.mi,
            KS.gas.ni,
            KS.gas.me,
            KS.gas.ne,
            KS.gas.Kn[1],
        )
    end
    g = mixture_maxwellian(KS.vs.u, prim)

    # BGK term
    Mu, Mv, Mw, MuL, MuR = mixture_gauss_moments(prim, KS.gas.K)
    for k in axes(cell.h0, 2)
        @. cell.h0[:, k] = (cell.h0[:, k] + dt / tau[k] * g[:, k]) / (1.0 + dt / tau[k])
        @. cell.h1[:, k] =
            (cell.h1[:, k] + dt / tau[k] * Mv[1, k] * g[:, k]) / (1.0 + dt / tau[k])
        @. cell.h2[:, k] =
            (cell.h2[:, k] + dt / tau[k] * Mw[1, k] * g[:, k]) / (1.0 + dt / tau[k])
        @. cell.h3[:, k] =
            (cell.h3[:, k] + dt / tau[k] * (Mv[2, k] + Mw[2, k]) * g[:, k]) /
            (1.0 + dt / tau[k])
    end

    #--- record residuals ---#
    @. RES += (w_old - cell.w)^2
    @. AVG += abs(cell.w)

end

"""
$(SIGNATURES)

1D3F2V
"""
function step!(
    KS::T,
    cell::ControlVolume1D3F,
    faceL::Interface1D3F,
    faceR::Interface1D3F,
    dx,
    dt,
    RES,
    AVG,
    collision = :bgk::Symbol,
    isMHD = true::Bool,
) where {T<:AbstractSolverSet}

    #--- update conservative flow variables: step 1 ---#
    # w^n
    w_old = deepcopy(cell.w)
    prim_old = deepcopy(cell.prim)

    # flux -> w^{n+1}
    @. cell.w += (faceL.fw - faceR.fw) / dx
    cell.prim .= mixture_conserve_prim(cell.w, KS.gas.γ)

    # temperature protection
    if cell.prim[end, 1] < 0
        @warn ("ion temperature update is negative")
        cell.w .= w_old
        cell.prim .= prim_old
    elseif cell.prim[end, 2] < 0
        @warn ("electron temperature update is negative")
        cell.w .= w_old
        cell.prim .= prim_old
    end

    # source -> w^{n+1}
    if isMHD == false
        #=
        # DifferentialEquations.jl
        tau = get_tau(cell.prim, KS.gas.mi, KS.gas.ni, KS.gas.me, KS.gas.ne, KS.gas.Kn[1])
        for j in axes(wRan, 2)
        prob = ODEProblem( mixture_source,
            vcat(cell.w[1:5,j,1], cell.w[1:5,j,2]),
            dt,
            (tau[1], tau[2], KS.gas.mi, KS.gas.ni, KS.gas.me, KS.gas.ne, KS.gas.Kn[1], KS.gas.γ) )
        sol = solve(prob, Rosenbrock23())

        cell.w[1:5,j,1] .= sol[end][1:5]
        cell.w[1:5,j,2] .= sol[end][6:10]
        for k=1:2
        cell.prim[:,j,k] .= conserve_prim(cell.w[:,j,k], KS.gas.γ)
        end
        end
        =#

        # explicit
        tau = aap_hs_collision_time(
            cell.prim,
            KS.gas.mi,
            KS.gas.ni,
            KS.gas.me,
            KS.gas.ne,
            KS.gas.Kn[1],
        )
        mprim = aap_hs_prim(
            cell.prim,
            tau,
            KS.gas.mi,
            KS.gas.ni,
            KS.gas.me,
            KS.gas.ne,
            KS.gas.Kn[1],
        )
        mw = mixture_prim_conserve(mprim, KS.gas.γ)
        for k in axes(cell.w, 2)
            @. cell.w[:, k] += (mw[:, k] - w_old[:, k]) * dt / tau[k]
        end
        cell.prim .= mixture_conserve_prim(cell.w, KS.gas.γ)
    end

    #--- update electromagnetic variables ---#
    # flux -> E^{n+1} & B^{n+1}
    cell.E[1] -= dt * (faceL.femR[1] + faceR.femL[1]) / dx
    cell.E[2] -= dt * (faceL.femR[2] + faceR.femL[2]) / dx
    cell.E[3] -= dt * (faceL.femR[3] + faceR.femL[3]) / dx
    cell.B[1] -= dt * (faceL.femR[4] + faceR.femL[4]) / dx
    cell.B[2] -= dt * (faceL.femR[5] + faceR.femL[5]) / dx
    cell.B[3] -= dt * (faceL.femR[6] + faceR.femL[6]) / dx
    cell.ϕ -= dt * (faceL.femR[7] + faceR.femL[7]) / dx
    cell.ψ -= dt * (faceL.femR[8] + faceR.femL[8]) / dx

    for i = 1:3
        if 1 ∈ vcat(isnan.(cell.E), isnan.(cell.B))
            @warn "electromagnetic update is NaN"
        end
    end

    # source -> ϕ
    #@. cell.ϕ += dt * (cell.w[1,:,1] / KS.gas.mi - cell.w[1,:,2] / KS.gas.me) / (KS.gas.lD^2 * KS.gas.rL)

    # source -> U^{n+1}, E^{n+1} and B^{n+1}
    mr = KS.gas.mi / KS.gas.me
    A, b = em_coefficients(cell.prim, cell.E, cell.B, mr, KS.gas.lD, KS.gas.rL, dt)
    x = A \ b

    #--- calculate lorenz force ---#
    cell.lorenz[1, 1] =
        0.5 * (
            x[1] + cell.E[1] + (cell.prim[3, 1] + x[5]) * cell.B[3] -
            (cell.prim[4, 1] + x[6]) * cell.B[2]
        ) / KS.gas.rL
    cell.lorenz[2, 1] =
        0.5 * (
            x[2] + cell.E[2] + (cell.prim[4, 1] + x[6]) * cell.B[1] -
            (cell.prim[2, 1] + x[4]) * cell.B[3]
        ) / KS.gas.rL
    cell.lorenz[3, 1] =
        0.5 * (
            x[3] + cell.E[3] + (cell.prim[2, 1] + x[4]) * cell.B[2] -
            (cell.prim[3, 1] + x[5]) * cell.B[1]
        ) / KS.gas.rL
    cell.lorenz[1, 2] =
        -0.5 *
        (
            x[1] + cell.E[1] + (cell.prim[3, 2] + x[8]) * cell.B[3] -
            (cell.prim[4, 2] + x[9]) * cell.B[2]
        ) *
        mr / KS.gas.rL
    cell.lorenz[2, 2] =
        -0.5 *
        (
            x[2] + cell.E[2] + (cell.prim[4, 2] + x[9]) * cell.B[1] -
            (cell.prim[2, 2] + x[7]) * cell.B[3]
        ) *
        mr / KS.gas.rL
    cell.lorenz[3, 2] =
        -0.5 *
        (
            x[3] + cell.E[3] + (cell.prim[2, 2] + x[7]) * cell.B[2] -
            (cell.prim[3, 2] + x[8]) * cell.B[1]
        ) *
        mr / KS.gas.rL

    cell.E[1] = x[1]
    cell.E[2] = x[2]
    cell.E[3] = x[3]

    #--- update conservative flow variables: step 2 ---#
    cell.prim[2, 1] = x[4]
    cell.prim[3, 1] = x[5]
    cell.prim[4, 1] = x[6]
    cell.prim[2, 2] = x[7]
    cell.prim[3, 2] = x[8]
    cell.prim[4, 2] = x[9]

    cell.w .= mixture_prim_conserve(cell.prim, KS.gas.γ)

    #--- update particle distribution function ---#
    # flux -> f^{n+1}
    @. cell.h0 += (faceL.fh0 - faceR.fh0) / dx
    @. cell.h1 += (faceL.fh1 - faceR.fh1) / dx
    @. cell.h2 += (faceL.fh2 - faceR.fh2) / dx

    # force -> f^{n+1} : step 1
    for j in axes(cell.h0, 3) # component
        for i in axes(cell.h0, 2) # v
            _h0 = @view cell.h0[:, i, j]
            _h1 = @view cell.h1[:, i, j]
            _h2 = @view cell.h2[:, i, j]

            shift_pdf!(_h0, cell.lorenz[1, j], KS.vs.du[1, i, j], dt)
            shift_pdf!(_h1, cell.lorenz[1, j], KS.vs.du[1, i, j], dt)
            shift_pdf!(_h2, cell.lorenz[1, j], KS.vs.du[1, i, j], dt)
        end
    end

    for j in axes(cell.h0, 3) # component
        for i in axes(cell.h0, 1) # u
            _h0 = @view cell.h0[i, :, j]
            _h1 = @view cell.h1[i, :, j]
            _h2 = @view cell.h2[i, :, j]

            shift_pdf!(_h0, cell.lorenz[2, j], KS.vs.dv[i, 1, j], dt)
            shift_pdf!(_h1, cell.lorenz[2, j], KS.vs.dv[i, 1, j], dt)
            shift_pdf!(_h2, cell.lorenz[2, j], KS.vs.dv[i, 1, j], dt)
        end
    end

    # force -> f^{n+1} : step 2
    for k in axes(cell.h1, 3)
        @. cell.h2[:, :, k] +=
            2.0 * dt * cell.lorenz[3, k] * cell.h1[:, :, k] +
            (dt * cell.lorenz[3, k])^2 * cell.h0[:, :, k]
        @. cell.h1[:, :, k] += dt * cell.lorenz[3, k] * cell.h0[:, :, k]
    end

    # source -> f^{n+1}
    tau = aap_hs_collision_time(
        cell.prim,
        KS.gas.mi,
        KS.gas.ni,
        KS.gas.me,
        KS.gas.ne,
        KS.gas.Kn[1],
    )

    # interspecies interaction
    if isMHD == true
        prim = deepcopy(cell.prim)
    else
        prim = aap_hs_prim(
            cell.prim,
            tau,
            KS.gas.mi,
            KS.gas.ni,
            KS.gas.me,
            KS.gas.ne,
            KS.gas.Kn[1],
        )
    end

    H0 = similar(KS.vs.u)
    H1 = similar(H0)
    H2 = similar(H0)
    for k in axes(H0, 3)
        H0[:, :, k] .= maxwellian(KS.vs.u[:, :, k], KS.vs.v[:, :, k], prim[:, k])
        @. H1[:, :, k] = H0[:, :, k] * prim[4, k]
        @. H2[:, :, k] = H0[:, :, k] * (prim[4, k]^2 + 1.0 / (2.0 * prim[5, k]))
    end

    # BGK term
    for k in axes(cell.h0, 3)
        @. cell.h0[:, :, k] =
            (cell.h0[:, :, k] + dt / tau[k] * H0[:, :, k]) / (1.0 + dt / tau[k])
        @. cell.h1[:, :, k] =
            (cell.h1[:, :, k] + dt / tau[k] * H1[:, :, k]) / (1.0 + dt / tau[k]) # NOTICE the h1 here is h2 in 1d4f case
        @. cell.h2[:, :, k] =
            (cell.h2[:, :, k] + dt / tau[k] * H2[:, :, k]) / (1.0 + dt / tau[k]) # NOTICE the h2 here is h3 in 1d4f case
    end

    #--- record residuals ---#
    @. RES += (w_old - cell.w)^2
    @. AVG += abs(cell.w)

end

"""
$(SIGNATURES)

2D0F @ triangle
"""
function step!(
    w::T1,
    prim::T1,
    fw1::T1,
    fw2::T1,
    fw3::T1,
    γ,
    Δs,
    dirc::T2,
    RES,
    AVG,
) where {T1<:AV{<:FN},T2<:AV{<:Real}}

    #--- store W^n and calculate shakhov term ---#
    w_old = deepcopy(w)

    #--- update W^{n+1} ---#
    @. w -= (fw1 * dirc[1] + fw2 * dirc[2] + fw3 * dirc[3]) / Δs
    prim .= conserve_prim(w, γ)

    #--- record residuals ---#
    @. RES += (w - w_old)^2
    @. AVG += abs(w)

end

"""
$(SIGNATURES)

2D0F @ quadrilateral
"""
function step!(
    w::T1,
    prim::T1,
    fwL::T1,
    fwR::T1,
    fwD::T1,
    fwU::T1,
    γ,
    Δs,
    RES,
    AVG,
    collision = :bgk,
) where {T1<:AA{<:FN,1}}

    #--- store W^n and calculate shakhov term ---#
    w_old = deepcopy(w)

    #--- update W^{n+1} ---#
    @. w += (fwL - fwR + fwD - fwU) / Δs
    prim .= conserve_prim(w, γ)

    #--- record residuals ---#
    @. RES += (w - w_old)^2
    @. AVG += abs(w)

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
    collision = :bgk,
) where {T1<:AA{<:FN,1},T2<:AA{<:FN,2},T3<:AA{<:FN,2}}

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
    for j in axes(v, 2), i in axes(u, 1)
        h[i, j] =
            (
                h[i, j] +
                (fhL[i, j] - fhR[i, j] + fhD[i, j] - fhU[i, j]) / Δs +
                dt / τ * MH[i, j]
            ) / (1.0 + dt / τ)
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
    dirc::T4,
    dt,
    RES,
    AVG,
    collision = :bgk,
) where {T1<:AV{<:FN},T2<:AA{<:FN,2},T3<:AA{<:FN,2},T4<:AV{<:Real}}

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
    for j in axes(v, 2), i in axes(u, 1)
        f[i, j] =
            (
                f[i, j] -
                (ff1[i, j] * dirc[1] + ff2[i, j] * dirc[2] + ff3[i, j] * dirc[3]) / Δs +
                dt / τ * M[i, j]
            ) / (1.0 + dt / τ)
    end

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
    collision = :bgk,
) where {T1<:AA{<:FN,1},T2<:AA{<:FN,2},T3<:AA{<:FN,2}}

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
    for j in axes(v, 2), i in axes(u, 1)
        h[i, j] =
            (
                h[i, j] +
                (fhL[i, j] - fhR[i, j] + fhD[i, j] - fhU[i, j]) / Δs +
                dt / τ * MH[i, j]
            ) / (1.0 + dt / τ)
        b[i, j] =
            (
                b[i, j] +
                (fbL[i, j] - fbR[i, j] + fbD[i, j] - fbU[i, j]) / Δs +
                dt / τ * MB[i, j]
            ) / (1.0 + dt / τ)
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
    dirc::T4,
    dt,
    RES,
    AVG,
    collision = :bgk,
) where {T1<:AV{<:FN},T2<:AA{<:FN,2},T3<:AA{<:FN,2},T4<:AV{<:Real}}

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
    for j in axes(v, 2), i in axes(u, 1)
        h[i, j] =
            (
                h[i, j] -
                (fh1[i, j] * dirc[1] + fh2[i, j] * dirc[2] + fh3[i, j] * dirc[3]) / Δs +
                dt / τ * MH[i, j]
            ) / (1.0 + dt / τ)
        b[i, j] =
            (
                b[i, j] -
                (fb1[i, j] * dirc[1] + fb2[i, j] * dirc[2] + fb3[i, j] * dirc[3]) / Δs +
                dt / τ * MB[i, j]
            ) / (1.0 + dt / τ)
    end

end

"""
$(SIGNATURES)
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
    collision = :fsm::Symbol,
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

end
