# ============================================================
# Gas Kinetic Flux Functions
# ============================================================

"""
$(SIGNATURES)

Gas kinetic flux
"""
flux_gks!(KS::AbstractSolverSet, face, ctrL, ctrR, args...) =
    flux_gks!(face, ctrL, ctrR, KS.gas, args...)

"""
$(SIGNATURES)

Scalar

`flux_gks` is called since there is no in-place operation for scalar
"""
function flux_gks!(face::Interface, ctrL, ctrR, gas::Scalar, p, dt = 1.0)

    dxL, dxR = p[1:2]

    face.fw = flux_gks(
        ctrL.w + dxL * ctrL.sw,
        ctrR.w - dxR * ctrR.sw,
        gas.μᵣ,
        dt,
        dxL,
        dxR,
        gas.a,
        ctrL.sw,
        ctrR.sw,
    )

    return nothing

end

"""
$(SIGNATURES)

Gas
"""
function flux_gks!(face::Interface, ctrL, ctrR, gas::Gas, p, dt = 1.0)

    dxL, dxR = p[1:2]

    if size(ctrL.w, 1) == 3
        flux_gks!(
            face.fw,
            ctrL.w .+ dxL .* ctrL.sw,
            ctrR.w .- dxR .* ctrR.sw,
            gas.γ,
            gas.K,
            gas.μᵣ,
            gas.ω,
            dt,
            dxL,
            dxR,
            ctrL.sw,
            ctrR.sw,
        )
    else
    end

    return nothing

end

"""
$(SIGNATURES)

Mixture
"""
function flux_gks!(face::Interface, ctrL, ctrR, gas::Mixture, p, dt = 1.0)

    dxL, dxR = p[1:2]

    if size(ctrL.w, 1) == 3
        flux_gks!(
            face.fw,
            ctrL.w .+ dxL .* ctrL.sw,
            ctrR.w .- dxR .* ctrR.sw,
            gas.γ,
            gas.K,
            gas.mi,
            gas.ni,
            gas.me,
            gas.ne,
            gas.Kn[1],
            dt,
            dxL,
            dxR,
            ctrL.sw,
            ctrR.sw,
        )
    else
    end

    return nothing

end

# ------------------------------------------------------------
# Low-level backends
# ------------------------------------------------------------

"""
$(SIGNATURES)

Gas kinetic scalar flux
"""
function flux_gks(u, μ, dt, a = 0, su = 0.0)
    prim = ifelse(a == 0, conserve_prim(u), conserve_prim(u, a))

    Mu = gauss_moments(prim)[1]
    tau = 2.0 * μ

    fa = pdf_slope(u, su)
    Δ = -prim[1] * moments_conserve_slope(fa, Mu, 1)
    faT = pdf_slope(u, Δ)

    # time-integration constants
    Mt = zeros(5)
    Mt[4] = dt
    Mt[5] = -tau * dt * exp(-dt / tau) + tau * Mt[4]

    # flux related to upwind distribution
    Muv = moments_conserve(Mu, 1)
    Mau = moments_conserve_slope(fa, Mu, 2)
    MauT = moments_conserve_slope(faT, Mu, 1)

    fw =
        Mt[4] * prim[1] * Muv - (Mt[5] + tau * Mt[4]) * prim[1] * Mau -
        tau * Mt[4] * prim[1] * MauT

    return fw / dt
end

"""
$(SIGNATURES)
"""
function flux_gks(
    uL::T,
    uR::T,
    μ,
    dt,
    dxL,
    dxR,
    a = 0,
    suL = 0.0,
    suR = 0.0,
) where {T<:Real}

    primL = ifelse(a == 0, conserve_prim(uL), conserve_prim(uL, a))
    primR = ifelse(a == 0, conserve_prim(uR), conserve_prim(uR, a))

    Mu1, MuL1, MuR1 = gauss_moments(primL)
    Mu2, MuL2, MuR2 = gauss_moments(primR)

    u = primL[1] * moments_conserve(MuL1, 0) + primR[1] * moments_conserve(MuR2, 0)
    prim = ifelse(a == 0, conserve_prim(u), conserve_prim(u, a))
    tau = 2.0 * μ + 1.0 * abs(uL - uR) / (abs(uL) + abs(uR)) * dt

    faL = pdf_slope(uL, suL)
    Δ = -primL[1] * moments_conserve_slope(faL, Mu1, 1)
    faTL = pdf_slope(uL, Δ)

    faR = pdf_slope(uR, suR)
    Δ = -primR[1] * moments_conserve_slope(faR, Mu2, 1)
    faTR = pdf_slope(uR, Δ)

    Mu, MuL, MuR = gauss_moments(prim)
    sw0L = (u - (uL - suL * dxL)) / dxL
    sw0R = ((uR + suR * dxR) - u) / dxR
    gaL = pdf_slope(u, sw0L)
    gaR = pdf_slope(u, sw0R)
    Δ =
        -prim[1] *
        (moments_conserve_slope(gaL, MuL, 1) + moments_conserve_slope(gaR, MuR, 1))
    #sw = (uR - uL) / (dxL + dxR)
    #ga = pdf_slope(u, sw)
    #Δ = -prim[1] .* moments_conserve_slope(ga, Mu, 1)
    gaT = pdf_slope(u, Δ)

    # time-integration constants
    Mt = zeros(5)
    Mt[4] = tau * (1.0 - exp(-dt / tau))
    Mt[5] = -tau * dt * exp(-dt / tau) + tau * Mt[4]
    Mt[1] = dt - Mt[4]
    Mt[2] = -tau * Mt[1] + Mt[5]
    Mt[3] = 0.5 * dt^2 - tau * Mt[1]

    # flux related to central distribution
    Muv = moments_conserve(Mu, 1)
    MauL = moments_conserve_slope(gaL, MuL, 2)
    MauR = moments_conserve_slope(gaR, MuR, 2)
    #Mau = moments_conserve_slope(ga, MuR, 2)
    MauT = moments_conserve_slope(gaT, Mu, 1)

    fw = Mt[1] * prim[1] * Muv + Mt[2] * prim[1] * (MauL + MauR) + Mt[3] * prim[1] * MauT
    #fw = Mt[1] * prim[1] * Muv + Mt[2] * prim[1] * Mau + Mt[3] * prim[1] * MauT

    # flux related to upwind distribution
    MuvL = moments_conserve(MuL1, 1)
    MauL = moments_conserve_slope(faL, MuL1, 2)
    MauLT = moments_conserve_slope(faTL, MuL1, 1)

    MuvR = moments_conserve(MuR2, 1)
    MauR = moments_conserve_slope(faR, MuR2, 2)
    MauRT = moments_conserve_slope(faTR, MuR2, 1)

    fw +=
        Mt[4] * primL[1] * MuvL - (Mt[5] + tau * Mt[4]) * primL[1] * MauL -
        tau * Mt[4] * primL[1] * MauLT + Mt[4] * primR[1] * MuvR -
        (Mt[5] + tau * Mt[4]) * primR[1] * MauR - tau * Mt[4] * primR[1] * MauRT
    # fw += Mt[4] * primL[1] * MuvL + Mt[4] * primR[1] * MuvR

    return fw

end


"""
$(SIGNATURES)

Gas kinetic Navier-Stokes flux

Continuous case
"""
function flux_gks!(fw::AV, w::Y, inK, γ, μᵣ, ω, sw = zero(w)::Y) where {Y<:AV}

    prim = conserve_prim(w, γ)
    mus = gauss_moments(prim, inK)
    tau = vhs_collision_time(prim, μᵣ, ω)

    ax = pdf_slope(prim, sw, inK)
    if length(prim) == 3
        ∂ft = -prim[1] .* moments_conserve_slope(ax, mus[1], mus[2], 1)
        at = pdf_slope(prim, ∂ft, inK)

        Muv = moments_conserve(mus[1], mus[2], 1, 0)
        Mau = moments_conserve_slope(ax, mus[1], mus[2], 2)
        Mtu = moments_conserve_slope(at, mus[1], mus[2], 1)
    elseif length(prim) == 4
        ∂ft = -prim[1] .* moments_conserve_slope(ax, mus[1], mus[2], mus[3], 1, 0)
        at = pdf_slope(prim, ∂ft, inK)

        Muv = moments_conserve(mus[1], mus[2], mus[3], 1, 0, 0)
        Mau = moments_conserve_slope(ax, mus[1], mus[2], mus[3], 2, 0)
        Mtu = moments_conserve_slope(at, mus[1], mus[2], mus[3], 1, 0)
    elseif length(prim) == 5
        ∂ft = -prim[1] .* moments_conserve_slope(ax, mus[1], mus[2], mus[3], 1, 0, 0)
        at = pdf_slope(prim, ∂ft, inK)

        Muv = moments_conserve(mus[1], mus[2], mus[3], 1, 0, 0)
        Mau = moments_conserve_slope(ax, mus[1], mus[2], mus[3], 2, 0, 0)
        Mtu = moments_conserve_slope(at, mus[1], mus[2], mus[3], 1, 0, 0)
    end

    @. fw = prim[1] * (Muv - tau * Mau - tau * Mtu)

    return nothing

end

"""
$(SIGNATURES)

1D
"""
function flux_gks!(
    fw::AV,
    wL::Y,
    wR::Y,
    inK,
    γ,
    μᵣ,
    ω,
    dt,
    dxL,
    dxR,
    swL::Y,
    swR::Y,
) where {Y<:AV}

    primL = conserve_prim(wL, γ)
    primR = conserve_prim(wR, γ)

    Mu1, Mxi1, MuL1, MuR1 = gauss_moments(primL, inK)
    Mu2, Mxi2, MuL2, MuR2 = gauss_moments(primR, inK)

    w =
        primL[1] .* moments_conserve(MuL1, Mxi1, 0, 0) .+
        primR[1] .* moments_conserve(MuR2, Mxi2, 0, 0)
    prim = conserve_prim(w, γ)
    tau =
        vhs_collision_time(prim, μᵣ, ω) +
        2.0 * dt * abs(primL[1] / primL[end] - primR[1] / primR[end]) /
        (primL[1] / primL[end] + primR[1] / primR[end])

    faL = pdf_slope(primL, swL, inK)
    sw = -primL[1] .* moments_conserve_slope(faL, Mu1, Mxi1, 1)
    faTL = pdf_slope(primL, sw, inK)

    faR = pdf_slope(primR, swR, inK)
    sw = -primR[1] .* moments_conserve_slope(faR, Mu2, Mxi2, 1)
    faTR = pdf_slope(primR, sw, inK)

    Mu, Mxi, MuL, MuR = gauss_moments(prim, inK)
    sw0L = (w .- (wL .- swL .* dxL)) ./ dxL
    sw0R = ((wR .+ swR .* dxR) .- w) ./ dxR
    gaL = pdf_slope(prim, sw0L, inK)
    gaR = pdf_slope(prim, sw0R, inK)
    sw =
        -prim[1] .* (
            moments_conserve_slope(gaL, MuL, Mxi, 1) .+
            moments_conserve_slope(gaR, MuR, Mxi, 1)
        )
    # ga = pdf_slope(prim, sw, inK)
    # sw = -prim[1] .* moments_conserve_slope(ga, Mu, Mxi, 1)
    gaT = pdf_slope(prim, sw, inK)

    # time-integration constants
    Mt = zeros(5)
    Mt[4] = tau * (1.0 - exp(-dt / tau))
    Mt[5] = -tau * dt * exp(-dt / tau) + tau * Mt[4]
    Mt[1] = dt - Mt[4]
    Mt[2] = -tau * Mt[1] + Mt[5]
    Mt[3] = 0.5 * dt^2 - tau * Mt[1]

    # flux related to central distribution
    Muv = moments_conserve(Mu, Mxi, 1, 0)
    MauL = moments_conserve_slope(gaL, MuL, Mxi, 2)
    MauR = moments_conserve_slope(gaR, MuR, Mxi, 2)
    # Mau = moments_conserve_slope(ga, MuR, Mxi, 2)
    MauT = moments_conserve_slope(gaT, Mu, Mxi, 1)

    fw .=
        Mt[1] .* prim[1] .* Muv .+ Mt[2] .* prim[1] .* (MauL .+ MauR) .+
        Mt[3] .* prim[1] .* MauT
    # fw .= Mt[1] .* prim[1] .* Muv .+ Mt[2] .* prim[1] .* Mau .+ Mt[3] .* prim[1] .* MauT

    # flux related to upwind distribution
    MuvL = moments_conserve(MuL1, Mxi1, 1, 0)
    MauL = moments_conserve_slope(faL, MuL1, Mxi1, 2)
    MauLT = moments_conserve_slope(faTL, MuL1, Mxi1, 1)

    MuvR = moments_conserve(MuR2, Mxi2, 1, 0)
    MauR = moments_conserve_slope(faR, MuR2, Mxi2, 2)
    MauRT = moments_conserve_slope(faTR, MuR2, Mxi2, 1)

    @. fw +=
        Mt[4] * primL[1] * MuvL - (Mt[5] + tau * Mt[4]) * primL[1] * MauL -
        tau * Mt[4] * primL[1] * MauLT + Mt[4] * primR[1] * MuvR -
        (Mt[5] + tau * Mt[4]) * primR[1] * MauR - tau * Mt[4] * primR[1] * MauRT
    # @. fw += Mt[4] * primL[1] * MuvL + Mt[4] * primR[1] * MuvR

    return nothing

end

"""
$(SIGNATURES)

Gas mixture
"""
function flux_gks!(
    fw::AM,
    wL::Y,
    wR::Y,
    inK,
    γ,
    mi,
    ni,
    me,
    ne,
    Kn,
    dt,
    dxL,
    dxR,
    swL::Y,
    swR::Y,
) where {Y<:AM}

    primL = mixture_conserve_prim(wL, γ)
    primR = mixture_conserve_prim(wR, γ)

    Mu1, Mxi1, MuL1, MuR1 = mixture_gauss_moments(primL, inK)
    Mu2, Mxi2, MuL2, MuR2 = mixture_gauss_moments(primR, inK)
    Muv1 = mixture_moments_conserve(MuL1, Mxi1, 0, 0)
    Muv2 = mixture_moments_conserve(MuR2, Mxi2, 0, 0)

    faL = mixture_pdf_slope(primL, swL, inK)
    mm = mixture_moments_conserve_slope(faL, Mu1, Mxi1, 1)
    sw = similar(swL)
    for j in axes(sw, 2)
        @. sw[:, j] = -primL[1, j] * mm[:, j]
    end
    faTL = mixture_pdf_slope(primL, sw, inK)

    faR = mixture_pdf_slope(primR, swR, inK)
    mm = mixture_moments_conserve_slope(faR, Mu2, Mxi2, 1)
    for j in axes(sw, 2)
        @. sw[:, j] = -primR[1, j] * mm[:, j]
    end
    faTR = mixture_pdf_slope(primR, sw, inK)

    w = similar(wL)
    for j in axes(w, 2)
        @. w[:, j] = primL[1, j] * Muv1[:, j] + primR[1, j] * Muv2[:, j]
    end
    prim = mixture_conserve_prim(w, γ)
    tau = aap_hs_collision_time(prim, mi, ni, me, ne, Kn)
    prim = aap_hs_prim(prim, tau, mi, ni, me, ne, Kn)

    Mu, Mxi, MuL, MuR = mixture_gauss_moments(prim, inK)
    sw0L = (w .- (wL .- swL .* dxL)) ./ dxL
    sw0R = ((wR .+ swR .* dxR) .- w) ./ dxR
    gaL = mixture_pdf_slope(prim, sw0L, inK)
    gaR = mixture_pdf_slope(prim, sw0R, inK)
    mmL = mixture_moments_conserve_slope(gaL, MuL, Mxi, 1)
    mmR = mixture_moments_conserve_slope(gaR, MuR, Mxi, 1)
    for j in axes(sw, 2)
        @. sw[:, j] = -prim[1, j] * (mmL[:, j] + mmR[:, j])
    end
    gaT = mixture_pdf_slope(prim, sw, inK)

    Mt = zeros(5, 2)
    for j in axes(Mt, 2)
        Mt[4, j] = tau[j] * (1.0 - exp(-dt / tau[j]))
        Mt[5, j] = -tau[j] * dt * exp(-dt / tau[j]) + tau[j] * Mt[4]
        Mt[1, j] = dt - Mt[4, j]
        Mt[2, j] = -tau[j] * Mt[1, j] + Mt[5, j]
        Mt[3, j] = 0.5 * dt^2 - tau[j] * Mt[1, j]
    end

    # flux related to central distribution
    Muv = mixture_moments_conserve(Mu, Mxi, 1, 0)
    MauL = mixture_moments_conserve_slope(gaL, MuL, Mxi, 2)
    MauR = mixture_moments_conserve_slope(gaR, MuR, Mxi, 2)
    MauT = mixture_moments_conserve_slope(gaT, Mu, Mxi, 1)
    for j in axes(fw, 2)
        @. fw[:, j] =
            Mt[1, j] * prim[1, j] * Muv[:, j] +
            Mt[2, j] * prim[1, j] * (MauL[:, j] + MauR[:, j]) +
            Mt[3, j] * prim[1, j] * MauT[:, j]
    end

    # flux related to upwind distribution
    MuvL = mixture_moments_conserve(MuL1, Mxi1, 1, 0)
    MauL = mixture_moments_conserve_slope(faL, MuL1, Mxi1, 2)
    MauLT = mixture_moments_conserve_slope(faTL, MuL1, Mxi1, 1)

    MuvR = mixture_moments_conserve(MuR2, Mxi2, 1, 0)
    MauR = mixture_moments_conserve_slope(faR, MuR2, Mxi2, 2)
    MauRT = mixture_moments_conserve_slope(faTR, MuR2, Mxi2, 1)

    for j in axes(fw, 2)
        @. fw[:, j] +=
            Mt[4, j] * primL[1, j] * MuvL[:, j] -
            (Mt[5, j] + tau[j] * Mt[4, j]) * primL[1, j] * MauL[:, j] -
            tau[j] * Mt[4, j] * primL[1, j] * MauLT[:, j] +
            Mt[4, j] * primR[1, j] * MuvR[:, j] -
            (Mt[5, j] + tau[j] * Mt[4, j]) * primR[1, j] * MauR[:, j] -
            tau[j] * Mt[4, j] * primR[1, j] * MauRT[:, j]
    end

    return nothing

end

"""
$(SIGNATURES)

1D1F1V
"""
function flux_gks!(
    fw::AV,
    fh::AV,
    wL::T3,
    wR::T3,
    u::AV,
    inK,
    γ,
    visRef,
    visIdx,
    dt,
    dxL,
    dxR,
    swL::T3,
    swR::T3,
) where {T3<:AV}

    primL = conserve_prim(wL, γ)
    primR = conserve_prim(wR, γ)

    Mu1, Mxi1, MuL1, _ = gauss_moments(primL, inK)
    Mu2, Mxi2, _, MuR2 = gauss_moments(primR, inK)

    w =
        primL[1] .* moments_conserve(MuL1, Mxi1, 0, 0) .+
        primR[1] .* moments_conserve(MuR2, Mxi2, 0, 0)
    prim = conserve_prim(w, γ)
    tau =
        vhs_collision_time(prim, visRef, visIdx) +
        2.0 * dt * abs(primL[1] / primL[end] - primR[1] / primR[end]) /
        (primL[1] / primL[end] + primR[1] / primR[end])

    faL = pdf_slope(primL, swL, inK)
    sw = -primL[1] .* moments_conserve_slope(faL, Mu1, Mxi1, 1)
    faTL = pdf_slope(primL, sw, inK)

    faR = pdf_slope(primR, swR, inK)
    sw = -primR[1] .* moments_conserve_slope(faR, Mu2, Mxi2, 1)
    faTR = pdf_slope(primR, sw, inK)

    Mu, Mxi, MuL, MuR = gauss_moments(prim, inK)
    sw0L = (w .- (wL .- swL .* dxL)) ./ dxL
    sw0R = ((wR .+ swR .* dxR) .- w) ./ dxR
    gaL = pdf_slope(prim, sw0L, inK)
    gaR = pdf_slope(prim, sw0R, inK)
    sw =
        -prim[1] .* (
            moments_conserve_slope(gaL, MuL, Mxi, 1) .+
            moments_conserve_slope(gaR, MuR, Mxi, 1)
        )
    # ga = pdf_slope(prim, sw, inK)
    # sw = -prim[1] .* moments_conserve_slope(ga, Mu, Mxi, 1)
    gaT = pdf_slope(prim, sw, inK)

    # time-integration constants
    Mt = zeros(5)
    Mt[4] = tau * (1.0 - exp(-dt / tau))
    Mt[5] = -tau * dt * exp(-dt / tau) + tau * Mt[4]
    Mt[1] = dt - Mt[4]
    Mt[2] = -tau * Mt[1] + Mt[5]
    Mt[3] = 0.5 * dt^2 - tau * Mt[1]

    # flux related to central distribution
    Muv = moments_conserve(Mu, Mxi, 1, 0)
    MauL = moments_conserve_slope(gaL, MuL, Mxi, 2)
    MauR = moments_conserve_slope(gaR, MuR, Mxi, 2)
    # Mau = moments_conserve_slope(ga, MuR, Mxi, 2)
    MauT = moments_conserve_slope(gaT, Mu, Mxi, 1)

    fw .=
        Mt[1] .* prim[1] .* Muv .+ Mt[2] .* prim[1] .* (MauL .+ MauR) .+
        Mt[3] .* prim[1] .* MauT
    # fw .= Mt[1] .* prim[1] .* Muv .+ Mt[2] .* prim[1] .* Mau .+ Mt[3] .* prim[1] .* MauT

    # flux related to upwind distribution
    MuvL = moments_conserve(MuL1, Mxi1, 1, 0)
    MauL = moments_conserve_slope(faL, MuL1, Mxi1, 2)
    MauLT = moments_conserve_slope(faTL, MuL1, Mxi1, 1)

    MuvR = moments_conserve(MuR2, Mxi2, 1, 0)
    MauR = moments_conserve_slope(faR, MuR2, Mxi2, 2)
    MauRT = moments_conserve_slope(faTR, MuR2, Mxi2, 1)

    @. fw +=
        Mt[4] * primL[1] * MuvL - (Mt[5] + tau * Mt[4]) * primL[1] * MauL -
        tau * Mt[4] * primL[1] * MauLT + Mt[4] * primR[1] * MuvR -
        (Mt[5] + tau * Mt[4]) * primR[1] * MauR - tau * Mt[4] * primR[1] * MauRT
    # @. fw += Mt[4] * primL[1] * MuvL + Mt[4] * primR[1] * MuvR

    #--- fluxes of distribution functions ---#
    δ = heaviside.(u)
    HL = maxwellian(u, primL)
    HR = maxwellian(u, primR)
    H = maxwellian(u, prim)

    @. fh =
        Mt[1] * u * H +
        Mt[2] * u^2 * (gaL[1] * H + gaL[2] * u * H + 0.5 * gaL[3] * (u^2 * H)) * δ +
        Mt[2] * u^2 * (gaR[1] * H + gaR[2] * u * H + 0.5 * gaR[3] * (u^2 * H)) * (1.0 - δ) +
        Mt[3] * u * (gaT[1] * H + gaT[2] * u * H + 0.5 * gaT[3] * (u^2 * H)) +
        Mt[4] * u * HL * δ -
        (Mt[5] + tau * Mt[4]) *
        u^2 *
        (faL[1] * HL + faL[2] * u * HL + 0.5 * faL[3] * (u^2 * HL)) *
        δ + Mt[4] * u * HR * (1.0 - δ) -
        (Mt[5] + tau * Mt[4]) *
        u^2 *
        (faR[1] * HR + faR[2] * u * HR + 0.5 * faR[3] * (u^2 * HR)) *
        (1.0 - δ) -
        tau *
        Mt[4] *
        u *
        (faTL[1] * HL + faTL[2] * u * HL + 0.5 * faTL[3] * (u^2 * HL)) *
        δ -
        tau *
        Mt[4] *
        u *
        (faTR[1] * HR + faTR[2] * u * HR + 0.5 * faTR[3] * (u^2 * HR)) *
        (1.0 - δ)

    return nothing

end

"""
$(SIGNATURES)

1D2F1V
"""
function flux_gks!(
    fw::AV,
    fh::T2,
    fb::T2,
    wL::T3,
    wR::T3,
    u::AV,
    inK,
    γ,
    visRef,
    visIdx,
    dt,
    dxL,
    dxR,
    swL::T3,
    swR::T3,
) where {T2<:AV,T3<:AV}

    primL = conserve_prim(wL, γ)
    primR = conserve_prim(wR, γ)

    Mu1, Mxi1, MuL1, MuR1 = gauss_moments(primL, inK)
    Mu2, Mxi2, MuL2, MuR2 = gauss_moments(primR, inK)

    w =
        primL[1] .* moments_conserve(MuL1, Mxi1, 0, 0) .+
        primR[1] .* moments_conserve(MuR2, Mxi2, 0, 0)
    prim = conserve_prim(w, γ)
    tau =
        vhs_collision_time(prim, visRef, visIdx) +
        2.0 * dt * abs(primL[1] / primL[end] - primR[1] / primR[end]) /
        (primL[1] / primL[end] + primR[1] / primR[end])

    faL = pdf_slope(primL, swL, inK)
    sw = -primL[1] .* moments_conserve_slope(faL, Mu1, Mxi1, 1)
    faTL = pdf_slope(primL, sw, inK)

    faR = pdf_slope(primR, swR, inK)
    sw = -primR[1] .* moments_conserve_slope(faR, Mu2, Mxi2, 1)
    faTR = pdf_slope(primR, sw, inK)

    Mu, Mxi, MuL, MuR = gauss_moments(prim, inK)
    sw0L = (w .- (wL .- swL .* dxL)) ./ dxL
    sw0R = ((wR .+ swR .* dxR) .- w) ./ dxR
    gaL = pdf_slope(prim, sw0L, inK)
    gaR = pdf_slope(prim, sw0R, inK)
    sw =
        -prim[1] .* (
            moments_conserve_slope(gaL, MuL, Mxi, 1) .+
            moments_conserve_slope(gaR, MuR, Mxi, 1)
        )
    # ga = pdf_slope(prim, sw, inK)
    # sw = -prim[1] .* moments_conserve_slope(ga, Mu, Mxi, 1)
    gaT = pdf_slope(prim, sw, inK)

    # time-integration constants
    Mt = zeros(5)
    Mt[4] = tau * (1.0 - exp(-dt / tau))
    Mt[5] = -tau * dt * exp(-dt / tau) + tau * Mt[4]
    Mt[1] = dt - Mt[4]
    Mt[2] = -tau * Mt[1] + Mt[5]
    Mt[3] = 0.5 * dt^2 - tau * Mt[1]

    # flux related to central distribution
    Muv = moments_conserve(Mu, Mxi, 1, 0)
    MauL = moments_conserve_slope(gaL, MuL, Mxi, 2)
    MauR = moments_conserve_slope(gaR, MuR, Mxi, 2)
    # Mau = moments_conserve_slope(ga, MuR, Mxi, 2)
    MauT = moments_conserve_slope(gaT, Mu, Mxi, 1)

    fw .=
        Mt[1] .* prim[1] .* Muv .+ Mt[2] .* prim[1] .* (MauL .+ MauR) .+
        Mt[3] .* prim[1] .* MauT
    # fw .= Mt[1] .* prim[1] .* Muv .+ Mt[2] .* prim[1] .* Mau .+ Mt[3] .* prim[1] .* MauT

    # flux related to upwind distribution
    MuvL = moments_conserve(MuL1, Mxi1, 1, 0)
    MauL = moments_conserve_slope(faL, MuL1, Mxi1, 2)
    MauLT = moments_conserve_slope(faTL, MuL1, Mxi1, 1)

    MuvR = moments_conserve(MuR2, Mxi2, 1, 0)
    MauR = moments_conserve_slope(faR, MuR2, Mxi2, 2)
    MauRT = moments_conserve_slope(faTR, MuR2, Mxi2, 1)

    @. fw +=
        Mt[4] * primL[1] * MuvL - (Mt[5] + tau * Mt[4]) * primL[1] * MauL -
        tau * Mt[4] * primL[1] * MauLT + Mt[4] * primR[1] * MuvR -
        (Mt[5] + tau * Mt[4]) * primR[1] * MauR - tau * Mt[4] * primR[1] * MauRT
    # @. fw += Mt[4] * primL[1] * MuvL + Mt[4] * primR[1] * MuvR

    #--- fluxes of distribution functions ---#
    δ = heaviside.(u)

    HL = maxwellian(u, primL)
    BL = HL .* inK ./ (2.0 * primL[end])
    HR = maxwellian(u, primR)
    BR = HR .* inK ./ (2.0 * primR[end])

    H = maxwellian(u, prim)
    B = H .* inK ./ (2.0 * prim[end])

    @. fh =
        Mt[1] * u * H +
        Mt[2] * u^2 * (gaL[1] * H + gaL[2] * u * H + 0.5 * gaL[3] * (u^2 * H + B)) * δ +
        Mt[2] *
        u^2 *
        (gaR[1] * H + gaR[2] * u * H + 0.5 * gaR[3] * (u^2 * H + B)) *
        (1.0 - δ) +
        Mt[3] * u * (gaT[1] * H + gaT[2] * u * H + 0.5 * gaT[3] * (u^2 * H + B)) +
        Mt[4] * u * HL * δ -
        (Mt[5] + tau * Mt[4]) *
        u^2 *
        (faL[1] * HL + faL[2] * u * HL + 0.5 * faL[3] * (u^2 * HL + BL)) *
        δ + Mt[4] * u * HR * (1.0 - δ) -
        (Mt[5] + tau * Mt[4]) *
        u^2 *
        (faR[1] * HR + faR[2] * u * HR + 0.5 * faR[3] * (u^2 * HR + BR)) *
        (1.0 - δ) -
        tau *
        Mt[4] *
        u *
        (faTL[1] * HL + faTL[2] * u * HL + 0.5 * faTL[3] * (u^2 * HL + BL)) *
        δ -
        tau *
        Mt[4] *
        u *
        (faTR[1] * HR + faTR[2] * u * HR + 0.5 * faTR[3] * (u^2 * HR + BR)) *
        (1.0 - δ)

    @. fb =
        Mt[1] * u * B +
        Mt[2] *
        u^2 *
        (gaL[1] * B + gaL[2] * u * B + 0.5 * gaL[3] * (u^2 * B + Mxi[2] * H)) *
        δ +
        Mt[2] *
        u^2 *
        (gaR[1] * B + gaR[2] * u * B + 0.5 * gaR[3] * (u^2 * B + Mxi[2] * H)) *
        (1.0 - δ) +
        Mt[3] * u * (gaT[1] * H + gaT[2] * u * H + 0.5 * gaT[3] * (u^2 * B + Mxi[2] * H)) +
        Mt[4] * u * BL * δ -
        (Mt[5] + tau * Mt[4]) *
        u^2 *
        (faL[1] * BL + faL[2] * u * BL + 0.5 * faL[3] * (u^2 * BL + Mxi[2] * HL)) *
        δ + Mt[4] * u * BR * (1.0 - δ) -
        (Mt[5] + tau * Mt[4]) *
        u^2 *
        (faR[1] * BR + faR[2] * u * BR + 0.5 * faR[3] * (u^2 * BR + Mxi[2] * HR)) *
        (1.0 - δ) -
        tau *
        Mt[4] *
        u *
        (faTL[1] * BL + faTL[2] * u * BL + 0.5 * faTL[3] * (u^2 * BL + Mxi[2] * HL)) *
        δ -
        tau *
        Mt[4] *
        u *
        (faTR[1] * BR + faTR[2] * u * BR + 0.5 * faTR[3] * (u^2 * BR + Mxi[2] * HR)) *
        (1.0 - δ)

    return nothing

end

# ------------------------------------------------------------
# 2D & 3D
# ------------------------------------------------------------

"""
$(SIGNATURES)

2D & 2D
"""
function flux_gks!(
    fw::AV,
    wL::Y,
    wR::Y,
    inK,
    γ,
    μᵣ,
    ω,
    dt,
    dxL,
    dxR,
    dy,
    swL::Y,
    swR::Y,
) where {Y<:AV}

    primL = conserve_prim(wL, γ)
    primR = conserve_prim(wR, γ)

    Mu1, Mv1, Mw1, MuL1, MuR1 = gauss_moments(primL, inK)
    Mu2, Mv2, Mw2, MuL2, MuR2 = gauss_moments(primR, inK)

    w =
        primL[1] .* moments_conserve(MuL1, Mv1, Mw1, 0, 0, 0) .+
        primR[1] .* moments_conserve(MuR2, Mv2, Mw2, 0, 0, 0)
    prim = conserve_prim(w, γ)
    tau =
        vhs_collision_time(prim, μᵣ, ω) +
        2.0 * dt * abs(primL[1] / primL[end] - primR[1] / primR[end]) /
        (primL[1] / primL[end] + primR[1] / primR[end])

    faL = pdf_slope(primL, swL, inK)
    sw = -primL[1] .* moments_conserve_slope(faL, Mu1, Mv1, Mw1, 1, 0, 0)
    faTL = pdf_slope(primL, sw, inK)

    faR = pdf_slope(primR, swR, inK)
    sw = -primR[1] .* moments_conserve_slope(faR, Mu2, Mv1, Mw2, 1, 0, 0)
    faTR = pdf_slope(primR, sw, inK)

    Mu, Mv, Mw, MuL, MuR = gauss_moments(prim, inK)
    sw0L = (w .- (wL .- swL .* dxL)) ./ dxL
    sw0R = ((wR .+ swR .* dxR) .- w) ./ dxR
    gaL = pdf_slope(prim, sw0L, inK)
    gaR = pdf_slope(prim, sw0R, inK)
    sw =
        -prim[1] .* (
            moments_conserve_slope(gaL, MuL, Mv, Mw, 1, 0, 0) .+
            moments_conserve_slope(gaR, MuR, Mv, Mw, 1, 0, 0)
        )
    # ga = pdf_slope(prim, sw, inK)
    # sw =  -prim[1] .* moments_conserve_slope(ga, Mu, Mv, Mw, 1, 0, 0)
    gaT = pdf_slope(prim, sw, inK)

    # time-integration constants
    Mt = zeros(5)
    Mt[4] = tau * (1.0 - exp(-dt / tau))
    Mt[5] = -tau * dt * exp(-dt / tau) + tau * Mt[4]
    Mt[1] = dt - Mt[4]
    Mt[2] = -tau * Mt[1] + Mt[5]
    Mt[3] = 0.5 * dt^2 - tau * Mt[1]

    # flux related to central distribution
    Muv = moments_conserve(Mu, Mv, Mw, 1, 0, 0)
    MauL = moments_conserve_slope(gaL, MuL, Mv, Mw, 2, 0, 0)
    MauR = moments_conserve_slope(gaR, MuR, Mv, Mw, 2, 0, 0)
    # Mau = moments_conserve_slope(ga, Mu, Mv, Mw, 2, 0, 0)
    MauT = moments_conserve_slope(gaT, Mu, Mv, Mw, 1, 0, 0)

    fw .=
        Mt[1] .* prim[1] .* Muv .+ Mt[2] .* prim[1] .* (MauL .+ MauR) .+
        Mt[3] .* prim[1] .* MauT
    # fw .= Mt[1] .* prim[1] .* Muv .+ Mt[2] .* prim[1] .* Mau .+ Mt[3] .* prim[1] .* MauT

    # flux related to upwind distribution
    MuvL = moments_conserve(MuL1, Mv1, Mw1, 1, 0, 0)
    MauL = moments_conserve_slope(faL, MuL1, Mv1, Mw1, 2, 0, 0)
    MauLT = moments_conserve_slope(faTL, MuL1, Mv1, Mw1, 1, 0, 0)

    MuvR = moments_conserve(MuR2, Mv2, Mw2, 1, 0, 0)
    MauR = moments_conserve_slope(faR, MuR2, Mv2, Mw2, 2, 0, 0)
    MauRT = moments_conserve_slope(faTR, MuR2, Mv2, Mw2, 1, 0, 0)

    @. fw +=
        Mt[4] * primL[1] * MuvL - (Mt[5] + tau * Mt[4]) * primL[1] * MauL -
        tau * Mt[4] * primL[1] * MauLT + Mt[4] * primR[1] * MuvR -
        (Mt[5] + tau * Mt[4]) * primR[1] * MauR - tau * Mt[4] * primR[1] * MauRT
    # @. fw += Mt[4] * primL[1] * MuvL + Mt[4] * primR[1] * MuvR

    fw .*= dy

    return nothing

end

"""
$(SIGNATURES)

Mixture
"""
function flux_gks!(
    fw::AM,
    wL::Y,
    wR::Y,
    inK,
    γ,
    mi,
    ni,
    me,
    ne,
    Kn,
    dt,
    dxL,
    dxR,
    len,
    swL::Y,
    swR::Y,
) where {Y<:AM}

    primL = mixture_conserve_prim(wL, γ)
    primR = mixture_conserve_prim(wR, γ)

    Mu1, Mv1, Mw1, MuL1, MuR1 = mixture_gauss_moments(primL, inK)
    Mu2, Mv2, Mw2, MuL2, MuR2 = mixture_gauss_moments(primR, inK)

    w =
        primL[1] .* mixture_moments_conserve(MuL1, Mv1, Mw1, 0, 0, 0) .+
        primR[1] .* mixture_moments_conserve(MuR2, Mv2, Mw2, 0, 0, 0)
    prim = mixture_conserve_prim(w, γ)
    tau = aap_hs_collision_time(prim, mi, ni, me, ne, Kn)
    for i in eachindex(tau)
        tau[i] +=
            2.0 * dt * abs(primL[1, i] / primL[end, i] - primR[1, i] / primR[end, i]) /
            (primL[1, i] / primL[end, i] + primR[1, i] / primR[end, i])
    end
    prim = aap_hs_prim(prim, tau, mi, ni, me, ne, Kn) # pseudo primitive variables

    faL = mixture_pdf_slope(primL, swL, inK)
    sw = mixture_moments_conserve_slope(faL, Mu1, Mv1, Mw1, 1, 0)
    for j in axes(sw, 2)
        sw[:, j] .*= -primL[1, j]
    end
    faTL = mixture_pdf_slope(primL, sw, inK)
    faR = mixture_pdf_slope(primR, swR, inK)
    sw = mixture_moments_conserve_slope(faR, Mu2, Mv2, Mw2, 1, 0)
    for j in axes(sw, 2)
        sw[:, j] .*= -primR[1, j]
    end
    faTR = mixture_pdf_slope(primR, sw, inK)

    Mu, Mv, Mw, MuL, MuR = mixture_gauss_moments(prim, inK)
    sw0L = (w .- (wL .- swL .* dxL)) ./ dxL
    sw0R = ((wR .+ swR .* dxR) .- w) ./ dxR
    gaL = mixture_pdf_slope(prim, sw0L, inK)
    gaR = mixture_pdf_slope(prim, sw0R, inK)
    sw =
        mixture_moments_conserve_slope(gaL, MuL, Mv, Mw, 1, 0) .+
        mixture_moments_conserve_slope(gaR, MuR, Mv, Mw, 1, 0)
    for j in axes(sw, 2)
        sw[:, j] .*= -prim[1, j]
    end
    gaT = mixture_pdf_slope(prim, sw, inK)

    # time-integration constants
    Mt = zeros(5, axes(fw, 2))
    for j in axes(Mt, 2)
        Mt[4, j] = tau[j] * (1.0 - exp(-dt / tau[j]))
        Mt[5, j] = -tau[j] * dt * exp(-dt / tau[j]) + tau[j] * Mt[4]
        Mt[1, j] = dt - Mt[4, j]
        Mt[2, j] = -tau[j] * Mt[1, j] + Mt[5, j]
        Mt[3, j] = 0.5 * dt^2 - tau[j] * Mt[1, j]
    end

    # flux related to central distribution
    Muv = mixture_moments_conserve(Mu, Mv, Mw, 1, 0, 0)
    MauL = mixture_moments_conserve_slope(gaL, MuL, Mv, Mw, 2, 0)
    MauR = mixture_moments_conserve_slope(gaR, MuR, Mv, Mw, 2, 0)
    MauT = mixture_moments_conserve_slope(gaT, Mu, Mv, Mw, 1, 0)

    for j in axes(fw, 2)
        @. fw[:, j] =
            Mt[1, j] * prim[1, j] * Muv[:, j] +
            Mt[2, j] * prim[1, j] * (MauL[:, j] + MauR[:, j]) +
            Mt[3, j] * prim[1, j] * MauT[:, j]
    end

    # flux related to upwind distribution
    MuvL = mixture_moments_conserve(MuL1, Mv1, Mw1, 1, 0, 0)
    MauL = mixture_moments_conserve_slope(faL, MuL1, Mv1, Mw1, 2, 0)
    MauLT = mixture_moments_conserve_slope(faTL, MuL1, Mv1, Mw1, 1, 0)

    MuvR = mixture_moments_conserve(MuR2, Mv2, Mw2, 1, 0, 0)
    MauR = mixture_moments_conserve_slope(faR, MuR2, Mv2, Mw2, 2, 0)
    MauRT = mixture_moments_conserve_slope(faTR, MuR2, Mv2, Mw2, 1, 0)

    for j in axes(fw, 2)
        @. fw[:, j] +=
            Mt[4, j] * primL[1, j] * MuvL[:, j] -
            (Mt[5, j] + tau[j] * Mt[4, j]) * primL[1, j] * MauL[:, j] -
            tau[j] * Mt[4, j] * primL[1, j] * MauLT[:, j] +
            Mt[4, j] * primR[1, j] * MuvR[:, j] -
            (Mt[5, j] + tau[j] * Mt[4, j]) * primR[1, j] * MauR[:, j] -
            tau[j] * Mt[4, j] * primR[1, j] * MauRT[:, j]
    end

    @. fw .* len

    return nothing

end

"""
$(SIGNATURES)

2D1F2V
"""
function flux_gks!(
    fw::AV,
    ff::AM,
    wL::T3,
    wR::T3,
    u::T4,
    v::T4,
    inK,
    γ,
    μᵣ,
    ω,
    dt,
    dxL,
    dxR,
    dy,
    swL::T3,
    swR::T3,
) where {T3<:AV,T4<:AM}

    primL = conserve_prim(wL, γ)
    primR = conserve_prim(wR, γ)

    Mu1, Mv1, Mxi1, MuL1, MuR1 = gauss_moments(primL, inK)
    Mu2, Mv2, Mxi2, MuL2, MuR2 = gauss_moments(primR, inK)

    w =
        primL[1] .* moments_conserve(MuL1, Mv1, Mxi1, 0, 0, 0) .+
        primR[1] .* moments_conserve(MuR2, Mv2, Mxi2, 0, 0, 0)
    prim = conserve_prim(w, γ)
    tau =
        vhs_collision_time(prim, μᵣ, ω) +
        2.0 * dt * abs(primL[1] / primL[end] - primR[1] / primR[end]) /
        (primL[1] / primL[end] + primR[1] / primR[end])

    faL = pdf_slope(primL, swL, inK)
    sw = -primL[1] .* moments_conserve_slope(faL, Mu1, Mv1, Mxi1, 1, 0)
    faTL = pdf_slope(primL, sw, inK)

    faR = pdf_slope(primR, swR, inK)
    sw = -primR[1] .* moments_conserve_slope(faR, Mu2, Mv1, Mxi2, 1, 0)
    faTR = pdf_slope(primR, sw, inK)

    Mu, Mv, Mxi, MuL, MuR = gauss_moments(prim, inK)
    sw0L = (w .- (wL .- swL .* dxL)) ./ dxL
    sw0R = ((wR .+ swR .* dxR) .- w) ./ dxR
    gaL = pdf_slope(prim, sw0L, inK)
    gaR = pdf_slope(prim, sw0R, inK)
    sw =
        -prim[1] .* (
            moments_conserve_slope(gaL, MuL, Mv, Mxi, 1, 0) .+
            moments_conserve_slope(gaR, MuR, Mv, Mxi, 1, 0)
        )
    # ga = pdf_slope(prim, sw, inK)
    # sw = -prim[1] .* moments_conserve_slope(ga, Mu, Mv, Mxi, 1, 0)
    gaT = pdf_slope(prim, sw, inK)

    # time-integration constants
    Mt = zeros(5)
    Mt[4] = tau * (1.0 - exp(-dt / tau))
    Mt[5] = -tau * dt * exp(-dt / tau) + tau * Mt[4]
    Mt[1] = dt - Mt[4]
    Mt[2] = -tau * Mt[1] + Mt[5]
    Mt[3] = 0.5 * dt^2 - tau * Mt[1]

    # flux related to central distribution
    Muv = moments_conserve(Mu, Mv, Mxi, 1, 0, 0)
    MauL = moments_conserve_slope(gaL, MuL, Mv, Mxi, 2, 0)
    MauR = moments_conserve_slope(gaR, MuR, Mv, Mxi, 2, 0)
    # Mau = moments_conserve_slope(ga, Mu, Mv, Mxi, 2, 0)
    MauT = moments_conserve_slope(gaT, Mu, Mv, Mxi, 1, 0)

    fw .=
        Mt[1] .* prim[1] .* Muv .+ Mt[2] .* prim[1] .* (MauL .+ MauR) .+
        Mt[3] .* prim[1] .* MauT
    # fw .= Mt[1] .* prim[1] .* Muv .+ Mt[2] .* prim[1] .* Mau .+ Mt[3] .* prim[1] .* MauT

    # flux related to upwind distribution
    MuvL = moments_conserve(MuL1, Mv1, Mxi1, 1, 0, 0)
    MauL = moments_conserve_slope(faL, MuL1, Mv1, Mxi1, 2, 0)
    MauLT = moments_conserve_slope(faTL, MuL1, Mv1, Mxi1, 1, 0)

    MuvR = moments_conserve(MuR2, Mv2, Mxi2, 1, 0, 0)
    MauR = moments_conserve_slope(faR, MuR2, Mv2, Mxi2, 2, 0)
    MauRT = moments_conserve_slope(faTR, MuR2, Mv2, Mxi2, 1, 0)

    @. fw +=
        Mt[4] * primL[1] * MuvL - (Mt[5] + tau * Mt[4]) * primL[1] * MauL -
        tau * Mt[4] * primL[1] * MauLT + Mt[4] * primR[1] * MuvR -
        (Mt[5] + tau * Mt[4]) * primR[1] * MauR - tau * Mt[4] * primR[1] * MauRT
    # @. fw += Mt[4] * primL[1] * MuvL + Mt[4] * primR[1] * MuvR

    fw .*= dy

    #--- fluxes of distribution functions ---#
    δ = heaviside.(u)
    HL = maxwellian(u, v, primL)
    HR = maxwellian(u, v, primR)
    H = maxwellian(u, v, prim)

    @. ff =
        Mt[1] * u * H +
        Mt[2] *
        u^2 *
        (gaL[1] * H + gaL[2] * u * H + gaL[3] * v * H + 0.5 * gaL[4] * (u^2 + v^2) * H) *
        δ +
        Mt[2] *
        u^2 *
        (gaR[1] * H + gaR[2] * u * H + gaR[3] * v * H + 0.5 * gaR[4] * (u^2 + v^2) * H) *
        (1.0 - δ) +
        Mt[3] *
        u *
        (gaT[1] * H + gaT[2] * u * H + gaT[3] * v * H + 0.5 * gaT[4] * (u^2 + v^2) * H) +
        Mt[4] * u * HL * δ -
        (Mt[5] + tau * Mt[4]) *
        u^2 *
        (
            faL[1] * HL +
            faL[2] * u * HL +
            faL[3] * v * HL +
            0.5 * faL[4] * (u^2 + v^2) * HL
        ) *
        δ + Mt[4] * u * HR * (1.0 - δ) -
        (Mt[5] + tau * Mt[4]) *
        u^2 *
        (
            faR[1] * HR +
            faR[2] * u * HR +
            faR[3] * v * HL +
            0.5 * faR[4] * (u^2 + v^2) * HR
        ) *
        (1.0 - δ) -
        tau *
        Mt[4] *
        u *
        (
            faTL[1] * HL +
            faTL[2] * u * HL +
            faTL[3] * v * HL +
            0.5 * faTL[4] * (u^2 + v^2) * HL
        ) *
        δ -
        tau *
        Mt[4] *
        u *
        (
            faTR[1] * HR +
            faTR[2] * u * HR +
            faTR[3] * v * HR +
            0.5 * faTR[4] * (u^2 + v^2) * HR
        ) *
        (1.0 - δ)

    ff .*= dy

    return nothing

end

"""
$(SIGNATURES)

2D2F2V
"""
function flux_gks!(
    fw::AV,
    fh::T2,
    fb::T2,
    wL::T3,
    wR::T3,
    u::T4,
    v::T4,
    inK,
    γ,
    μᵣ,
    ω,
    dt,
    dxL,
    dxR,
    dy,
    swL::T3,
    swR::T3,
) where {T2<:AM,T3<:AV,T4<:AM}

    primL = conserve_prim(wL, γ)
    primR = conserve_prim(wR, γ)

    Mu1, Mv1, Mxi1, MuL1, MuR1 = gauss_moments(primL, inK)
    Mu2, Mv2, Mxi2, MuL2, MuR2 = gauss_moments(primR, inK)

    w =
        primL[1] .* moments_conserve(MuL1, Mv1, Mxi1, 0, 0, 0) .+
        primR[1] .* moments_conserve(MuR2, Mv2, Mxi2, 0, 0, 0)
    prim = conserve_prim(w, γ)
    tau =
        vhs_collision_time(prim, μᵣ, ω) +
        2.0 * dt * abs(primL[1] / primL[end] - primR[1] / primR[end]) /
        (primL[1] / primL[end] + primR[1] / primR[end])

    faL = pdf_slope(primL, swL, inK)
    sw = -primL[1] .* moments_conserve_slope(faL, Mu1, Mv1, Mxi1, 1, 0)
    faTL = pdf_slope(primL, sw, inK)

    faR = pdf_slope(primR, swR, inK)
    sw = -primR[1] .* moments_conserve_slope(faR, Mu2, Mv1, Mxi2, 1, 0)
    faTR = pdf_slope(primR, sw, inK)

    Mu, Mv, Mxi, MuL, MuR = gauss_moments(prim, inK)
    sw0L = (w .- (wL .- swL .* dxL)) ./ dxL
    sw0R = ((wR .+ swR .* dxR) .- w) ./ dxR
    gaL = pdf_slope(prim, sw0L, inK)
    gaR = pdf_slope(prim, sw0R, inK)
    sw =
        -prim[1] .* (
            moments_conserve_slope(gaL, MuL, Mv, Mxi, 1, 0) .+
            moments_conserve_slope(gaR, MuR, Mv, Mxi, 1, 0)
        )
    # ga = pdf_slope(prim, sw, inK)
    # sw = -prim[1] .* moments_conserve_slope(ga, Mu, Mv, Mxi, 1, 0)
    gaT = pdf_slope(prim, sw, inK)

    # time-integration constants
    Mt = zeros(5)
    Mt[4] = tau * (1.0 - exp(-dt / tau))
    Mt[5] = -tau * dt * exp(-dt / tau) + tau * Mt[4]
    Mt[1] = dt - Mt[4]
    Mt[2] = -tau * Mt[1] + Mt[5]
    Mt[3] = 0.5 * dt^2 - tau * Mt[1]

    # flux related to central distribution
    Muv = moments_conserve(Mu, Mv, Mxi, 1, 0, 0)
    MauL = moments_conserve_slope(gaL, MuL, Mv, Mxi, 2, 0)
    MauR = moments_conserve_slope(gaR, MuR, Mv, Mxi, 2, 0)
    # Mau = moments_conserve_slope(ga, Mu, Mv, Mxi, 2, 0)
    MauT = moments_conserve_slope(gaT, Mu, Mv, Mxi, 1, 0)

    fw .=
        Mt[1] .* prim[1] .* Muv .+ Mt[2] .* prim[1] .* (MauL .+ MauR) .+
        Mt[3] .* prim[1] .* MauT
    # fw .= Mt[1] .* prim[1] .* Muv .+ Mt[2] .* prim[1] .* Mau .+ Mt[3] .* prim[1] .* MauT

    # flux related to upwind distribution
    MuvL = moments_conserve(MuL1, Mv1, Mxi1, 1, 0, 0)
    MauL = moments_conserve_slope(faL, MuL1, Mv1, Mxi1, 2, 0)
    MauLT = moments_conserve_slope(faTL, MuL1, Mv1, Mxi1, 1, 0)

    MuvR = moments_conserve(MuR2, Mv2, Mxi2, 1, 0, 0)
    MauR = moments_conserve_slope(faR, MuR2, Mv2, Mxi2, 2, 0)
    MauRT = moments_conserve_slope(faTR, MuR2, Mv2, Mxi2, 1, 0)

    @. fw +=
        Mt[4] * primL[1] * MuvL - (Mt[5] + tau * Mt[4]) * primL[1] * MauL -
        tau * Mt[4] * primL[1] * MauLT + Mt[4] * primR[1] * MuvR -
        (Mt[5] + tau * Mt[4]) * primR[1] * MauR - tau * Mt[4] * primR[1] * MauRT
    # @. fw += Mt[4] * primL[1] * MuvL + Mt[4] * primR[1] * MuvR

    fw .*= dy

    #--- fluxes of distribution functions ---#
    δ = heaviside.(u)
    HL = maxwellian(u, v, primL)
    HR = maxwellian(u, v, primR)
    H = maxwellian(u, v, prim)
    BL = HL .* inK ./ (2.0 * primL[end])
    BR = HR .* inK ./ (2.0 * primR[end])
    B = H .* inK ./ (2.0 * prim[end])

    @. fh =
        Mt[1] * u * H +
        Mt[2] *
        u^2 *
        (
            gaL[1] * H +
            gaL[2] * u * H +
            gaL[3] * v * H +
            0.5 * gaL[4] * ((u^2 + v^2) * H + B)
        ) *
        δ +
        Mt[2] *
        u^2 *
        (
            gaR[1] * H +
            gaR[2] * u * H +
            gaR[3] * v * H +
            0.5 * gaR[4] * ((u^2 + v^2) * H + B)
        ) *
        (1.0 - δ) +
        Mt[3] *
        u *
        (
            gaT[1] * H +
            gaT[2] * u * H +
            gaT[3] * v * H +
            0.5 * gaT[4] * ((u^2 + v^2) * H + B)
        ) +
        Mt[4] * u * HL * δ -
        (Mt[5] + tau * Mt[4]) *
        u^2 *
        (
            faL[1] * HL +
            faL[2] * u * HL +
            faL[3] * v * HL +
            0.5 * faL[4] * ((u^2 + v^2) * HL + BL)
        ) *
        δ + Mt[4] * u * HR * (1.0 - δ) -
        (Mt[5] + tau * Mt[4]) *
        u^2 *
        (
            faR[1] * HR +
            faR[2] * u * HR +
            faR[3] * v * HR +
            0.5 * faR[4] * ((u^2 + v^2) * HR + BR)
        ) *
        (1.0 - δ) -
        tau *
        Mt[4] *
        u *
        (
            faTL[1] * HL +
            faTL[2] * u * HL +
            faTL[3] * v * HL +
            0.5 * faTL[4] * ((u^2 + v^2) * HL + BL)
        ) *
        δ -
        tau *
        Mt[4] *
        u *
        (
            faTR[1] * HR +
            faTR[2] * u * HR +
            faTR[3] * v * HR +
            0.5 * faTR[4] * ((u^2 + v^2) * HR + BR)
        ) *
        (1.0 - δ)
    @. fb =
        Mt[1] * u * B +
        Mt[2] *
        u^2 *
        (
            gaL[1] * B +
            gaL[2] * u * B +
            gaL[3] * v * B +
            0.5 * gaL[4] * ((u^2 + v^2) * B + Mxi[2] * H)
        ) *
        δ +
        Mt[2] *
        u^2 *
        (
            gaR[1] * B +
            gaR[2] * u * B +
            gaR[3] * v * B +
            0.5 * gaR[4] * ((u^2 + v^2) * B + Mxi[2] * H)
        ) *
        (1.0 - δ) +
        Mt[3] *
        u *
        (
            gaT[1] * B +
            gaT[2] * u * B +
            gaT[3] * v * B +
            0.5 * gaT[4] * ((u^2 + v^2) * B + Mxi[2] * H)
        ) +
        Mt[4] * u * BL * δ -
        (Mt[5] + tau * Mt[4]) *
        u^2 *
        (
            faL[1] * BL +
            faL[2] * u * BL +
            faL[3] * v * BL +
            0.5 * faL[4] * ((u^2 + v^2) * BL + Mxi[2] * HL)
        ) *
        δ + Mt[4] * u * BR * (1.0 - δ) -
        (Mt[5] + tau * Mt[4]) *
        u^2 *
        (
            faR[1] * BR +
            faR[2] * u * BR +
            faR[3] * v * BR +
            0.5 * faR[4] * ((u^2 + v^2) * BR + Mxi[2] * HR)
        ) *
        (1.0 - δ) -
        tau *
        Mt[4] *
        u *
        (
            faTL[1] * BL +
            faTL[2] * u * BL +
            faTL[3] * v * BL +
            0.5 * faTL[4] * ((u^2 + v^2) * BL + Mxi[2] * HL)
        ) *
        δ -
        tau *
        Mt[4] *
        u *
        (
            faTR[1] * BR +
            faTR[2] * u * BR +
            faTR[3] * v * BR +
            0.5 * faTR[4] * ((u^2 + v^2) * BR + Mxi[2] * HR)
        ) *
        (1.0 - δ)

    fh .*= dy
    fb .*= dy

    return nothing

end

"""
$(SIGNATURES)

Flux reconstruction
"""
function flux_gks!(
    fw::AV,
    wL::Y,
    wR::Y,
    inK,
    γ,
    μᵣ,
    ω,
    dt,
    swL = zero(wL)::Y,
    swR = zero(wR)::Y,
) where {Y<:AV}
    primL = conserve_prim(wL, γ)
    primR = conserve_prim(wR, γ)

    Mu1, Mv1, Mxi1, MuL1, MuR1 = gauss_moments(primL, inK)
    Mu2, Mv2, Mxi2, MuL2, MuR2 = gauss_moments(primR, inK)

    w =
        primL[1] .* moments_conserve(MuL1, Mv1, Mxi1, 0, 0, 0) .+
        primR[1] .* moments_conserve(MuR2, Mv2, Mxi2, 0, 0, 0)
    prim = conserve_prim(w, γ)
    tau =
        vhs_collision_time(prim, μᵣ, ω) +
        2.0 * dt * abs(primL[1] / primL[end] - primR[1] / primR[end]) /
        (primL[1] / primL[end] + primR[1] / primR[end])

    faL = pdf_slope(primL, swL, inK)
    sw = -primL[1] .* moments_conserve_slope(faL, Mu1, Mv1, Mxi1, 1, 0)
    faTL = pdf_slope(primL, sw, inK)

    faR = pdf_slope(primR, swR, inK)
    sw = -primR[1] .* moments_conserve_slope(faR, Mu2, Mv1, Mxi2, 1, 0)
    faTR = pdf_slope(primR, sw, inK)

    MuvL = moments_conserve(MuL1, Mv1, Mxi1, 1, 0, 0)
    MauL = moments_conserve_slope(faL, MuL1, Mv1, Mxi1, 2, 0)
    MauLT = moments_conserve_slope(faTL, MuL1, Mv1, Mxi1, 1, 0)

    MuvR = moments_conserve(MuR2, Mv2, Mxi2, 1, 0, 0)
    MauR = moments_conserve_slope(faR, MuR2, Mv2, Mxi2, 2, 0)
    MauRT = moments_conserve_slope(faTR, MuR2, Mv2, Mxi2, 1, 0)

    @. fw =
        primL[1] * (MuvL - tau * MauL - tau * MauLT) +
        primR[1] * (MuvR - tau * MauR - tau * MauRT)

    return nothing
end
