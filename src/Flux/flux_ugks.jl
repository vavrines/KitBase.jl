"""
$(SIGNATURES)

Unified gas kinetic scheme (UGKS)

1D2F1V
"""
function flux_ugks!(
    fw::AV,
    fh::T2,
    fb::T2,
    wL::T3,
    hL::T4,
    bL::T4,
    wR::T3,
    hR::T4,
    bR::T4,
    u::T5,
    ω::T5,
    inK,
    γ,
    visRef,
    visIdx,
    pr,
    dt,
    dxL,
    dxR,
    shL = zeros(eltype(hL), axes(hL))::T4,
    sbL = zeros(eltype(bL), axes(bL))::T4,
    shR = zeros(eltype(hR), axes(hR))::T4,
    sbR = zeros(eltype(bR), axes(bR))::T4,
) where {T2<:AV,T3<:AV,T4<:AV,T5<:AV} # 1D2F flux

    #--- reconstruct initial distribution ---#
    δ = heaviside.(u)

    h = @. hL * δ + hR * (1.0 - δ)
    b = @. bL * δ + bR * (1.0 - δ)
    sh = @. shL * δ + shR * (1.0 - δ)
    sb = @. sbL * δ + sbR * (1.0 - δ)

    #--- construct interface variables ---#
    #w = moments_conserve(h, b, u, v, ω)
    #prim = conserve_prim(w, γ)

    primL = conserve_prim(wL, γ)
    primR = conserve_prim(wR, γ)
    Mu1, Mxi1, MuL1, MuR1 = gauss_moments(primL, inK)
    Muv1 = moments_conserve(MuL1, Mxi1, 0, 0)
    Mu2, Mxi2, MuL2, MuR2 = gauss_moments(primR, inK)
    Muv2 = moments_conserve(MuR2, Mxi2, 0, 0)

    w = @. primL[1] * Muv1 + primR[1] * Muv2
    prim = conserve_prim(w, γ)
    aL = pdf_slope(prim, (w .- wL) ./ dxL, inK)
    aR = pdf_slope(prim, (wR .- w) ./ dxR, inK)

    Mu, Mxi, MuL, MuR = gauss_moments(prim, inK)
    MauL = moments_conserve_slope(aL, MuL, Mxi, 1)
    MauR = moments_conserve_slope(aR, MuR, Mxi, 1)
    aT = pdf_slope(prim, -prim[1] .* (MauL .+ MauR), inK)

    #--- calculate integral time constants ---#
    tau = vhs_collision_time(prim, visRef, visIdx)

    Mt = zeros(5)
    Mt[4] = tau * (1.0 - exp(-dt / tau)) # f0
    Mt[5] = -tau * dt * exp(-dt / tau) + tau * Mt[4]
    Mt[1] = dt - Mt[4] # M0
    Mt[2] = -tau * Mt[1] + Mt[5]
    Mt[3] = dt^2 / 2.0 - tau * Mt[1]

    #--- calculate flux from M0 ---#
    Muv = moments_conserve(Mu, Mxi, 1, 0)
    MauL = moments_conserve_slope(aL, MuL, Mxi, 2)
    MauR = moments_conserve_slope(aR, MuR, Mxi, 2)
    MauT = moments_conserve_slope(aT, Mu, Mxi, 1)

    @. fw = Mt[1] * prim[1] * Muv + Mt[2] * prim[1] * (MauL + MauR) + Mt[3] * prim[1] * MauT

    #--- calculate flux from f0 ---#
    H = maxwellian(u, prim)
    B = H .* inK ./ (2.0 * prim[end])

    fw[1] += Mt[4] * sum(ω .* u .* h) - Mt[5] * sum(ω .* u .^ 2 .* sh)
    fw[2] += Mt[4] * sum(ω .* u .^ 2 .* h) - Mt[5] * sum(ω .* u .^ 3 .* sh)
    fw[3] +=
        Mt[4] * 0.5 * (sum(ω .* u .^ 3 .* h) + sum(ω .* u .* b)) -
        Mt[5] * 0.5 * (sum(ω .* u .^ 4 .* sh) + sum(ω .* u .^ 2 .* sb))

    @. fh =
        Mt[1] * u * H +
        Mt[2] * u^2 * (aL[1] * H + aL[2] * u * H + 0.5 * aL[3] * (u^2 * H + B)) * δ +
        Mt[2] *
        u^2 *
        (aR[1] * H + aR[2] * u * H + 0.5 * aR[3] * (u^2 * H + B)) *
        (1.0 - δ) +
        Mt[3] * u * (aT[1] * H + aT[2] * u * H + 0.5 * aT[3] * (u^2 * H + B)) +
        Mt[4] * u * h - Mt[5] * u^2 * sh
    @. fb =
        Mt[1] * u * B +
        Mt[2] *
        u^2 *
        (aL[1] * B + aL[2] * u * B + 0.5 * aL[3] * (u^2 * B + Mxi[2] * H)) *
        δ +
        Mt[2] *
        u^2 *
        (aR[1] * B + aR[2] * u * B + 0.5 * aR[3] * (u^2 * B + Mxi[2] * H)) *
        (1.0 - δ) +
        Mt[3] * u * (aT[1] * B + aT[2] * u * B + 0.5 * aT[3] * (u^2 * B + Mxi[2] * H)) +
        Mt[4] * u * b - Mt[5] * u^2 * sb

    return nothing

end

"""
$(SIGNATURES)

2D2F2V
"""
function flux_ugks!(
    fw::AV,
    fh::T2,
    fb::T2,
    wL::T3,
    hL::T4,
    bL::T4,
    wR::T3,
    hR::T4,
    bR::T4,
    u::T5,
    v::T5,
    ω::T5,
    inK,
    γ,
    visRef,
    visIdx,
    pr,
    dt,
    dxL,
    dxR,
    len,
    shL = zeros(eltype(hL), axes(hL))::T4,
    sbL = zeros(eltype(bL), axes(bL))::T4,
    shR = zeros(eltype(hR), axes(hR))::T4,
    sbR = zeros(eltype(bR), axes(bR))::T4,
) where {T2<:AM,T3<:AV,T4<:AM,T5<:AM} # 2D2F flux

    #--- reconstruct initial distribution ---#
    δ = heaviside.(u)

    h = @. hL * δ + hR * (1.0 - δ)
    b = @. bL * δ + bR * (1.0 - δ)
    sh = @. shL * δ + shR * (1.0 - δ)
    sb = @. sbL * δ + sbR * (1.0 - δ)

    #--- construct interface variables ---#
    #w = moments_conserve(h, b, u, v, ω)
    #prim = conserve_prim(w, γ)

    primL = conserve_prim(wL, γ)
    primR = conserve_prim(wR, γ)
    Mu1, Mv1, Mxi1, MuL1, MuR1 = gauss_moments(primL, inK)
    Muv1 = moments_conserve(MuL1, Mv1, Mxi1, 0, 0, 0)
    Mu2, Mv2, Mxi2, MuL2, MuR2 = gauss_moments(primR, inK)
    Muv2 = moments_conserve(MuR2, Mv2, Mxi2, 0, 0, 0)
    w = @. primL[1] * Muv1 + primR[1] * Muv2
    prim = conserve_prim(w, γ)

    aL = pdf_slope(prim, (w .- wL) ./ dxL, inK)
    aR = pdf_slope(prim, (wR .- w) ./ dxR, inK)

    Mu, Mv, Mxi, MuL, MuR = gauss_moments(prim, inK)
    MauL = moments_conserve_slope(aL, MuL, Mv, Mxi, 1, 0)
    MauR = moments_conserve_slope(aR, MuR, Mv, Mxi, 1, 0)
    aT = pdf_slope(prim, -prim[1] .* (MauL .+ MauR), inK)

    #--- calculate integral time constants ---#
    tau = vhs_collision_time(prim, visRef, visIdx)

    Mt = zeros(5)
    Mt[4] = tau * (1.0 - exp(-dt / tau)) # f0
    Mt[5] = -tau * dt * exp(-dt / tau) + tau * Mt[4]
    Mt[1] = dt - Mt[4] # M0
    Mt[2] = -tau * Mt[1] + Mt[5]
    Mt[3] = dt^2 / 2.0 - tau * Mt[1]

    # --- calculate flux from M0 ---#
    Muv = moments_conserve(Mu, Mv, Mxi, 1, 0, 0)
    MauL = moments_conserve_slope(aL, MuL, Mv, Mxi, 2, 0)
    MauR = moments_conserve_slope(aR, MuR, Mv, Mxi, 2, 0)
    MauT = moments_conserve_slope(aT, Mu, Mv, Mxi, 1, 0)

    @. fw = Mt[1] * prim[1] * Muv + Mt[2] * prim[1] * (MauL + MauR) + Mt[3] * prim[1] * MauT

    # --- calculate flux from f0 ---#
    H = maxwellian(u, v, prim)
    B = H .* inK ./ (2.0 * prim[end])

    fw[1] += Mt[4] * sum(ω .* u .* h) - Mt[5] * sum(ω .* u .^ 2 .* sh)
    fw[2] += Mt[4] * sum(ω .* u .^ 2 .* h) - Mt[5] * sum(ω .* u .^ 3 .* sh)
    fw[3] += Mt[4] * sum(ω .* v .* u .* h) - Mt[5] * sum(ω .* v .* u .^ 2 .* sh)
    fw[4] +=
        Mt[4] * 0.5 * (sum(ω .* u .* (u .^ 2 .+ v .^ 2) .* h) + sum(ω .* u .* b)) -
        Mt[5] *
        0.5 *
        (sum(ω .* u .^ 2 .* (u .^ 2 .+ v .^ 2) .* sh) + sum(ω .* u .^ 2 .* sb))

    @. fh =
        Mt[1] * u * H +
        Mt[2] *
        u^2 *
        (aL[1] * H + aL[2] * u * H + aL[3] * v * H + 0.5 * aL[4] * ((u^2 + v^2) * H + B)) *
        δ +
        Mt[2] *
        u^2 *
        (aR[1] * H + aR[2] * u * H + aR[3] * v * H + 0.5 * aR[4] * ((u^2 + v^2) * H + B)) *
        (1.0 - δ) +
        Mt[3] *
        u *
        (aT[1] * H + aT[2] * u * H + aT[3] * v * H + 0.5 * aT[4] * ((u^2 + v^2) * H + B)) +
        Mt[4] * u * h - Mt[5] * u^2 * sh
    @. fb =
        Mt[1] * u * B +
        Mt[2] *
        u^2 *
        (
            aL[1] * B +
            aL[2] * u * B +
            aL[3] * v * B +
            0.5 * aL[4] * ((u^2 + v^2) * B + Mxi[2] * H)
        ) *
        δ +
        Mt[2] *
        u^2 *
        (
            aR[1] * B +
            aR[2] * u * B +
            aR[3] * v * B +
            0.5 * aR[4] * ((u^2 + v^2) * B + Mxi[2] * H)
        ) *
        (1.0 - δ) +
        Mt[3] *
        u *
        (
            aT[1] * B +
            aT[2] * u * B +
            aT[3] * v * B +
            0.5 * aT[4] * ((u^2 + v^2) * B + Mxi[2] * H)
        ) +
        Mt[4] * u * b - Mt[5] * u^2 * sb

    # multiply interface length
    fw .*= len
    fh .*= len
    fb .*= len

    return nothing

end

"""
$(SIGNATURES)

3F2V with AAP model
"""
function flux_ugks!(
    fw::AM,
    fh0::T2,
    fh1::T2,
    fh2::T2,
    wL::T3,
    h0L::T4,
    h1L::T4,
    h2L::T4,
    wR::T3,
    h0R::T4,
    h1R::T4,
    h2R::T4,
    u::T5,
    v::T5,
    ω::T5,
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
    sh0L = zeros(eltype(h0L), axes(h0L))::T4,
    sh1L = zeros(eltype(h1L), axes(h1L))::T4,
    sh2L = zeros(eltype(h2L), axes(h2L))::T4,
    sh0R = zeros(eltype(h0R), axes(h0R))::T4,
    sh1R = zeros(eltype(h1R), axes(h1R))::T4,
    sh2R = zeros(eltype(h2R), axes(h2R))::T4,
) where {T2<:AA3,T3<:AM,T4<:AA3,T5<:AA3}

    #--- reconstruct initial distribution ---#
    δ = heaviside.(u)

    h0 = @. h0L * δ + h0R * (1.0 - δ)
    h1 = @. h1L * δ + h1R * (1.0 - δ)
    h2 = @. h2L * δ + h2R * (1.0 - δ)

    sh0 = @. sh0L * δ + sh0R * (1.0 - δ)
    sh1 = @. sh1L * δ + sh1R * (1.0 - δ)
    sh2 = @. sh2L * δ + sh2R * (1.0 - δ)

    primL = mixture_conserve_prim(wL, γ)
    primR = mixture_conserve_prim(wR, γ)

    #--- construct interface distribution ---#
    Mu1, Mv1, Mxi1, MuL1, MuR1 = mixture_gauss_moments(primL, inK)
    Muv1 = mixture_moments_conserve(MuL1, Mv1, Mxi1, 0, 0, 0)
    Mu2, Mv2, Mxi2, MuL2, MuR2 = mixture_gauss_moments(primR, inK)
    Muv2 = mixture_moments_conserve(MuR2, Mv2, Mxi2, 0, 0, 0)

    w = similar(wL)
    for j in axes(w, 2)
        @. w[:, j] = primL[1, j] * Muv1[:, j] + primR[1, j] * Muv2[:, j]
    end
    prim = mixture_conserve_prim(w, γ)
    tau = aap_hs_collision_time(prim, mi, ni, me, ne, Kn)
    prim = aap_hs_prim(prim, tau, mi, ni, me, ne, Kn)

    a = mixture_pdf_slope(prim, (wR .- wL) ./ (dxR + dxL), inK)
    Mu, Mv, Mxi, MuL, MuR = mixture_gauss_moments(prim, inK)
    Mau = mixture_moments_conserve_slope(a, Mu, Mv, Mxi, 1, 0, 0)
    aT = mixture_pdf_slope(prim, -prim[1] .* Mau, inK)

    Mt = zeros(2, 2)
    @. Mt[2, :] = tau * (1.0 - exp(-dt / tau)) # f0
    @. Mt[1, :] = dt - Mt[2, :] # M0

    Mt = zeros(5, axes(w, 2))
    for j in axes(Mt, 2)
        Mt[4, j] = tau[j] * (1.0 - exp(-dt / tau[j])) # f0
        Mt[5, j] = -tau[j] * dt * exp(-dt / tau[j]) + tau[j] * Mt[4, j]
        Mt[1, j] = dt - Mt[4, j] # M0
        Mt[2, j] = -tau[j] * Mt[1, j] + Mt[5, j]
        Mt[3, j] = dt^2 / 2.0 - tau[j] * Mt[1, j]
    end

    #--- calculate interface flux ---#
    # flux from M0
    Muv = mixture_moments_conserve(Mu, Mv, Mxi, 1, 0, 0)
    Mau = mixture_moments_conserve_slope(a, Mu, Mv, Mxi, 2, 0, 0)
    MauT = mixture_moments_conserve_slope(aT, Mu, Mv, Mxi, 1, 0, 0)
    for j in axes(fw, 2)
        @. fw[:, j] =
            Mt[1, j] * prim[1, j] * Muv[:, j] +
            Mt[2, j] * prim[1, j] * Mau[:, j] +
            Mt[3, j] * prim[1, j] * MauT[:, j]
    end

    # flux from f0
    H0 = mixture_maxwellian(u, v, prim)
    H1 = similar(H0)
    H2 = similar(H0)
    for j in axes(H0, 3)
        H1[:, :, j] = H0[:, :, j] .* prim[4, j]
        H2[:, :, j] .= H0[:, :, j] .* (prim[4, j]^2 + 1.0 / (2.0 * prim[5, j]))
    end

    for j in axes(fw, 2)
        fw[1, j] +=
            Mt[4, j] * sum(ω[:, :, j] .* u[:, :, j] .* h0[:, :, j]) -
            Mt[5, j] * sum(ω[:, :, j] .* u[:, :, j] .^ 2 .* sh0[:, :, j])
        fw[2, j] +=
            Mt[4, j] * sum(ω[:, :, j] .* u[:, :, j] .^ 2 .* h0[:, :, j]) -
            Mt[5, j] * sum(ω[:, :, j] .* u[:, :, j] .^ 3 .* sh0[:, :, j])
        fw[3, j] +=
            Mt[4, j] * sum(ω[:, :, j] .* v[:, :, j] .* u[:, :, j] .* h0[:, :, j]) -
            Mt[5, j] * sum(ω[:, :, j] .* v[:, :, j] .* u[:, :, j] .^ 2 .* sh0[:, :, j])
        fw[4, j] +=
            Mt[4, j] * sum(ω[:, :, j] .* u[:, :, j] .* h1[:, :, j]) -
            Mt[5, j] * sum(ω[:, :, j] .* u[:, :, j] .^ 2 .* sh1[:, :, j])
        fw[5, j] +=
            Mt[4, j] *
            0.5 *
            (
                sum(
                    ω[:, :, j] .* u[:, :, j] .* (u[:, :, j] .^ 2 .+ v[:, :, j] .^ 2) .*
                    h0[:, :, j],
                ) + sum(ω[:, :, j] .* u[:, :, j] .* h2[:, :, j])
            ) -
            Mt[5, j] *
            0.5 *
            (
                sum(
                    ω[:, :, j] .* u[:, :, j] .^ 2 .* (u[:, :, j] .^ 2 .+ v[:, :, j] .^ 2) .*
                    sh0[:, :, j],
                ) + sum(ω[:, :, j] .* u[:, :, j] .^ 2 .* sh2[:, :, j])
            )

        @. fh0[:, :, j] =
            Mt[1, j] * u[:, :, j] * H0[:, :, j] +
            Mt[2, j] *
            u[:, :, j]^2 *
            (
                a[1, j] * H0[:, :, j] +
                a[2, j] * u[:, :, j] * H0[:, :, j] +
                a[3, j] * v[:, :, j] * H0[:, :, j] +
                a[4, j] * u[:, :, j] * H1[:, :, j] +
                0.5 * a[5, j] * ((u[:, :, j]^2 + v[:, :, j]^2) * H0[:, :, j] + H2[:, :, j])
            ) +
            Mt[3, j] *
            u[:, :, j] *
            (
                aT[1, j] * H0[:, :, j] +
                aT[2, j] * u[:, :, j] * H0[:, :, j] +
                aT[3, j] * v[:, :, j] * H0[:, :, j] +
                aT[4, j] * u[:, :, j] * H1[:, :, j] +
                0.5 * aT[5, j] * ((u[:, :, j]^2 + v[:, :, j]^2) * H0[:, :, j] + H2[:, :, j])
            ) +
            Mt[4, j] * u[:, :, j] * h0[:, :, j] - Mt[5, j] * u[:, :, j]^2 * sh0[:, :, j]
        @. fh1[:, :, j] =
            Mt[1, j] * u[:, :, j] * H1[:, :, j] +
            Mt[2, j] *
            u[:, :, j]^2 *
            (
                a[1, j] * H1[:, :, j] +
                a[2, j] * u[:, :, j] * H1[:, :, j] +
                a[3, j] * v[:, :, j] * H1[:, :, j] +
                a[4, j] * u[:, :, j] * H2[:, :, j] +
                0.5 *
                a[5, j] *
                ((u[:, :, j]^2 + v[:, :, j]^2) * H1[:, :, j] + Mxi[3, j] * H0[:, :, j])
            ) +
            Mt[3, j] *
            u[:, :, j] *
            (
                aT[1, j] * H1[:, :, j] +
                aT[2, j] * u[:, :, j] * H1[:, :, j] +
                aT[3, j] * v[:, :, j] * H1[:, :, j] +
                aT[4, j] * u[:, :, j] * H2[:, :, j] +
                0.5 *
                aT[5, j] *
                ((u[:, :, j]^2 + v[:, :, j]^2) * H1[:, :, j] + Mxi[3, j] * H0[:, :, j])
            ) +
            Mt[4, j] * u[:, :, j] * h1[:, :, j] - Mt[5, j] * u[:, :, j]^2 * sh1[:, :, j]
        @. fh2[:, :, j] =
            Mt[1, j] * u[:, :, j] * H2[:, :, j] +
            Mt[2, j] *
            u[:, :, j]^2 *
            (
                a[1, j] * H2[:, :, j] +
                a[2, j] * u[:, :, j] * H2[:, :, j] +
                a[3, j] * v[:, :, j] * H2[:, :, j] +
                a[4, j] * u[:, :, j] * Mxi[3, j] * H0[:, :, j] +
                0.5 *
                a[5, j] *
                ((u[:, :, j]^2 + v[:, :, j]^2) * H2[:, :, j] + Mxi[4, j] * H0[:, :, j])
            ) +
            Mt[3, j] *
            u[:, :, j] *
            (
                aT[1, j] * H2[:, :, j] +
                aT[2, j] * u[:, :, j] * H2[:, :, j] +
                aT[3, j] * v[:, :, j] * H2[:, :, j] +
                aT[4, j] * u[:, :, j] * Mxi[3, j] * H0[:, :, j] +
                0.5 *
                aT[5, j] *
                ((u[:, :, j]^2 + v[:, :, j]^2) * H2[:, :, j] + Mxi[4, j] * H0[:, :, j])
            ) +
            Mt[4, j] * u[:, :, j] * h2[:, :, j] - Mt[5, j] * u[:, :, j]^2 * sh2[:, :, j]
    end

    @. fw *= len
    @. fh0 *= len
    @. fh1 *= len
    @. fh2 *= len

    return nothing

end
