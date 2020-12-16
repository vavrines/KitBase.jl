"""
Equilibrium part of gas kinetic flux

works with fluid-kinetic simulations

"""
function flux_equilibrium!(
    fw::T1,
    wL::T2,
    wR::T2,
    inK::Real,
    γ::Real,
    μᵣ::Real,
    ω::Real,
    Pr::Real,
    dt::Real,
    dxL::Real,
    dxR::Real,
    swL = zeros(eltype(fw), axes(wL))::T1,
    swR = zeros(eltype(fw), axes(wR))::T1,
) where {T1<:AbstractArray{<:AbstractFloat,1},T2<:AbstractArray{<:Real,1}} # 1D

    primL = conserve_prim(wL, γ)
    primR = conserve_prim(wR, γ)

    Mu1, Mxi1, MuL1, MuR1 = gauss_moments(primL, inK)
    Mu2, Mxi2, MuL2, MuR2 = gauss_moments(primR, inK)

    w =
        primL[1] .* moments_conserve(MuL1, Mxi1, 0, 0) .+
        primR[1] .* moments_conserve(MuR2, Mxi2, 0, 0)
    prim = conserve_prim(w, γ)

    # positivity check
    if prim[end] < 0.0
        primL = conserve_prim(wL .- swL .* dxL, γ)
        primR = conserve_prim(wR .+ swR .* dxR, γ)

        Mu1, Mxi1, MuL1, MuR1 = gauss_moments(primL, inK)
        Mu2, Mxi2, MuL2, MuR2 = gauss_moments(primR, inK)

        w =
            primL[1] .* moments_conserve(MuL1, Mxi1, 0, 0) .+
            primR[1] .* moments_conserve(MuR2, Mxi2, 0, 0)
        prim = conserve_prim(w, γ)

        tau =
            vhs_collision_time(prim, μᵣ, ω) #+
            #2.0 * dt * abs(primL[1] / primL[end] - primR[1] / primR[end]) /
            #(primL[1] / primL[end] + primR[1] / primR[end])

        # time-integration constants
        Mt = zeros(2)
        Mt[2] = tau * (1.0 - exp(-dt / tau))
        Mt[1] = dt - Mt[2]

        # flux related to central distribution
        Mu, Mxi, MuL, MuR = gauss_moments(prim, inK)
        Muv = moments_conserve(Mu, Mxi, 1, 0)
        fw .= Mt[1] .* prim[1] .* Muv

    else
        tau =
            vhs_collision_time(prim, μᵣ, ω) +
            2.0 * dt * abs(primL[1] / primL[end] - primR[1] / primR[end]) /
            (primL[1] / primL[end] + primR[1] / primR[end])

        Mu, Mxi, MuL, MuR = gauss_moments(prim, inK)
        sw0L = (w .- wL) ./ dxL
        sw0R = (wR .- w) ./ dxR
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
    end
    
    return nothing

end

#--- 2D ---#
function flux_equilibrium!(
    fw::T1,
    wL::T2,
    wR::T2,
    inK::Real,
    γ::Real,
    μᵣ::Real,
    ω::Real,
    Pr::Real,
    dt::Real,
    dxL::Real,
    dxR::Real,
    dy::Real,
    swL = zeros(eltype(fw), axes(wL))::T1,
    swR = zeros(eltype(fw), axes(wR))::T1,
) where {T1<:AbstractArray{<:AbstractFloat,1},T2<:AbstractArray{<:Real,1}}

    Mu1, Mv1, Mxi1, MuL1, MuR1 = gauss_moments(primL, inK)
    Mu2, Mv2, Mxi2, MuL2, MuR2 = gauss_moments(primR, inK)

    w =
        primL[1] .* moments_conserve(MuL1, Mv1, Mxi1, 0, 0, 0) .+
        primR[1] .* moments_conserve(MuR2, Mv2, Mxi2, 0, 0, 0)
    prim = conserve_prim(w, γ)

    # positivity check
    if prim[end] < 0.0
        primL = conserve_prim(wL .- swL .* dxL, γ)
        primR = conserve_prim(wR .+ swR .* dxR, γ)

        Mu1, Mv1, Mxi1, MuL1, MuR1 = gauss_moments(primL, inK)
        Mu2, Mv2, Mxi2, MuL2, MuR2 = gauss_moments(primR, inK)

        w =
            primL[1] .* moments_conserve(MuL1, Mv1, Mxi1, 0, 0, 0) .+
            primR[1] .* moments_conserve(MuR2, Mv2, Mxi2, 0, 0, 0)
        prim = conserve_prim(w, γ)

        swL .= 0.0
        swR .= 0.0
    end

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
    sw0L = (w .- wL) ./ dxL
    sw0R = (wR .- w) ./ dxR
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

    @. fw =
        (Mt[1] * prim[1] * Muv + Mt[2] * prim[1] * (MauL + MauR) + Mt[3] * prim[1] * MauT) *
        dy
    # fw .= (Mt[1] .* prim[1] .* Muv .+ Mt[2] .* prim[1] .* Mau .+ Mt[3] .* prim[1] .* MauT) .* dy

    return nothing

end
