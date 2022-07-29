"""
$(SIGNATURES)

Scalar
"""
function step!(w::Number, prim, fwL, fwR, a, dx, RES, AVG)
    # store W^n
    w_old = deepcopy(w)

    # update W^{n+1}
    w += (fwL - fwR) / dx
    prim .= conserve_prim(w, a)

    # record residuals
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
) where {X<:AV,Y<:AV}

    w_old = deepcopy(w)

    @. w += (fwL - fwR) / dx
    prim .= conserve_prim(w, γ)

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
) where {T1<:AM,T2<:AM}

    # W^n
    w_old = deepcopy(w)
    prim_old = deepcopy(prim)

    # flux -> W^{n+1}
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

    # source -> W^{n+1}
    tau = aap_hs_collision_time(prim, mi, ni, me, ne, Kn)
    mprim = aap_hs_prim(prim, tau, mi, ni, me, ne, Kn)
    mw = mixture_prim_conserve(mprim, γ)
    for k in axes(w, 2)
        @. w[:, k] += (mw[:, k] - w_old[:, k]) * dt / tau[k]
    end
    prim .= mixture_conserve_prim(w, γ)

    @. RES += (w_old - w)^2
    @. AVG += abs(w)

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
    collision,
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
