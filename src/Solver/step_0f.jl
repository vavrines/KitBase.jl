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
function step!(w::Y, prim::Y, fwL::X, fwR::X, γ, dx, RES, AVG) where {X<:AV,Y<:AV}
    w_old = deepcopy(w)

    @. w += (fwL - fwR) / dx
    prim .= conserve_prim(w, γ)

    @. RES += (w - w_old)^2
    @. AVG += abs(w)

    return nothing
end

"""
$(SIGNATURES)

1D0F with Workspace
"""
function step!(w::Y, prim::Y, fwL::X, fwR::X, γ, dx, RES, AVG, ws::StepWorkspace) where {X<:AV,Y<:AV}
    w_old = ws.w_old

    w_old .= w

    @. w += (fwL - fwR) / dx
    conserve_prim!(prim, w, γ)

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

1D0F Mixture with workspace
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
    ws::StepWorkspaceM,
) where {T1<:AM,T2<:AM}
    w_old, prim_old, mw, mprim, tau = 
        ws.w_old, ws.prim_old, ws.mw, ws.mprim, ws.tau

    # W^n
    w_old .= w
    prim_old .= prim

    # flux -> W^{n+1}
    @. w += (fwL - fwR) / dx
    mixture_conserve_prim!(prim, w, γ)

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
    aap_hs_collision_time!(tau, prim, mi, ni, me, ne, Kn)
    aap_hs_prim!(mprim, prim, tau, mi, ni, me, ne, Kn)
    mixture_prim_conserve!(mw, mprim, γ)
    for k in axes(w, 2)
        @. w[:, k] += (mw[:, k] - w_old[:, k]) * dt / tau[k]
    end
    mixture_conserve_prim!(prim, w, γ)

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
    dirc::AV,
    RES,
    AVG,
) where {T1<:AV}

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

2D0F with Workspace @ triangle
"""
function step!(
    w::T1,
    prim::T1,
    fw1::T1,
    fw2::T1,
    fw3::T1,
    γ,
    Δs,
    dirc::AV,
    RES,
    AVG,
    ws::StepWorkspace,
) where {T1<:AV}
    w_old = ws.w_old

    #--- store W^n and calculate shakhov term ---#
    w_old = deepcopy(w)

    #--- update W^{n+1} ---#
    @. w -= (fw1 * dirc[1] + fw2 * dirc[2] + fw3 * dirc[3]) / Δs
    conserve_prim!(prim, w, γ)

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
) where {T1<:AV}

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

2D0F with Workspace @ quadrilateral
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
    ws::StepWorkspace,
    collision,
) where {T1<:AV}
    w_old = ws.w_old

    #--- store W^n and calculate shakhov term ---#
    w_old .= w

    #--- update W^{n+1} ---#
    @. w += (fwL - fwR + fwD - fwU) / Δs
    conserve_prim!(prim, w, γ)

    #--- record residuals ---#
    @. RES += (w - w_old)^2
    @. AVG += abs(w)
end

"""
$(SIGNATURES)

3D0F
"""
function step!(
    w::T1,
    prim::T1,
    fwXL::T1,
    fwXR::T1,
    fwYL::T1,
    fwYR::T1,
    fwZL::T1,
    fwZR::T1,
    γ,
    Δv,
    RES,
    AVG,
    collision,
) where {T1<:AV}

    #--- store W^n and calculate shakhov term ---#
    w_old = deepcopy(w)

    #--- update W^{n+1} ---#
    @. w += (fwXL - fwXR + fwYL - fwYR + fwZL - fwZR) / Δv
    prim = conserve_prim(w, γ)

    #--- record residuals ---#
    @. RES += (w - w_old)^2
    @. AVG += abs(w)
end

"""
$(SIGNATURES)

3D0F with Workspace
"""
function step!(
    w::T1,
    prim::T1,
    fwXL::T1,
    fwXR::T1,
    fwYL::T1,
    fwYR::T1,
    fwZL::T1,
    fwZR::T1,
    γ,
    Δv,
    RES,
    AVG,
    ws::StepWorkspace,
    collision,
) where {T1<:AV}

    #--- store W^n and calculate shakhov term ---#
    w_old = ws.w_old
    w_old .= w

    #--- update W^{n+1} ---#
    @. w += (fwXL - fwXR + fwYL - fwYR + fwZL - fwZR) / Δv
    conserve_prim!(prim, w, γ)

    #--- record residuals ---#
    @. RES += (w - w_old)^2
    @. AVG += abs(w)
end