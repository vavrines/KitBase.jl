"""
$(SIGNATURES)

Kinetic flux vector splitting (KFVS) flux
"""
function flux_kfvs!(
    face::Interface1F,
    ctrL::T,
    ctrR::T,
    vs::VSpace1D,
    p,
    dt = 1.0,
) where {T<:ControlVolume1F}

    dxL, dxR = p

    flux_kfvs!(
        face.fw,
        face.ff,
        ctrL.f + ctrL.sf * dxL,
        ctrR.f - ctrR.sf * dxR,
        vs.u,
        vs.weights,
        dt,
        ctrL.sf,
        ctrR.sf,
    )

    return nothing

end

function flux_kfvs!(
    face::Interface2F,
    ctrL::T,
    ctrR::T,
    vs::VSpace1D,
    p,
    dt = 1.0,
) where {T<:ControlVolume2F}

    dxL, dxR = p

    flux_kfvs!(
        face.fw,
        face.fh,
        face.fb,
        ctrL.h .+ ctrL.sb .* dxL,
        ctrL.b .+ ctrL.sb .* dxL,
        ctrR.h .- ctrR.sb .* dxR,
        ctrR.b .- ctrR.sb .* dxR,
        vs.u,
        vs.weights,
        dt,
        ctrL.sh,
        ctrL.sb,
        ctrR.sh,
        ctrR.sb,
    )

    return nothing

end

function flux_kfvs!(
    face::Interface1F,
    ctrL::T,
    ctrR::T,
    vs::VSpace2D,
    p,
    dt = 1.0,
) where {T<:ControlVolume1F}

    dxL, dxR, len = p

    flux_kfvs!(
        face.fw,
        face.ff,
        ctrL.f + ctrL.sf * dxL,
        ctrR.f - ctrR.sf * dxR,
        vs.u,
        vs.v,
        vs.weights,
        dt,
        len,
        ctrL.sf,
        ctrR.sf,
    )

    return nothing

end

function flux_kfvs!(
    face::Interface2F,
    ctrL::T,
    ctrR::T,
    vs::VSpace2D,
    p,
    dt = 1.0,
) where {T<:ControlVolume2F}

    dxL, dxR, len = p

    flux_kfvs!(
        face.fw,
        face.fh,
        face.fb,
        ctrL.h .+ ctrL.sb .* dxL,
        ctrL.b .+ ctrL.sb .* dxL,
        ctrR.h .- ctrR.sb .* dxR,
        ctrR.b .- ctrR.sb .* dxR,
        vs.u,
        vs.v,
        vs.weights,
        dt,
        len,
        ctrL.sh,
        ctrL.sb,
        ctrR.sh,
        ctrR.sb,
    )

    return nothing

end

# ------------------------------------------------------------
# Low-level backends
# ------------------------------------------------------------

"""
$(SIGNATURES)

Kinetic flux vector splitting (KFVS) flux
"""
function flux_kfvs(
    fL::Y,
    fR::Y,
    u::AV{<:FN},
    dt,
    sfL = zero(fL)::Y,
    sfR = zero(fR)::Y,
) where {Y<:AV{<:FN}}

    ff = similar(fL)
    flux_kfvs!(ff, fL, fR, u, dt, sfL, sfR)

    return ff

end

"""
$(SIGNATURES)

Kinetic flux vector splitting (KFVS) flux

- @args: particle distribution functions and their left/right slopes
- @args: particle velocity quadrature points and weights
- @args: time step and cell size

1F1V for pure DOM
"""
function flux_kfvs!(
    ff::AV{<:FN},
    fL::Y,
    fR::Y,
    u::AV{<:FN},
    dt,
    sfL = zero(fL)::Y,
    sfR = zero(fR)::Y,
) where {Y<:AV{<:FN}} # 1F1V flux for pure DOM

    # upwind reconstruction
    δ = heaviside.(u)

    f = @. fL * δ + fR * (1.0 - δ)
    sf = @. sfL * δ + sfR * (1.0 - δ)

    # calculate flux
    @. ff = dt * u * f - 0.5 * dt^2 * u^2 * sf

    return nothing

end

"""
$(SIGNATURES)

1F1V
"""
function flux_kfvs!(
    fw::AV{<:FN},
    ff::AV{<:FN},
    fL::Z,
    fR::Z,
    u::A,
    ω::A,
    dt,
    sfL = zero(fL)::Z,
    sfR = zero(fR)::Z,
) where {Z<:AV{<:FN},A<:AV{<:FN}}

    # upwind reconstruction
    δ = heaviside.(u)

    f = @. fL * δ + fR * (1.0 - δ)
    sf = @. sfL * δ + sfR * (1.0 - δ)

    # calculate fluxes
    fw[1] = dt * sum(ω .* u .* f) - 0.5 * dt^2 * sum(ω .* u .^ 2 .* sf)
    fw[2] = dt * sum(ω .* u .^ 2 .* f) - 0.5 * dt^2 * sum(ω .* u .^ 3 .* sf)
    fw[3] = dt * 0.5 * sum(ω .* u .^ 3 .* f) - 0.5 * dt^2 * 0.5 * sum(ω .* u .^ 4 .* sf)

    @. ff = dt * u * f - 0.5 * dt^2 * u^2 * sf

    return nothing

end

"""
$(SIGNATURES)

Mixture
"""
function flux_kfvs!(
    fw::AM{<:FN},
    ff::AM{<:FN},
    fL::Z,
    fR::Z,
    u::A,
    ω::A,
    dt,
    sfL = zero(fL)::Z,
    sfR = zero(fR)::Z,
) where {Z<:AM{<:FN},A<:AM{<:FN}}

    for j in axes(fw, 2)
        _fw = @view fw[:, j]
        _ff = @view ff[:, j]
        _fL = @view fL[:, j]
        _fR = @view fR[:, j]
        _u = @view u[:, j]
        _ω = @view ω[:, j]
        _sfL = @view sfL[:, j]
        _sfR = @view sfR[:, j]

        flux_kfvs!(_fw, _ff, _fL, _fR, _u, _ω, dt, _sfL, _sfR)
    end

    return nothing

end

"""
$(SIGNATURES)

2F1V
"""
function flux_kfvs!(
    fw::AV{<:FN},
    fh::Y,
    fb::Y,
    hL::Z,
    bL::Z,
    hR::Z,
    bR::Z,
    u::A,
    ω::A,
    dt,
    shL = zero(hL)::Z,
    sbL = zero(bL)::Z,
    shR = zero(hR)::Z,
    sbR = zero(bR)::Z,
) where {Y<:AV{<:FN},Z<:AV{<:FN},A<:AV{<:FN}}

    # upwind reconstruction
    δ = heaviside.(u)

    h = @. hL * δ + hR * (1.0 - δ)
    b = @. bL * δ + bR * (1.0 - δ)

    sh = @. shL * δ + shR * (1.0 - δ)
    sb = @. sbL * δ + sbR * (1.0 - δ)

    # calculate fluxes
    fw[1] = dt * sum(ω .* u .* h) - 0.5 * dt^2 * sum(ω .* u .^ 2 .* sh)
    fw[2] = dt * sum(ω .* u .^ 2 .* h) - 0.5 * dt^2 * sum(ω .* u .^ 3 .* sh)
    fw[3] =
        dt * 0.5 * (sum(ω .* u .^ 3 .* h) + sum(ω .* u .* b)) -
        0.5 * dt^2 * 0.5 * (sum(ω .* u .^ 4 .* sh) + sum(ω .* u .^ 2 .* sb))

    @. fh = dt * u * h - 0.5 * dt^2 * u^2 * sh
    @. fb = dt * u * b - 0.5 * dt^2 * u^2 * sb

    return nothing

end

"""
$(SIGNATURES)

Mixture
"""
function flux_kfvs!(
    fw::AM{<:FN},
    fh::Y,
    fb::Y,
    hL::Z,
    bL::Z,
    hR::Z,
    bR::Z,
    u::A,
    ω::A,
    dt,
    shL = zero(hL)::Z,
    sbL = zero(bL)::Z,
    shR = zero(hR)::Z,
    sbR = zero(bR)::Z,
) where {Y<:AM{<:FN},Z<:AM{<:FN},A<:AM{<:FN}}

    for j in axes(fw, 2)
        _fw = @view fw[:, j]
        _fh = @view fh[:, j]
        _fb = @view fb[:, j]
        _hL = @view hL[:, j]
        _bL = @view bL[:, j]
        _hR = @view hR[:, j]
        _bR = @view bR[:, j]
        _u = @view u[:, j]
        _ω = @view ω[:, j]
        _shL = @view shL[:, j]
        _sbL = @view sbL[:, j]
        _shR = @view shR[:, j]
        _sbR = @view sbR[:, j]

        flux_kfvs!(_fw, _fh, _fb, _hL, _bL, _hR, _bR, _u, _ω, dt, _shL, _sbL, _shR, _sbR)
    end

    return nothing

end

"""
$(SIGNATURES)

3F1V @ Rykov
"""
function flux_kfvs!(
    fw::AV{<:FN},
    fh::Y,
    fb::Y,
    fr::Y,
    hL::Z,
    bL::Z,
    rL::Z,
    hR::Z,
    bR::Z,
    rR::Z,
    u::A,
    ω::A,
    dt,
    shL = zero(hL)::Z,
    sbL = zero(bL)::Z,
    srL = zero(rL)::Z,
    shR = zero(hR)::Z,
    sbR = zero(bR)::Z,
    srR = zero(rR)::Z,
) where {Y<:AV{<:FN},Z<:AV{<:FN},A<:AV{<:FN}}

    # upwind reconstruction
    δ = heaviside.(u)

    h = @. hL * δ + hR * (1.0 - δ)
    b = @. bL * δ + bR * (1.0 - δ)
    r = @. rL * δ + rR * (1.0 - δ)

    sh = @. shL * δ + shR * (1.0 - δ)
    sb = @. sbL * δ + sbR * (1.0 - δ)
    sr = @. srL * δ + srR * (1.0 - δ)

    # macro fluxes
    fw[1] = dt * sum(ω .* u .* h) - 0.5 * dt^2 * sum(ω .* u .^ 2 .* sh)
    fw[2] = dt * sum(ω .* u .^ 2 .* h) - 0.5 * dt^2 * sum(ω .* u .^ 3 .* sh)
    fw[3] =
        dt * 0.5 * (sum(ω .* u .^ 3 .* h) + sum(ω .* u .* b)) -
        0.5 * dt^2 * 0.5 * (sum(ω .* u .^ 4 .* sh) + sum(ω .* u .^ 2 .* sb))
    fw[4] = dt * 0.5 * sum(ω .* u .* r) - 0.5 * dt^2 * sum(ω .* u .^ 2 .* sr)

    # micro fluxes
    @. fh = dt * u * h - 0.5 * dt^2 * u^2 * sh
    @. fb = dt * u * b - 0.5 * dt^2 * u^2 * sb
    @. fr = dt * u * r - 0.5 * dt^2 * u^2 * sr

    return nothing

end

"""
$(SIGNATURES)

1F3V
"""
function flux_kfvs!(
    fw::AV{<:FN},
    ff::AA{<:FN,3},
    fL::Z,
    fR::Z,
    u::A,
    v::A,
    w::A,
    ω::A,
    dt::Real,
    len = 1.0::Real,
    sfL = zero(fL)::Z,
    sfR = zero(fR)::Z,
) where {Z<:AA{<:FN,3},A<:AA{<:FN,3}}

    # upwind reconstruction
    δ = heaviside.(u)

    f = @. fL * δ + fR * (1.0 - δ)
    sf = @. sfL * δ + sfR * (1.0 - δ)

    # calculate fluxes
    fw[1] = dt * sum(ω .* u .* f) - 0.5 * dt^2 * sum(ω .* u .^ 2 .* sf)
    fw[2] = dt * sum(ω .* u .^ 2 .* f) - 0.5 * dt^2 * sum(ω .* u .^ 3 .* sf)
    fw[3] = dt * sum(ω .* u .* v .* f) - 0.5 * dt^2 * sum(ω .* u .^ 2 .* v .* sf)
    fw[4] = dt * sum(ω .* u .* w .* f) - 0.5 * dt^2 * sum(ω .* u .^ 2 .* w .* sf)
    fw[5] =
        dt * 0.5 * sum(ω .* u .* (u .^ 2 .+ v .^ 2 .+ w .^ 2) .* f) -
        0.5 * dt^2 * 0.5 * sum(ω .* u .^ 2 .* (u .^ 2 .+ v .^ 2 .+ w .^ 2) .* sf)

    @. fw *= len
    @. ff = (dt * u * f - 0.5 * dt^2 * u^2 * sf) * len

    return nothing

end

"""
$(SIGNATURES)

4F1V
"""
function flux_kfvs!(
    fw::AV{<:FN},
    fh0::Y,
    fh1::Y,
    fh2::Y,
    fh3::Y,
    h0L::Z,
    h1L::Z,
    h2L::Z,
    h3L::Z,
    h0R::Z,
    h1R::Z,
    h2R::Z,
    h3R::Z,
    u::A,
    ω::A,
    dt,
    sh0L = zero(h0L)::Z,
    sh1L = zero(h1L)::Z,
    sh2L = zero(h2L)::Z,
    sh3L = zero(h3L)::Z,
    sh0R = zero(h0R)::Z,
    sh1R = zero(h1R)::Z,
    sh2R = zero(h2R)::Z,
    sh3R = zero(h3R)::Z,
) where {Y<:AV{<:FN},Z<:AV{<:FN},A<:AV{<:FN}}

    # upwind reconstruction
    δ = heaviside.(u)

    h0 = @. h0L * δ + h0R * (1.0 - δ)
    h1 = @. h1L * δ + h1R * (1.0 - δ)
    h2 = @. h2L * δ + h2R * (1.0 - δ)
    h3 = @. h3L * δ + h3R * (1.0 - δ)

    sh0 = @. sh0L * δ + sh0R * (1.0 - δ)
    sh1 = @. sh1L * δ + sh1R * (1.0 - δ)
    sh2 = @. sh2L * δ + sh2R * (1.0 - δ)
    sh3 = @. sh3L * δ + sh3R * (1.0 - δ)

    # calculate fluxes
    fw[1] = dt * sum(ω .* u .* h0) - 0.5 * dt^2 * sum(ω .* u .^ 2 .* sh0)
    fw[2] = dt * sum(ω .* u .^ 2 .* h0) - 0.5 * dt^2 * sum(ω .* u .^ 3 .* sh0)
    fw[3] = dt * sum(ω .* u .* h1) - 0.5 * dt^2 * sum(ω .* u .^ 2 .* sh1)
    fw[4] = dt * sum(ω .* u .* h2) - 0.5 * dt^2 * sum(ω .* u .^ 2 .* sh2)
    fw[5] =
        dt * 0.5 * (sum(ω .* u .^ 3 .* h0) + sum(ω .* u .* h3)) -
        0.5 * dt^2 * 0.5 * (sum(ω .* u .^ 4 .* sh0) + sum(ω .* u .^ 2 .* sh3))

    @. fh0 = dt * u * h0 - 0.5 * dt^2 * u^2 * sh0
    @. fh1 = dt * u * h1 - 0.5 * dt^2 * u^2 * sh1
    @. fh2 = dt * u * h2 - 0.5 * dt^2 * u^2 * sh2
    @. fh3 = dt * u * h3 - 0.5 * dt^2 * u^2 * sh3

    return nothing

end

"""
$(SIGNATURES)

Mixture
"""
function flux_kfvs!(
    fw::AM{<:FN},
    fh0::Y,
    fh1::Y,
    fh2::Y,
    fh3::Y,
    h0L::Z,
    h1L::Z,
    h2L::Z,
    h3L::Z,
    h0R::Z,
    h1R::Z,
    h2R::Z,
    h3R::Z,
    u::A,
    ω::A,
    dt,
    sh0L = zero(h0L)::Z,
    sh1L = zero(h1L)::Z,
    sh2L = zero(h2L)::Z,
    sh3L = zero(h3L)::Z,
    sh0R = zero(h0R)::Z,
    sh1R = zero(h1R)::Z,
    sh2R = zero(h2R)::Z,
    sh3R = zero(h3R)::Z,
) where {Y<:AM{<:FN},Z<:AM{<:FN},A<:AM{<:FN}}

    for j in axes(fw, 2)
        _fw = @view fw[:, j]
        _fh0 = @view fh0[:, j]
        _fh1 = @view fh1[:, j]
        _fh2 = @view fh2[:, j]
        _fh3 = @view fh3[:, j]

        flux_kfvs!(
            _fw,
            _fh0,
            _fh1,
            _fh2,
            _fh3,
            h0L[:, j],
            h1L[:, j],
            h2L[:, j],
            h3L[:, j],
            h0R[:, j],
            h1R[:, j],
            h2R[:, j],
            h3R[:, j],
            u[:, j],
            ω[:, j],
            dt,
            sh0L[:, j],
            sh1L[:, j],
            sh2L[:, j],
            sh3L[:, j],
            sh0R[:, j],
            sh1R[:, j],
            sh2R[:, j],
            sh3R[:, j],
        )
    end

    return nothing

end

"""
$(SIGNATURES)

1F2V
"""
function flux_kfvs!(
    fw::AV{<:FN},
    ff::Union{AV{<:FN},AM{<:FN}},
    fL::Z,
    fR::Z,
    u::A,
    v::A,
    ω::A,
    dt,
    len,
    sfL = zero(fL)::Z,
    sfR = zero(fR)::Z,
) where {Z<:Union{AV{<:FN},AM{<:FN}},A<:Union{AV{<:FN},AM{<:FN}}}

    # --- upwind reconstruction ---#
    δ = heaviside.(u)

    f = @. fL * δ + fR * (1.0 - δ)
    sf = @. sfL * δ + sfR * (1.0 - δ)

    # --- calculate fluxes ---#
    fw[1] = dt * sum(ω .* u .* f) - 0.5 * dt^2 * sum(ω .* u .^ 2 .* sf)
    fw[2] = dt * sum(ω .* u .^ 2 .* f) - 0.5 * dt^2 * sum(ω .* u .^ 3 .* sf)
    fw[3] = dt * sum(ω .* v .* u .* f) - 0.5 * dt^2 * sum(ω .* v .* u .^ 2 .* sf)
    fw[4] =
        dt * 0.5 * sum(ω .* u .* (u .^ 2 .+ v .^ 2) .* f) -
        0.5 * dt^2 * 0.5 * sum(ω .* u .^ 2 .* (u .^ 2 .+ v .^ 2) .* sf)
    fw .*= len

    @. ff = (dt * u * f - 0.5 * dt^2 * u^2 * sf) * len

    return nothing

end

"""
$(SIGNATURES)

2F2V
"""
function flux_kfvs!(
    fw::AV{<:FN},
    fh::Y,
    fb::Y,
    hL::Z,
    bL::Z,
    hR::Z,
    bR::Z,
    u::A,
    v::A,
    ω::A,
    dt,
    len,
    shL = zero(hL)::Z,
    sbL = zero(bL)::Z,
    shR = zero(hR)::Z,
    sbR = zero(bR)::Z,
) where {
    Y<:Union{AV{<:FN},AM{<:FN}},
    Z<:Union{AV{<:FN},AM{<:FN}},
    A<:Union{AV{<:FN},AM{<:FN}},
}

    # --- upwind reconstruction ---#
    δ = heaviside.(u)

    h = @. hL * δ + hR * (1.0 - δ)
    b = @. bL * δ + bR * (1.0 - δ)
    sh = @. shL * δ + shR * (1.0 - δ)
    sb = @. sbL * δ + sbR * (1.0 - δ)

    # --- calculate fluxes ---#
    fw[1] = dt * sum(ω .* u .* h) - 0.5 * dt^2 * sum(ω .* u .^ 2 .* sh)
    fw[2] = dt * sum(ω .* u .^ 2 .* h) - 0.5 * dt^2 * sum(ω .* u .^ 3 .* sh)
    fw[3] = dt * sum(ω .* v .* u .* h) - 0.5 * dt^2 * sum(ω .* v .* u .^ 2 .* sh)
    fw[4] =
        dt * 0.5 * (sum(ω .* u .* (u .^ 2 .+ v .^ 2) .* h) + sum(ω .* u .* b)) -
        0.5 *
        dt^2 *
        0.5 *
        (sum(ω .* u .^ 2 .* (u .^ 2 .+ v .^ 2) .* sh) + sum(ω .* u .^ 2 .* sb))
    fw .*= len

    @. fh = (dt * u * h - 0.5 * dt^2 * u^2 * sh) * len
    @. fb = (dt * u * b - 0.5 * dt^2 * u^2 * sb) * len

    return nothing

end

"""
$(SIGNATURES)

3F2V
"""
function flux_kfvs!(
    fw::AV{<:FN},
    fh0::Y,
    fh1::Y,
    fh2::Y,
    h0L::Z,
    h1L::Z,
    h2L::Z,
    h0R::Z,
    h1R::Z,
    h2R::Z,
    u::A,
    v::A,
    ω::A,
    dt,
    len,
    sh0L = zero(h0L)::Z,
    sh1L = zero(h1L)::Z,
    sh2L = zero(h2L)::Z,
    sh0R = zero(h0R)::Z,
    sh1R = zero(h1R)::Z,
    sh2R = zero(h2R)::Z,
) where {Y<:AM{<:FN},Z<:AM{<:FN},A<:AM{<:FN}}

    #--- upwind reconstruction ---#
    δ = heaviside.(u)

    h0 = @. h0L * δ + h0R * (1.0 - δ)
    h1 = @. h1L * δ + h1R * (1.0 - δ)
    h2 = @. h2L * δ + h2R * (1.0 - δ)
    sh0 = @. sh0L * δ + sh0R * (1.0 - δ)
    sh1 = @. sh1L * δ + sh1R * (1.0 - δ)
    sh2 = @. sh2L * δ + sh2R * (1.0 - δ)

    #--- calculate fluxes ---#
    fw[1] = dt * sum(ω .* u .* h0) - 0.5 * dt^2 * sum(ω .* u .^ 2 .* sh0)
    fw[2] = dt * sum(ω .* u .^ 2 .* h0) - 0.5 * dt^2 * sum(ω .* u .^ 3 .* sh0)
    fw[3] = dt * sum(ω .* u .* v .* h0) - 0.5 * dt^2 * sum(ω .* u .^ 2 .* v .* sh0)
    fw[4] = dt * sum(ω .* u .* h1) - 0.5 * dt^2 * sum(ω .* u .^ 2 .* sh1)
    fw[5] =
        dt * 0.5 * (sum(ω .* u .* (u .^ 2 .+ v .^ 2) .* h0) + sum(ω .* u .* h2)) -
        0.5 *
        dt^2 *
        0.5 *
        (sum(ω .* u .^ 2 .* (u .^ 2 .+ v .^ 2) .* sh0) + sum(ω .* u .^ 2 .* sh2))

    fw .*= len
    @. fh0 = (dt * u * h0 - 0.5 * dt^2 * u^2 * sh0) * len
    @. fh1 = (dt * u * h1 - 0.5 * dt^2 * u^2 * sh1) * len
    @. fh2 = (dt * u * h2 - 0.5 * dt^2 * u^2 * sh2) * len

    return nothing

end

"""
$(SIGNATURES)

Mixture
"""
function flux_kfvs!(
    fw::AM{<:FN},
    fh0::Y,
    fh1::Y,
    fh2::Y,
    h0L::Z,
    h1L::Z,
    h2L::Z,
    h0R::Z,
    h1R::Z,
    h2R::Z,
    u::A,
    v::A,
    ω::A,
    dt,
    len,
    sh0L = zero(h0L)::Z,
    sh1L = zero(h1L)::Z,
    sh2L = zero(h2L)::Z,
    sh0R = zero(h0R)::Z,
    sh1R = zero(h1R)::Z,
    sh2R = zero(h2R)::Z,
) where {Y<:AA{<:FN,3},Z<:AA{<:FN,3},A<:AA{<:FN,3}}

    #--- reconstruct initial distribution ---#
    δ = heaviside.(u)

    h0 = @. h0L * δ + h0R * (1.0 - δ)
    h1 = @. h1L * δ + h1R * (1.0 - δ)
    h2 = @. h2L * δ + h2R * (1.0 - δ)

    sh0 = @. sh0L * δ + sh0R * (1.0 - δ)
    sh1 = @. sh1L * δ + sh1R * (1.0 - δ)
    sh2 = @. sh2L * δ + sh2R * (1.0 - δ)

    for j = 1:2
        fw[1, j] =
            dt * sum(ω[:, :, j] .* u[:, :, j] .* h0[:, :, j]) -
            0.5 * dt^2 * sum(ω[:, :, j] .* u[:, :, j] .^ 2 .* sh0[:, :, j])
        fw[2, j] =
            dt * sum(ω[:, :, j] .* u[:, :, j] .^ 2 .* h0[:, :, j]) -
            0.5 * dt^2 * sum(ω[:, :, j] .* u[:, :, j] .^ 3 .* sh0[:, :, j])
        fw[3, j] =
            dt * sum(ω[:, :, j] .* v[:, :, j] .* u[:, :, j] .* h0[:, :, j]) -
            0.5 * dt^2 * sum(ω[:, :, j] .* u[:, :, j] .^ 2 .* v[:, :, j] .* sh0[:, :, j])
        fw[4, j] =
            dt * sum(ω[:, :, j] .* u[:, :, j] .* h1[:, :, j]) -
            0.5 * dt^2 * sum(ω[:, :, j] .* u[:, :, j] .^ 2 .* sh1[:, :, j])
        fw[5, j] =
            dt *
            0.5 *
            (
                sum(
                    ω[:, :, j] .* u[:, :, j] .* (u[:, :, j] .^ 2 .+ v[:, :, j] .^ 2) .*
                    h0[:, :, j],
                ) + sum(ω[:, :, j] .* u[:, :, j] .* h2[:, :, j])
            ) -
            0.5 *
            dt^2 *
            0.5 *
            (
                sum(
                    ω[:, :, j] .* u[:, :, j] .^ 2 .* (u[:, :, j] .^ 2 .+ v[:, :, j] .^ 2) .*
                    sh0[:, :, j],
                ) + sum(ω[:, :, j] .* u[:, :, j] .^ 2 .* sh2[:, :, j])
            )

        @. fh0[:, :, j] =
            dt * u[:, :, j] * h0[:, :, j] - 0.5 * dt^2 * u[:, :, j]^2 * sh0[:, :, j]
        @. fh1[:, :, j] =
            dt * u[:, :, j] * h1[:, :, j] - 0.5 * dt^2 * u[:, :, j]^2 * sh1[:, :, j]
        @. fh2[:, :, j] =
            dt * u[:, :, j] * h2[:, :, j] - 0.5 * dt^2 * u[:, :, j]^2 * sh2[:, :, j]
    end

    @. fw *= len
    @. fh0 *= len
    @. fh1 *= len
    @. fh2 *= len

    return nothing

end
