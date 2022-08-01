# ============================================================
# Time Steppers
# Low-level methods can be found in the following files:
# step_0f.jl, step_1f.jl, step_2f.jl, step_3f.jl, step_4f.jl
# ============================================================

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
    coll = :bgk;
    st = step!,
) where {T<:AbstractInterface} = step!(KS, KS.vs, KS.gas, cell, faceL, faceR, p, coll; st = st)

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
    coll = :bgk;
    st = step!,
) where {T<:AbstractInterface} =
    step!(KS, KS.vs, KS.gas, cell, faceL, faceR, faceD, faceU, p, coll; st = st)

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
    coll = :bgk;
    st = step!,
) where {TC<:Union{ControlVolume,ControlVolume1D}}
    dt, dx, RES, AVG = p
    st(cell.w, cell.prim, faceL.fw, faceR.fw, gas.a, dx, RES, AVG)
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
    coll = :bgk;
    st = step!,
) where {TC<:Union{ControlVolume,ControlVolume1D}}
    dt, dx, RES, AVG = p
    st(cell.w, cell.prim, faceL.fw, faceR.fw, gas.γ, dx, RES, AVG)
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
    coll = :bgk;
    st = step!,
) where {TC<:Union{ControlVolume,ControlVolume1D}}
    dt, dx, RES, AVG = p
    st(
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
    coll = :bgk;
    st = step!,
) where {TC<:Union{ControlVolume1F,ControlVolume1D1F}}
    dt, dx, RES, AVG = p
    st(
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
    coll = :bgk;
    st = step!,
) where {TC<:Union{ControlVolume1F,ControlVolume1D1F}}
    dt, dx, RES, AVG = p
    if coll == :fsm
        st(
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
        st(
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
    coll = :bgk;
    st = step!,
) where {TC<:Union{ControlVolume2F,ControlVolume1D2F}}
    dt, dx, RES, AVG = p
    st(
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
    coll = :bgk;
    st = step!,
) where {TC<:Union{ControlVolume2F,ControlVolume1D2F}}
    dt, dx, RES, AVG = p
    st(
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
    coll = :bgk;
    st = step!,
) where {TC<:Union{ControlVolume,ControlVolume2D}}
    dt, Δs, RES, AVG = p
    st(
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
    coll = :bgk;
    st = step!,
) where {TC<:Union{ControlVolume1F,ControlVolume2D1F}}
    dt, Δs, RES, AVG = p
    st(
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
    coll = :bgk;
    st = step!,
) where {TC<:Union{ControlVolume2F,ControlVolume2D2F}}
    dt, Δs, RES, AVG = p
    st(
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
