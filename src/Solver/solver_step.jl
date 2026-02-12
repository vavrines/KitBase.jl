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
function step!(
    KS::AbstractSolverSet,
    cell::AbstractControlVolume,
    faceL::T,
    faceR::T,
    p,
    coll=:bgk;
    ws=nothing,
    kwargs...,
) where {T<:AbstractInterface}
    if isnothing(ws)
        step!(KS, KS.vs, KS.gas, cell, faceL, faceR, p, coll; kwargs...)
    else
        step!(KS, KS.vs, KS.gas, cell, faceL, faceR, p, ws, coll; kwargs...)
    end
end

"""
$(SIGNATURES)

2D
"""
function step!(
    KS::AbstractSolverSet,
    cell::AbstractControlVolume,
    faceL::T,
    faceR::T,
    faceD::T,
    faceU::T,
    p,
    coll=:bgk;
    ws=nothing,
    kwargs...,
) where {T<:AbstractInterface}
    if isnothing(ws)
        step!(KS, KS.vs, KS.gas, cell, faceL, faceR, faceD, faceU, p, coll; kwargs...)
    else
        step!(KS, KS.vs, KS.gas, cell, faceL, faceR, faceD, faceU, p, ws, coll; kwargs...)
    end
end

"""
$(SIGNATURES)

3D
"""
function step!(
    KS::AbstractSolverSet,
    cell::AbstractControlVolume,
    faceXL::T,
    faceXR::T,
    faceYL::T,
    faceYR::T,
    faceZL::T,
    faceZR::T,
    p,
    coll=:bgk;
    ws=nothing,
    kwargs...,
) where {T<:AbstractInterface}
    if isnothing(ws)
        step!(KS, KS.vs, KS.gas, cell, faceXL, faceXR, faceYL, faceYR, faceZL, faceZR, p, coll; kwargs...)
    else
        step!(KS, KS.vs, KS.gas, cell, faceXL, faceXR, faceYL, faceYR, faceZL, faceZR, p, ws, coll; kwargs...)
    end
end

"""
$(SIGNATURES)

1D scalar
"""
function step!(
    KS,
    vs,
    gas::Scalar,
    cell::TC,
    faceL::T,
    faceR::T,
    p,
    coll=:bgk;
    st=step!,
) where {TC<:Union{ControlVolume,ControlVolume1D}, T<:AbstractInterface}
    dt, dx, RES, AVG = p
    return st(cell.w, cell.prim, faceL.fw, faceR.fw, gas.a, dx, RES, AVG)
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
    faceL::T,
    faceR::T,
    p,
    coll=:bgk;
    st=step!,
) where {TC<:Union{ControlVolume,ControlVolume1D}, T<:AbstractInterface}
    dt, dx, RES, AVG = p
    return st(cell.w, cell.prim, faceL.fw, faceR.fw, gas.γ, dx, RES, AVG)
end

"""
$(SIGNATURES)

1D0F with Workspace
"""
function step!(
    KS,
    vs,
    gas::Gas,
    cell::TC,
    faceL::T,
    faceR::T,
    p,
    ws::AbstractWorkspace,
    coll=:bgk;
    st=step!,
) where {TC<:Union{ControlVolume,ControlVolume1D}, T<:AbstractInterface}
    dt, dx, RES, AVG = p
    return st(cell.w, cell.prim, faceL.fw, faceR.fw, gas.γ, dx, RES, AVG, ws)
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
    faceL::T,
    faceR::T,
    p,
    coll=:bgk;
    st=step!,
) where {TC<:Union{ControlVolume,ControlVolume1D}, T<:AbstractInterface}
    dt, dx, RES, AVG = p
    return st(
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
    faceL::T,
    faceR::T,
    p,
    coll=:bgk;
    st=step!,
) where {TC<:Union{ControlVolume1F,ControlVolume1D1F}, T<:AbstractInterface}
    dt, dx, RES, AVG = p
    return st(
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

1D1F1V with Workspace
"""
function step!(
    KS,
    vs::AbstractVelocitySpace1D,
    gas::Gas,
    cell::TC,
    faceL::T,
    faceR::T,
    p,
    ws::AbstractWorkspace,
    coll=:bgk;
    st=step!,
) where {TC<:Union{ControlVolume1F,ControlVolume1D1F}, T<:AbstractInterface}
    dt, dx, RES, AVG = p
    return st(
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
        ws,
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
    faceL::T,
    faceR::T,
    p,
    coll=:bgk;
    st=step!,
) where {TC<:Union{ControlVolume1F,ControlVolume1D1F}, T<:AbstractInterface}
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

1D1F3V
"""
function step!(
    KS,
    vs::AbstractVelocitySpace3D,
    gas::Gas,
    cell::TC,
    faceL::T,
    faceR::T,
    p,
    ws::AbstractWorkspace,
    coll=:bgk;
    st=step!,
) where {TC<:Union{ControlVolume1F,ControlVolume1D1F}, T<:AbstractInterface}
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
            ws,
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
            ws,
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
    faceL::T,
    faceR::T,
    p,
    coll=:bgk;
    st=step!,
) where {TC<:Union{ControlVolume2F,ControlVolume1D2F}, T<:AbstractInterface}
    dt, dx, RES, AVG = p
    return st(
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

1D2F1V
"""
function step!(
    KS,
    vs::AbstractVelocitySpace1D,
    gas::Gas,
    cell::TC,
    faceL::T,
    faceR::T,
    p,
    ws::AbstractWorkspace,
    coll=:bgk;
    st=step!,
) where {TC<:Union{ControlVolume2F,ControlVolume1D2F}, T<:AbstractInterface}
    dt, dx, RES, AVG = p
    return st(
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
        ws,
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
    faceL::T,
    faceR::T,
    p,
    coll=:bgk;
    st=step!,
) where {TC<:Union{ControlVolume2F,ControlVolume1D2F}, T<:AbstractInterface}
    dt, dx, RES, AVG = p
    return st(
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

1D2F1V mixture with workspace
"""
function step!(
    KS,
    vs::AbstractVelocitySpace1D,
    gas::Mixture,
    cell::TC,
    faceL::T,
    faceR::T,
    p,
    ws::AbstractWorkspace,
    coll=:bgk;
    st=step!,
) where {TC<:Union{ControlVolume2F,ControlVolume1D2F}, T<:AbstractInterface}
    dt, dx, RES, AVG = p
    return st(
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
        ws,
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
    faceL::T,
    faceR::T,
    faceD::T,
    faceU::T,
    p,
    coll=:bgk;
    st=step!,
) where {TC<:Union{ControlVolume,ControlVolume2D}, T<:AbstractInterface}
    dt, Δs, RES, AVG = p
    return st(
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

2D0F with Workspace
"""
function step!(
    KS,
    vs,
    gas::Gas,
    cell::TC,
    faceL::T,
    faceR::T,
    faceD::T,
    faceU::T,
    p,
    ws::AbstractWorkspace,
    coll=:bgk;
    st=step!,
) where {TC<:Union{ControlVolume,ControlVolume2D}, T<:AbstractInterface}
    dt, Δs, RES, AVG = p
    return st(
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
        ws,
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
    faceL::T,
    faceR::T,
    faceD::T,
    faceU::T,
    p,
    coll=:bgk;
    st=step!,
) where {TC<:Union{ControlVolume1F,ControlVolume2D1F}, T<:AbstractInterface}
    dt, Δs, RES, AVG = p
    return st(
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

2D1F2V with Workspace
"""
function step!(
    KS,
    vs::AbstractVelocitySpace2D,
    gas::Gas,
    cell::TC,
    faceL::T,
    faceR::T,
    faceD::T,
    faceU::T,
    p,
    ws::AbstractWorkspace,
    coll=:bgk;
    st=step!,
) where {TC<:Union{ControlVolume1F,ControlVolume2D1F}, T<:AbstractInterface}
    dt, Δs, RES, AVG = p
    return st(
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
        ws,
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
    faceL::T,
    faceR::T,
    faceD::T,
    faceU::T,
    p,
    coll=:bgk;
    st=step!,
) where {TC<:Union{ControlVolume2F,ControlVolume2D2F}, T<:AbstractInterface}
    dt, Δs, RES, AVG = p
    return st(
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

"""
$(SIGNATURES)

2D2F2V with Workspace
"""
function step!(
    KS,
    vs::AbstractVelocitySpace2D,
    gas::Gas,
    cell::TC,
    faceL::T,
    faceR::T,
    faceD::T,
    faceU::T,
    p,
    ws::AbstractWorkspace,
    coll=:bgk;
    st=step!,
) where {TC<:Union{ControlVolume2F,ControlVolume2D2F}, T<:AbstractInterface}
    dt, Δs, RES, AVG = p
    return st(
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
        ws,
        coll,
    )
end

"""
$(SIGNATURES)

3D0F
"""
function step!(
    KS,
    vs,
    gas::Gas,
    cell::TC,
    faceXL::T,
    faceXR::T,
    faceYL::T,
    faceYR::T,
    faceZL::T,
    faceZR::T,
    p,
    coll=:bgk;
    st=step!,
) where {TC<:Union{ControlVolume,ControlVolume3D}, T<:AbstractInterface}
    dt, Δv, RES, AVG = p
    return st(
        cell.w,
        cell.prim,
        faceXL.fw,
        faceXR.fw,
        faceYL.fw,
        faceYR.fw,
        faceZL.fw,
        faceZR.fw,
        gas.γ,
        Δv,
        RES,
        AVG,
        coll,
    )
end

"""
$(SIGNATURES)

3D0F with Workspace
"""
function step!(
    KS,
    vs,
    gas::Gas,
    cell::TC,
    faceXL::T,
    faceXR::T,
    faceYL::T,
    faceYR::T,
    faceZL::T,
    faceZR::T,
    p,
    ws::AbstractWorkspace,
    coll=:bgk;
    st=step!,
) where {TC<:Union{ControlVolume,ControlVolume3D}, T<:AbstractInterface}
    dt, Δv, RES, AVG = p
    return st(
        cell.w,
        cell.prim,
        faceXL.fw,
        faceXR.fw,
        faceYL.fw,
        faceYR.fw,
        faceZL.fw,
        faceZR.fw,
        gas.γ,
        Δv,
        RES,
        AVG,
        ws,
        coll,
    )
end