# ============================================================
# Structs of Particle
# ============================================================

"""
$(TYPEDEF)

Struct of arrays for particle simulation

## Fields

$(FIELDS)
"""
mutable struct Particle{T1,T2,T3,T4,T5} <: AbstractParticle

    m::T1
    x::T2
    v::T3
    e::T1
    idx::T4
    ref::T4
    flag::T5
    tc::T1

    function Particle(
        M,
        X,
        V,
        E,
        IDX,
        REF = zeros(eltype(IDX), axes(IDX, 1)),
        FLAG = zeros(eltype(IDX), axes(IDX, 1)),
        T = zero(M),
    )
        m = deepcopy(M)
        x = deepcopy(X)
        v = deepcopy(V)
        e = deepcopy(E)
        idx = deepcopy(IDX)
        ref = deepcopy(REF)
        flag = deepcopy(FLAG)
        tc = deepcopy(T)

        new{typeof(m),typeof(x),typeof(v),typeof(idx),typeof(flag)}(
            m,
            x,
            v,
            e,
            idx,
            ref,
            flag,
            tc,
        )
    end

end


"""
$(TYPEDEF)

1D particle

## Fields

$(FIELDS)
"""
mutable struct Particle1D{T1,T2,T3} <: AbstractParticle1D

    m::T1
    x::T1
    v::T2
    e::T1
    idx::T3
    flag::T3
    tc::T1

    function Particle1D(M, X, V, E, IDX::Integer, FLAG = zero(IDX), T = zero(M))
        m = deepcopy(M)
        x = deepcopy(X)
        v = deepcopy(V)
        e = deepcopy(E)
        idx = deepcopy(IDX)
        flag = deepcopy(FLAG)
        tc = deepcopy(T)

        new{typeof(m),typeof(v),typeof(idx)}(m, x, v, e, idx, flag, tc)
    end

end


"""
$(TYPEDEF)

2D particle

## Fields

$(FIELDS)
"""
mutable struct Particle2D{T1,T2,T3} <: AbstractParticle2D

    m::T1
    x::T1
    y::T1
    v::T2
    e::T1
    idx::T3
    idy::T3
    flag::T3
    tc::T1

    function Particle2D(
        M::FN,
        X::Real,
        Y::Real,
        V::AA,
        E,
        IDX::Integer,
        IDY::Integer,
        FLAG = zero(IDX),
        T = zero(M)::Real,
    )
        m = deepcopy(M)
        x = deepcopy(X)
        y = deepcopy(Y)
        v = deepcopy(V)
        e = deepcopy(E)
        idx = deepcopy(IDX)
        idy = deepcopy(IDY)
        flag = deepcopy(FLAG)
        tc = deepcopy(T)

        new{typeof(m),typeof(v),typeof(idx)}(m, x, y, v, e, idx, idy, flag, tc)
    end

end


"""
$(TYPEDEF)

1D control volume in correspondence with particle simulation

## Fields

$(FIELDS)
"""
mutable struct ControlVolumeParticle1D{F,A,I<:Integer} <: AbstractControlVolume1D

    x::F
    dx::F

    w::A
    prim::A
    sw::A

    wf::A
    wp::A
    τ::F

    np::I
    ip::I
    vrmax::F
    remainder::F

    function ControlVolumeParticle1D(
        X::Real,
        DX::Real,
        W::AA,
        PRIM::AA,
        TAU = 0.0::Real,
        NP = 0::Integer,
        IP = 0::Integer,
        VR = 3.0 * sqrt(PRIM[end])::Real,
        RE = 0::Integer,
    )
        x = deepcopy(X)
        dx = deepcopy(DX)

        w = deepcopy(W)
        prim = deepcopy(PRIM)
        sw = zeros(eltype(w), axes(w))

        wf = deepcopy(W)
        wp = zero(W)
        τ = deepcopy(TAU)

        np = deepcopy(NP)
        ip = deepcopy(IP)
        vrmax = deepcopy(VR)
        remainder = deepcopy(RE)

        new{typeof(x),typeof(w),typeof(np)}(
            x,
            dx,
            w,
            prim,
            sw,
            wf,
            wp,
            τ,
            np,
            ip,
            vrmax,
            remainder,
        )
    end

end


"""
$(TYPEDEF)

2D control volume in correspondence with particle simulation

## Fields

$(FIELDS)
"""
mutable struct ControlVolumeParticle2D{F,A,B,I} <: AbstractControlVolume2D

    x::F
    y::F
    dx::F
    dy::F

    w::A
    prim::A
    sw::B

    wf::A
    wp::A
    τ::F

    np::I
    ip::I
    vrmax::F
    remainder::F

    function ControlVolumeParticle2D(
        X::Real,
        DX::Real,
        Y::Real,
        DY::Real,
        W::AA,
        PRIM::AA,
        TAU = 0.0::Real,
        NP = 0::Integer,
        IP = 0::Integer,
        VR = 3.0 * sqrt(PRIM[end])::Real,
        RE = 0::Integer,
    )

        x = deepcopy(X)
        dx = deepcopy(DX)
        y = deepcopy(Y)
        dy = deepcopy(DY)

        w = deepcopy(W)
        prim = deepcopy(PRIM)
        sw = zeros(eltype(W), (axes(W)..., Base.OneTo(2)))

        wf = deepcopy(W)
        wp = zero(W)
        τ = deepcopy(TAU)

        np = deepcopy(NP)
        ip = deepcopy(IP)
        vrmax = deepcopy(VR)
        remainder = deepcopy(RE)

        new{typeof(x),typeof(w),typeof(sw),typeof(np)}(
            x,
            dx,
            y,
            dy,
            w,
            prim,
            sw,
            wf,
            wp,
            τ,
            np,
            ip,
            vrmax,
            remainder,
        )

    end

end
