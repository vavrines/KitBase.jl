# ============================================================
# Structs of Particle
# ============================================================

"""
particle

    Particle(M::AbstractFloat, X::Real, V::AbstractArray, IDX::Integer, T::Real)

- @vars: m, x, v, e, idx, flag, tc

"""
mutable struct Particle{T1,T2,T3} <: AbstractParticle1D

    m::T1
    x::T1
    v::T2
    e::T1
    idx::T3
    flag::T3
    tc::T1

    function Particle(M, X, V, E, IDX, FLAG = zero(IDX), T = zero(M))
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
1D particle

    Particle1D(M::AbstractFloat, X::Real, V::AbstractArray, IDX::Integer, T::Real)

- @vars: m, x, v, e, idx, tc

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
2D particle

    Particle2D(M::AbstractFloat, X::Real, V::AbstractArray, IDX::Integer, T::Real)

- @vars: m, x, v, idx, tb

"""
mutable struct Particle2D{T1,T2,T3} <: AbstractParticle2D

    m::T1
    x::T1
    y::T1
    v::T2
    idx::T3
    idy::T3
    flag::T3
    tc::T1

    function Particle2D(
        M::AbstractFloat,
        X::Real,
        Y::Real,
        V::AbstractArray,
        IDX::Integer,
        IDY::Integer,
        FLAG = zero(IDX),
        T = zero(M)::Real,
    )
        m = deepcopy(M)
        x = deepcopy(X)
        y = deepcopy(Y)
        v = deepcopy(V)
        idx = deepcopy(IDX)
        idy = deepcopy(IDY)
        flag = deepcopy(FLAG)
        tc = deepcopy(T)

        new{typeof(m),typeof(v),typeof(idx)}(m, x, y, v, idx, idy, flag, tc)
    end

end


"""
1D control volume in correspondence with particle simulation

    ControlVolume1D(X::Real, DX::Real, W::AbstractArray, PRIM::AbstractArray)

- @vars: x, dx, w, prim, sw, wg, τ

"""
mutable struct ControlVolumeParticle1D{F,A} <: AbstractControlVolume1D

    x::F
    dx::F

    w::A
    prim::A
    sw::A

    wf::A
    wp::A
    τ::F

    function ControlVolumeParticle1D(
        X::Real,
        DX::Real,
        W::AbstractArray,
        PRIM::AbstractArray,
        TAU = 0.0::Real,
    )
        x = deepcopy(X)
        dx = deepcopy(DX)

        w = deepcopy(W)
        prim = deepcopy(PRIM)
        sw = zeros(eltype(w), axes(w))

        wf = deepcopy(W)
        wp = zero(W)
        τ = deepcopy(TAU)

        new{typeof(x),typeof(w)}(x, dx, w, prim, sw, wf, wp, τ)
    end

end


"""
2D control volume in correspondence with particle simulation

    ControlVolume2D(X::Real, DX::Real, Y::Real, DY::Real, W::AbstractArray, PRIM::AbstractArray)

- @vars: x, y, dx, dy, w, prim, sw, wg, τ

"""
mutable struct ControlVolumeParticle2D{F,A,B} <: AbstractControlVolume2D

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

    function ControlVolumeParticle2D(
        X::Real,
        DX::Real,
        Y::Real,
        DY::Real,
        W::AbstractArray,
        PRIM::AbstractArray,
        TAU = 0.0::Real,
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

        new{typeof(x),typeof(w),typeof(sw)}(x, dx, y, dy, w, prim, sw, wf, wp, τ)

    end

end
