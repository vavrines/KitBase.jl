# ============================================================
# Structs of Particle
# ============================================================

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
    tc::T1

    function Particle1D(M, X, V, E, IDX::Integer, T = zero(M))
        m = deepcopy(M)
        x = deepcopy(X)
        v = deepcopy(V)
        e = deepcopy(E)
        idx = deepcopy(IDX)
        tb = deepcopy(T)

        new{typeof(m),typeof(v),typeof(idx)}(m, x, v, e, idx, tb)
    end

end


"""
2D particle

    Particle2D(M::AbstractFloat, X::Real, V::AbstractArray, IDX::Integer, T::Real)

- @vars: m, x, v, idx, tb

"""
mutable struct Particle2D{T1,T2,T3,T4,T5} <: AbstractParticle2D

    m::T1
    x::T2
    y::T2
    v::T3
    idx::T4
    idy::T4
    tc::T5

    function Particle2D(
        M::AbstractFloat,
        X::Real,
        Y::Real,
        V::AbstractArray,
        IDX::Integer,
        IDY::Integer,
        T = 0.0::Real,
    )
        m = deepcopy(M)
        x = deepcopy(X)
        y = deepcopy(Y)
        v = deepcopy(V)
        idx = deepcopy(IDX)
        idy = deepcopy(IDY)
        tb = deepcopy(T)

        new{typeof(m),typeof(x),typeof(v),typeof(idx),typeof(tb)}(m, x, y, v, idx, idy, tb)
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
