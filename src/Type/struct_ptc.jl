# ============================================================
# Structs of Particle
# ============================================================

"""
1D particle

    Particle1D(M::AbstractFloat, X::Real, V::AbstractArray, IDX::Integer, T::Real)

- @vars: m, x, v, idx, tb

"""
mutable struct Particle1D{T1,T2,T3,T4,T5} <: AbstractParticle1D

    m::T1
    x::T2
    v::T3
    idx::T4
    tb::T5

    function Particle1D(M::AbstractFloat, X::Real, V::AbstractArray, IDX::Integer, T::Real)
        m = deepcopy(M)
        x = deepcopy(X)
        v = deepcopy(V)
        idx = deepcopy(IDX)
        tb = deepcopy(T)

        new{typeof(m),typeof(x),typeof(v),typeof(idx),typeof(tb)}(m, x, v, idx, tb)
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
    tb::T5

    function Particle2D(
        M::AbstractFloat,
        X::Real,
        Y::Real,
        V::AbstractArray,
        IDX::Integer,
        IDY::Integer,
        T::Real,
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

    wg::A
    τ::F

    function ControlVolumeParticle1D(
        X::Real,
        DX::Real,
        W::AbstractArray,
        PRIM::AbstractArray,
    )
        x = deepcopy(X)
        dx = deepcopy(DX)

        w = deepcopy(W)
        prim = deepcopy(PRIM)
        sw = zeros(eltype(w), axes(w))

        wg = deepcopy(W)
        τ = zero(dx)

        new{typeof(x),typeof(w)}(x, dx, w, prim, sw, wg, τ)
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

    wg::A
    τ::F

    function ControlVolumeParticle2D(
        X::Real,
        DX::Real,
        Y::Real,
        DY::Real,
        W::AbstractArray,
        PRIM::AbstractArray,
    )

        x = deepcopy(X)
        dx = deepcopy(DX)
        y = deepcopy(Y)
        dy = deepcopy(DY)

        w = deepcopy(W)
        prim = deepcopy(PRIM)
        sw = zeros(eltype(W), (axes(W)..., Base.OneTo(2)))

        wg = deepcopy(W)
        τ = zero(dx)

        new{typeof(x),typeof(w),typeof(sw)}(x, dx, y, dy, w, prim, sw, wg, τ)

    end

end
