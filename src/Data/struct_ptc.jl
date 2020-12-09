# ============================================================
# Structs of Particle
# ============================================================

"""
1D particle

    Particle1D(M::AbstractFloat, X::Real, V::AbstractArray, IDX::Integer, T::Real)

- @vars: m, x, v, idx, tb

"""
mutable struct Particle1D{T1,T2,T3,T4,T5} <: AbstractControlVolume1D

    m::T1
    x::T2
    v::T3
    idx::T4
    tb::T5

    function Particle1D(
        M::AbstractFloat,
        X::Real,
        V::AbstractArray,
        IDX::Integer,
        T::Real,
    )
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
mutable struct Particle2D{T1,T2,T3,T4,T5} <: AbstractControlVolume1D

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