# ============================================================
# Structs of Control Volume
# Solver stores data in arrays of struct (AoS)
# ============================================================

# ------------------------------------------------------------
# 1D
# ------------------------------------------------------------

"""
    mutable struct ControlVolume1D{F,A,B} <: AbstractControlVolume1D
        x::F
        dx::F

        w::A
        prim::B
        sw::A
    end

1D control volume with no distribution function

"""
mutable struct ControlVolume1D{F,A,B} <: AbstractControlVolume1D
    x::F
    dx::F

    w::A
    prim::B
    sw::A
end

function ControlVolume1D(X::T1, DX::T1, W::T2, PRIM::T3) where {T1<:Real,T2,T3<:AbstractArray}
    x = deepcopy(X)
    dx = deepcopy(DX)

    w = deepcopy(W)
    prim = deepcopy(PRIM)
    sw = zero(W)

    return ControlVolume1D{typeof(x),typeof(w),typeof(prim)}(x, dx, w, prim, sw)
end

function Base.show(io::IO, ctr::ControlVolume1D{F,A,B}) where {F,A,B}
    print(
        io,
        "ControlVolume1D{$F,$A,$B}\n",
        "center: $(ctr.x)\n",
        "interval: $(ctr.dx)\n",
        "conservative vars: $(ctr.w)\n",
        "primitive vars: $(ctr.prim)\n",
        "conservative slopes: $(ctr.sw)\n",
    )
end


"""
    mutable struct ControlVolume1D1F{F,A,B} <: AbstractControlVolume1D
        x::F
        dx::F

        w::A
        prim::A
        sw::A

        f::B
        sf::B
    end

1D control volume with 1 distribution function

"""
mutable struct ControlVolume1D1F{F,A,B} <: AbstractControlVolume1D
    x::F
    dx::F

    w::A
    prim::A
    sw::A

    f::B
    sf::B
end

function ControlVolume1D1F(
    X,
    DX,
    W::T1,
    PRIM::T1,
    F::T2,
) where {T1<:AbstractArray,T2<:AbstractArray}
    x = deepcopy(X)
    dx = deepcopy(DX)

    w = deepcopy(W)
    prim = deepcopy(PRIM)
    sw = zero(w)

    f = deepcopy(F)
    sf = zero(f)

    return ControlVolume1D1F{typeof(x),typeof(w),typeof(f)}(x, dx, w, prim, sw, f, sf)
end

function Base.show(io::IO, ctr::ControlVolume1D1F{F,A,B}) where {F,A,B}
    print(
        io,
        "ControlVolume1D1F{$F,$A,$B}\n",
        "center: $(ctr.x)\n",
        "interval: $(ctr.dx)\n",
        "conservative vars: $(ctr.w)\n",
        "primitive vars: $(ctr.prim)\n",
        "conservative slopes: $(ctr.sw)\n",
        "pdf vars: f\n",
        "pdf slopes: sf\n",
    )
end


"""
    mutable struct ControlVolume1D2F{F,A,B} <: AbstractControlVolume1D
        x::F
        dx::F

        w::A
        prim::A
        sw::A

        h::B
        b::B
        sh::B
        sb::B
    end

1D control volume with 2 distribution functions

"""
mutable struct ControlVolume1D2F{F,A,B} <: AbstractControlVolume1D
    x::F
    dx::F

    w::A
    prim::A
    sw::A

    h::B
    b::B
    sh::B
    sb::B
end

function ControlVolume1D2F(
    X,
    DX,
    W::T1,
    PRIM::T1,
    H::T2,
    B::T2,
) where {T1<:AbstractArray,T2<:AbstractArray}

    x = deepcopy(X)
    dx = deepcopy(DX)

    w = deepcopy(W)
    prim = deepcopy(PRIM)
    sw = zero(w)

    h = deepcopy(H)
    b = deepcopy(B)
    sh = zero(h)
    sb = zero(b)

    return ControlVolume1D2F{typeof(x),typeof(w),typeof(h)}(
        x,
        dx,
        w,
        prim,
        sw,
        h,
        b,
        sh,
        sb,
    )

end

function Base.show(io::IO, ctr::ControlVolume1D2F{F,A,B}) where {F,A,B}
    print(
        io,
        "ControlVolume1D2F{$F,$A,$B}\n",
        "center: $(ctr.x)\n",
        "interval: $(ctr.dx)\n",
        "conservative vars: $(ctr.w)\n",
        "primitive vars: $(ctr.prim)\n",
        "conservative slopes: $(ctr.sw)\n",
        "pdf vars: h, b\n",
        "pdf slopes: sh, sb\n",
    )
end


"""
    mutable struct ControlVolume1D3F{F,A,B,C,D,E} <: AbstractControlVolume1D
        x::F
        dx::F

        w::A
        prim::A
        sw::A

        h0::B
        h1::B
        h2::B
        sh0::B
        sh1::B
        sh2::B

        E::C
        B::C
        ϕ::D
        ψ::D
        lorenz::E
    end

1D control volume with 3 distribution functions

"""
mutable struct ControlVolume1D3F{F,A,B,C,D,E} <: AbstractControlVolume1D
    x::F
    dx::F

    w::A
    prim::A
    sw::A

    h0::B
    h1::B
    h2::B
    sh0::B
    sh1::B
    sh2::B

    E::C
    B::C
    ϕ::D
    ψ::D
    lorenz::E
end

# deterministic
function ControlVolume1D3F(
    X::Real,
    DX::Real,
    W::AbstractArray{<:Real,2},
    PRIM::AbstractArray{<:Real,2},
    H0::AbstractArray{<:AbstractFloat,3},
    H1::AbstractArray{<:AbstractFloat,3},
    H2::AbstractArray{<:AbstractFloat,3},
    E0::AbstractArray{<:AbstractFloat,1},
    B0::AbstractArray{<:AbstractFloat,1},
    L::AbstractArray{<:AbstractFloat,2},
)
    x = deepcopy(X)
    dx = deepcopy(DX)

    w = deepcopy(W)
    prim = deepcopy(PRIM)
    sw = zeros(eltype(H0), axes(W))

    h0 = deepcopy(H0)
    h1 = deepcopy(H1)
    h2 = deepcopy(H2)
    sh0 = zeros(eltype(H0), axes(H0))
    sh1 = zeros(eltype(H1), axes(H1))
    sh2 = zeros(eltype(H2), axes(H2))

    E = deepcopy(E0)
    B = deepcopy(B0)
    ϕ = 0.0
    ψ = 0.0
    lorenz = deepcopy(L)

    return ControlVolume1D3F{
        typeof(x),
        typeof(w),
        typeof(h0),
        typeof(E),
        typeof(ϕ),
        typeof(lorenz),
    }(
        x,
        dx,
        w,
        prim,
        sw,
        h0,
        h1,
        h2,
        sh0,
        sh1,
        sh2,
        E,
        B,
        ϕ,
        ψ,
        lorenz,
    )
end

# Rykov
function ControlVolume1D3F(
    X::Real,
    DX::Real,
    W::AbstractArray{<:Real,1},
    PRIM::AbstractArray{<:Real,1},
    H0::AbstractArray{<:AbstractFloat,1},
    H1::AbstractArray{<:AbstractFloat,1},
    H2::AbstractArray{<:AbstractFloat,1},
)
    x = deepcopy(X)
    dx = deepcopy(DX)

    w = deepcopy(W)
    prim = deepcopy(PRIM)
    sw = zeros(eltype(H0), axes(W))

    h0 = deepcopy(H0)
    h1 = deepcopy(H1)
    h2 = deepcopy(H2)
    sh0 = zeros(eltype(H0), axes(H0))
    sh1 = zeros(eltype(H1), axes(H1))
    sh2 = zeros(eltype(H2), axes(H2))

    E = nothing
    B = nothing
    ϕ = nothing
    ψ = nothing
    lorenz = nothing

    return ControlVolume1D3F{
        typeof(x),
        typeof(w),
        typeof(h0),
        typeof(E),
        typeof(ϕ),
        typeof(lorenz),
    }(
        x,
        dx,
        w,
        prim,
        sw,
        h0,
        h1,
        h2,
        sh0,
        sh1,
        sh2,
        E,
        B,
        ϕ,
        ψ,
        lorenz,
    )
end

# stochastic
function ControlVolume1D3F(
    X::Real,
    DX::Real,
    W::AbstractArray{<:Real,3},
    PRIM::AbstractArray{<:Real,3},
    H0::AbstractArray{<:AbstractFloat,4},
    H1::AbstractArray{<:AbstractFloat,4},
    H2::AbstractArray{<:AbstractFloat,4},
    E0::AbstractArray{<:AbstractFloat,2},
    B0::AbstractArray{<:AbstractFloat,2},
    L::AbstractArray{<:AbstractFloat,3},
)
    x = deepcopy(X)
    dx = deepcopy(DX)

    w = deepcopy(W)
    prim = deepcopy(PRIM)
    sw = zeros(eltype(H0), axes(W))

    h0 = deepcopy(H0)
    h1 = deepcopy(H1)
    h2 = deepcopy(H2)
    sh0 = zeros(eltype(H0), axes(H0))
    sh1 = zeros(eltype(H1), axes(H1))
    sh2 = zeros(eltype(H2), axes(H2))

    E = deepcopy(E0)
    B = deepcopy(B0)
    ϕ = zeros(axes(E, 2)) # here is difference
    ψ = zeros(axes(B, 2))
    lorenz = deepcopy(L)

    return ControlVolume1D3F{
        typeof(x),
        typeof(w),
        typeof(h0),
        typeof(E),
        typeof(ϕ),
        typeof(lorenz),
    }(
        x,
        dx,
        w,
        prim,
        sw,
        h0,
        h1,
        h2,
        sh0,
        sh1,
        sh2,
        E,
        B,
        ϕ,
        ψ,
        lorenz,
    )
end

function Base.show(io::IO, ctr::ControlVolume1D3F{F,A,B,C,D,E}) where {F,A,B,C,D,E}
    print(
        io,
        "ControlVolume1D3F{$F,$A,$B,$C,$D,$E}\n",
        "center: $(ctr.x)\n",
        "interval: $(ctr.dx)\n",
        "conservative vars: $(ctr.w)\n",
        "primitive vars: $(ctr.prim)\n",
        "conservative slopes: $(ctr.sw)\n",
        "pdf vars: h0, h1, h2\n",
        "pdf slopes: sh0, sh1, sh2\n",
        "electric field: $(ctr.E)\n",
        "magnetic field: $(ctr.B)\n",
        "Lorenz force: $(ctr.lorenz)\n",
    )
end


"""
    mutable struct ControlVolume1D4F{F,A,B,C,D,E} <: AbstractControlVolume1D
        x::F
        dx::F

        w::A
        prim::A
        sw::A

        h0::B
        h1::B
        h2::B
        h3::B
        sh0::B
        sh1::B
        sh2::B
        sh3::B

        E::C
        B::C
        ϕ::D
        ψ::D
        lorenz::E
    end

1D control volume with 4 distribution functions

"""
mutable struct ControlVolume1D4F{F,A,B,C,D,E} <: AbstractControlVolume1D
    x::F
    dx::F

    w::A
    prim::A
    sw::A

    h0::B
    h1::B
    h2::B
    h3::B
    sh0::B
    sh1::B
    sh2::B
    sh3::B

    E::C
    B::C
    ϕ::D
    ψ::D
    lorenz::E
end

#--- deterministic ---#
function ControlVolume1D4F(
    X::Real,
    DX::Real,
    W::AbstractArray{<:Real,2},
    PRIM::AbstractArray{<:Real,2},
    H0::AbstractArray{<:AbstractFloat,2},
    H1::AbstractArray{Float64,2},
    H2::AbstractArray{Float64,2},
    H3::AbstractArray{Float64,2},
    E0::AbstractArray{Float64,1},
    B0::AbstractArray{Float64,1},
    L::AbstractArray{Float64,2},
)

    x = deepcopy(X)
    dx = deepcopy(DX)

    w = deepcopy(W)
    prim = deepcopy(PRIM)
    sw = zeros(typeof(W[1]), axes(W))

    h0 = deepcopy(H0)
    h1 = deepcopy(H1)
    h2 = deepcopy(H2)
    h3 = deepcopy(H3)
    sh0 = zeros(typeof(H0[1]), axes(H0))
    sh1 = zeros(typeof(H1[1]), axes(H1))
    sh2 = zeros(typeof(H2[1]), axes(H2))
    sh3 = zeros(typeof(H3[1]), axes(H3))

    E = deepcopy(E0)
    B = deepcopy(B0)
    ϕ = 0.0
    ψ = 0.0
    lorenz = deepcopy(L)

    return ControlVolume1D4F{
        typeof(x),
        typeof(w),
        typeof(h0),
        typeof(E),
        typeof(ϕ),
        typeof(lorenz),
    }(
        x,
        dx,
        w,
        prim,
        sw,
        h0,
        h1,
        h2,
        h3,
        sh0,
        sh1,
        sh2,
        sh3,
        E,
        B,
        ϕ,
        ψ,
        lorenz,
    )

end

#--- uncertainty quantification ---#
function ControlVolume1D4F(
    X::Real,
    DX::Real,
    W::AbstractArray{<:Real,3},
    PRIM::AbstractArray{<:Real,3},
    H0::AbstractArray{<:AbstractFloat,3},
    H1::AbstractArray{Float64,3},
    H2::AbstractArray{Float64,3},
    H3::AbstractArray{Float64,3},
    E0::AbstractArray{Float64,2},
    B0::AbstractArray{Float64,2},
    L::AbstractArray{Float64,3},
)

    x = deepcopy(X)
    dx = deepcopy(DX)

    w = deepcopy(W)
    prim = deepcopy(PRIM)
    sw = zeros(typeof(W[1]), axes(W))

    h0 = deepcopy(H0)
    h1 = deepcopy(H1)
    h2 = deepcopy(H2)
    h3 = deepcopy(H3)
    sh0 = zeros(typeof(H0[1]), axes(H0))
    sh1 = zeros(typeof(H1[1]), axes(H1))
    sh2 = zeros(typeof(H2[1]), axes(H2))
    sh3 = zeros(typeof(H3[1]), axes(H3))

    E = deepcopy(E0)
    B = deepcopy(B0)
    ϕ = zeros(axes(E, 2))
    ψ = zeros(axes(B, 2))
    lorenz = deepcopy(L)

    return ControlVolume1D4F{
        typeof(x),
        typeof(w),
        typeof(h0),
        typeof(E),
        typeof(ϕ),
        typeof(lorenz),
    }(
        x,
        dx,
        w,
        prim,
        sw,
        h0,
        h1,
        h2,
        h3,
        sh0,
        sh1,
        sh2,
        sh3,
        E,
        B,
        ϕ,
        ψ,
        lorenz,
    )

end

function Base.show(io::IO, ctr::ControlVolume1D4F{F,A,B,C,D,E}) where {F,A,B,C,D,E}
    print(
        io,
        "ControlVolume1D4F{$F,$A,$B,$C,$D,$E}\n",
        "center: $(ctr.x)\n",
        "interval: $(ctr.dx)\n",
        "conservative vars: $(ctr.w)\n",
        "primitive vars: $(ctr.prim)\n",
        "conservative slopes: $(ctr.sw)\n",
        "pdf vars: h0, h1, h2, h3\n",
        "pdf slopes: sh0, sh1, sh2, sh3\n",
        "electric field: $(ctr.E)\n",
        "magnetic field: $(ctr.B)\n",
        "Lorenz force: $(ctr.lorenz)\n",
    )
end

# ------------------------------------------------------------
# 2D
# ------------------------------------------------------------

"""
    mutable struct ControlVolume2D{F,A,B} <: AbstractControlVolume2D
        x::F
        y::F
        dx::F
        dy::F

        w::A
        prim::A
        sw::B
    end

2D control volume with no distribution function

"""
mutable struct ControlVolume2D{F,A,B} <: AbstractControlVolume2D
    x::F
    y::F
    dx::F
    dy::F

    w::A
    prim::A
    sw::B
end

function ControlVolume2D(
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

    return ControlVolume2D{typeof(x),typeof(w),typeof(sw)}(x, dx, y, dy, w, prim, sw)
end

function Base.show(io::IO, ctr::ControlVolume2D{F,A,B}) where {F,A,B}
    print(
        io,
        "ControlVolume2D{$F,$A,$B}\n",
        "center: ($(ctr.x),$(ctr.y))\n",
        "interval: ($(ctr.dx),$(ctr.dy))\n",
        "conservative vars: $(ctr.w)\n",
        "primitive vars: $(ctr.prim)\n",
        "conservative slopes: $(ctr.sw)\n",
    )
end


"""
    mutable struct ControlVolume2D1F{F,A,B,C,D} <: AbstractControlVolume2D
        x::F
        y::F
        dx::F
        dy::F

        w::A
        prim::A
        sw::B

        f::C
        sf::D
    end

2D control volume with 1 distribution function

"""
mutable struct ControlVolume2D1F{F,A,B,C,D} <: AbstractControlVolume2D
    x::F
    y::F
    dx::F
    dy::F

    w::A
    prim::A
    sw::B

    f::C
    sf::D
end

function ControlVolume2D1F(
    X::Real,
    DX::Real,
    Y::Real,
    DY::Real,
    W::AbstractArray,
    PRIM::AbstractArray,
    F::AbstractArray,
)

    x = deepcopy(X)
    dx = deepcopy(DX)
    y = deepcopy(Y)
    dy = deepcopy(DY)

    w = deepcopy(W)
    prim = deepcopy(PRIM)
    sw = zeros(eltype(W), (axes(W)..., Base.OneTo(2)))

    f = deepcopy(F)
    sf = zeros(eltype(F), (axes(F)..., Base.OneTo(2)))

    return ControlVolume2D1F{typeof(x),typeof(w),typeof(sw),typeof(f),typeof(sf)}(
        x,
        dx,
        y,
        dy,
        w,
        prim,
        sw,
        f,
        sf,
    )

end

function Base.show(io::IO, ctr::ControlVolume2D1F{F,A,B,C,D}) where {F,A,B,C,D}
    print(
        io,
        "ControlVolume2D1F{$F,$A,$B,$C,$D}\n",
        "center: ($(ctr.x),$(ctr.y))\n",
        "interval: ($(ctr.dx),$(ctr.dy))\n",
        "conservative vars: $(ctr.w)\n",
        "primitive vars: $(ctr.prim)\n",
        "conservative slopes: $(ctr.sw)\n",
        "pdf vars: f\n",
        "pdf slopes: sf\n",
    )
end


"""
    mutable struct ControlVolume2D2F{F,A,B,C,D} <: AbstractControlVolume2D
        x::F
        y::F
        dx::F
        dy::F

        w::A
        prim::A
        sw::B

        h::C
        b::C
        sh::D
        sb::D
    end

2D control volume with 2 distribution functions

"""
mutable struct ControlVolume2D2F{F,A,B,C,D} <: AbstractControlVolume2D
    x::F
    y::F
    dx::F
    dy::F

    w::A
    prim::A
    sw::B

    h::C
    b::C
    sh::D
    sb::D
end

function ControlVolume2D2F(
    X::Real,
    DX::Real,
    Y::Real,
    DY::Real,
    W::AbstractArray,
    PRIM::AbstractArray,
    H::AbstractArray,
    B::AbstractArray,
)

    x = deepcopy(X)
    dx = deepcopy(DX)
    y = deepcopy(Y)
    dy = deepcopy(DY)

    w = deepcopy(W)
    prim = deepcopy(PRIM)
    sw = zeros(eltype(W), (axes(W)..., Base.OneTo(2)))

    h = deepcopy(H)
    b = deepcopy(B)
    sh = zeros(eltype(H), (axes(H)..., Base.OneTo(2)))
    sb = zeros(eltype(B), (axes(B)..., Base.OneTo(2)))

    return ControlVolume2D2F{typeof(x),typeof(w),typeof(sw),typeof(h),typeof(sh)}(
        x,
        dx,
        y,
        dy,
        w,
        prim,
        sw,
        h,
        b,
        sh,
        sb,
    )

end

function Base.show(io::IO, ctr::ControlVolume2D2F{F,A,B,C,D}) where {F,A,B,C,D}
    print(
        io,
        "ControlVolume2D2F{$F,$A,$B,$C,$D}\n",
        "center: ($(ctr.x),$(ctr.y))\n",
        "interval: ($(ctr.dx),$(ctr.dy))\n",
        "conservative vars: $(ctr.w)\n",
        "primitive vars: $(ctr.prim)\n",
        "conservative slopes: $(ctr.sw)\n",
        "pdf vars: h, b\n",
        "pdf slopes: sh, sb\n",
    )
end


"""
    mutable struct ControlVolume2D3F{T1,T2,T3,T4,T5,T6,T7,T8} <: AbstractControlVolume2D
        x::T1
        y::T1
        dx::T1
        dy::T1

        w::T2
        prim::T2
        sw::T3

        h0::T4
        h1::T4
        h2::T4
        sh0::T5
        sh1::T5
        sh2::T5

        E::T6
        B::T6
        ϕ::T7
        ψ::T7
        lorenz::T8
    end

2D control volume with 3 distribution functions

"""
mutable struct ControlVolume2D3F{T1,T2,T3,T4,T5,T6,T7,T8} <: AbstractControlVolume2D
    x::T1
    y::T1
    dx::T1
    dy::T1

    w::T2
    prim::T2
    sw::T3

    h0::T4
    h1::T4
    h2::T4
    sh0::T5
    sh1::T5
    sh2::T5

    E::T6
    B::T6
    ϕ::T7
    ψ::T7
    lorenz::T8
end

#--- deterministic & stochastic ---#
function ControlVolume2D3F(
    X::Real,
    DX::Real,
    Y::Real,
    DY::Real,
    W::AbstractArray,
    PRIM::AbstractArray,
    H0::AbstractArray,
    H1::AbstractArray,
    H2::AbstractArray,
    E0::AbstractArray,
    B0::AbstractArray,
    L::AbstractArray,
)
    x = deepcopy(X)
    dx = deepcopy(DX)
    y = deepcopy(Y)
    dy = deepcopy(DY)

    w = deepcopy(W)
    prim = deepcopy(PRIM)
    sw = zeros(eltype(W), (axes(W)..., Base.OneTo(2))) # 2D

    h0 = deepcopy(H0)
    h1 = deepcopy(H1)
    h2 = deepcopy(H2)
    sh0 = zeros(eltype(H0), (axes(H0)..., Base.OneTo(2)))
    sh1 = zeros(eltype(H1), (axes(H1)..., Base.OneTo(2)))
    sh2 = zeros(eltype(H2), (axes(H2)..., Base.OneTo(2)))

    E = deepcopy(E0)
    B = deepcopy(B0)
    ϕ = 0.0
    ψ = 0.0
    lorenz = deepcopy(L)

    return ControlVolume2D3F{
        typeof(x),
        typeof(w),
        typeof(sw),
        typeof(h0),
        typeof(sh2),
        typeof(E),
        typeof(ϕ),
        typeof(lorenz),
    }(
        x,
        dx,
        y,
        dy,
        w,
        prim,
        sw,
        h0,
        h1,
        h2,
        sh0,
        sh1,
        sh2,
        E,
        B,
        ϕ,
        ψ,
        lorenz,
    )
end

#--- Rykov ---#
function ControlVolume2D3F(
    X::Real,
    DX::Real,
    Y::Real,
    DY::Real,
    W::AbstractArray{<:Real,1},
    PRIM::AbstractArray{<:Real,1},
    H0::AbstractArray{<:AbstractFloat,2},
    H1::AbstractArray{<:AbstractFloat,2},
    H2::AbstractArray{<:AbstractFloat,2},
)
    x = deepcopy(X)
    dx = deepcopy(DX)
    y = deepcopy(Y)
    dy = deepcopy(DY)

    w = deepcopy(W)
    prim = deepcopy(PRIM)
    sw = zeros(eltype(W), (axes(W)..., Base.OneTo(2))) # 2D

    h0 = deepcopy(H0)
    h1 = deepcopy(H1)
    h2 = deepcopy(H2)
    sh0 = zeros(eltype(H0), (axes(H0)..., Base.OneTo(2)))
    sh1 = zeros(eltype(H1), (axes(H1)..., Base.OneTo(2)))
    sh2 = zeros(eltype(H2), (axes(H2)..., Base.OneTo(2)))

    E = nothing
    B = nothing
    ϕ = nothing
    ψ = nothing
    lorenz = nothing

    return ControlVolume2D3F{
        typeof(x),
        typeof(w),
        typeof(sw),
        typeof(h0),
        typeof(sh2),
        typeof(E),
        typeof(ϕ),
        typeof(lorenz),
    }(
        x,
        dx,
        y,
        dy,
        w,
        prim,
        sw,
        h0,
        h1,
        h2,
        sh0,
        sh1,
        sh2,
        E,
        B,
        ϕ,
        ψ,
        lorenz,
    )
end

function Base.show(
    io::IO,
    ctr::ControlVolume2D3F{T1,T2,T3,T4,T5,T6,T7,T8},
) where {T1,T2,T3,T4,T5,T6,T7,T8}
    print(
        io,
        "ControlVolume2D3F{$T1,$T2,$T3,$T4,$T5,$T6,$T7,$T8}\n",
        "center: ($(ctr.x),$(ctr.y))\n",
        "interval: ($(ctr.dx),$(ctr.dy))\n",
        "conservative vars: $(ctr.w)\n",
        "primitive vars: $(ctr.prim)\n",
        "conservative slopes: $(ctr.sw)\n",
        "pdf vars: h0, h1, h2\n",
        "pdf slopes: sh0, sh1, sh2\n",
        "electric field: $(ctr.E)\n",
        "magnetic field: $(ctr.B)\n",
        "Lorenz force: $(ctr.lorenz)\n",
    )
end


# ------------------------------------------------------------
# Unstructured cell
# ------------------------------------------------------------

"""
    mutable struct ControlVolumeUS{F,A,B,C} <: AbstractControlVolume
        x::F
        dx::F

        w::A
        prim::B
        sw::C
    end

Unstructured control volume with no distribution function

"""
mutable struct ControlVolumeUS{F,A,B,C} <: AbstractControlVolume
    x::F
    dx::F

    w::A
    prim::B
    sw::C
end

function ControlVolumeUS(X::T, DX::T, W, PRIM) where {T<:Union{Real,AbstractVector}}
    x = deepcopy(X)
    dx = deepcopy(DX)

    w = deepcopy(W)
    prim = deepcopy(PRIM)
    sw = zeros(eltype(W), axes(W)..., axes(X)...)

    return ControlVolumeUS{typeof(x),typeof(w),typeof(prim),typeof(sw)}(x, dx, w, prim, sw)
end


"""
    mutable struct ControlVolumeUS1F{F,A,B,C,D} <: AbstractControlVolume
        x::F
        dx::F

        w::A
        prim::A
        sw::B

        f::C
        sf::D
    end

Unstructured control volume with 1 distribution function

"""
mutable struct ControlVolumeUS1F{F,A,B,C,D} <: AbstractControlVolume
    x::F
    dx::F

    w::A
    prim::A
    sw::B

    f::C
    sf::D
end

function ControlVolumeUS1F(
    X,
    DX,
    W,
    PRIM,
    F::T,
) where {T<:AbstractArray}
    x = deepcopy(X)
    dx = deepcopy(DX)

    w = deepcopy(W)
    prim = deepcopy(PRIM)
    sw = zeros(eltype(W), axes(W)..., axes(X)...)

    f = deepcopy(F)
    sf = zeros(eltype(F), axes(F)..., axes(X)...)

    return ControlVolumeUS1F{typeof(x),typeof(w),typeof(sw),typeof(f),typeof(sf)}(x, dx, w, prim, sw, f, sf)
end



"""
    mutable struct ControlVolumeUS2F{F,A,B,C,D} <: AbstractControlVolume
        x::F
        dx::F

        w::A
        prim::A
        sw::B

        h::C
        b::C
        sh::D
        sb::D
    end

Unstructured control volume with 2 distribution functions

"""
mutable struct ControlVolumeUS2F{F,A,B,C,D} <: AbstractControlVolume
    x::F
    dx::F

    w::A
    prim::A
    sw::B

    h::C
    b::C
    sh::D
    sb::D
end

function ControlVolumeUS2F(
    X,
    DX,
    W::T1,
    PRIM::T1,
    H::T2,
    B::T2,
) where {T1<:AbstractArray,T2<:AbstractArray}
    x = deepcopy(X)
    dx = deepcopy(DX)

    w = deepcopy(W)
    prim = deepcopy(PRIM)
    sw = zeros(eltype(W), axes(W)..., axes(X)...)

    h = deepcopy(H)
    b = deepcopy(B)
    sh = zeros(eltype(H), axes(H)..., axes(X)...)
    sb = zeros(eltype(B), axes(B)..., axes(X)...)

    return ControlVolumeUS2F{typeof(x),typeof(w),typeof(sw),typeof(h),typeof(sh)}(
        x,
        dx,
        w,
        prim,
        sw,
        h,
        b,
        sh,
        sb,
    )
end
