# ============================================================
# Structs of Control Volume
# Solver stores data in arrays of struct (AoS)
# ============================================================

# ------------------------------------------------------------
# General
# ------------------------------------------------------------

"""
$(TYPEDEF)

Control volume with no distribution function

# Fields

$(FIELDS)
"""
mutable struct ControlVolume{T1,T2,ND} <: AbstractControlVolume
    w::T1
    prim::T1
    sw::T2
end

function Base.show(io::IO, ctr::ControlVolume{A,B,N}) where {A,B,N}
    print(
        io,
        "ControlVolume$(N)D{$A,$B}\n",
        "conservative vars: $(ctr.w)\n",
        "primitive vars: $(ctr.prim)\n",
        "conservative slopes: $(ctr.sw)\n",
    )
end


"""
$(TYPEDEF)

Control volume with 1 distribution function

# Fields

$(FIELDS)
"""
struct ControlVolume1F{T1,T2,T3,T4,ND} <: AbstractControlVolume
    w::T1
    prim::T1
    sw::T2
    f::T3
    sf::T4
end

function Base.show(io::IO, ctr::ControlVolume1F{A,B,C,D,N}) where {A,B,C,D,N}
    print(
        io,
        "ControlVolume$(N)D1F{$A,$B,$C,$D}\n",
        "conservative vars: $(ctr.w)\n",
        "primitive vars: $(ctr.prim)\n",
        "conservative slopes: $(ctr.sw)\n",
        "pdf vars: f\n",
        "pdf slopes: sf\n",
    )
end


"""
$(TYPEDEF)

Control volume with 2 distribution functions

# Fields

$(FIELDS)
"""
struct ControlVolume2F{T1,T2,T3,T4,ND} <: AbstractControlVolume
    w::T1
    prim::T1
    sw::T2
    h::T3
    b::T3
    sh::T4
    sb::T4
end

function Base.show(io::IO, ctr::ControlVolume2F{A,B,C,D,N}) where {A,B,C,D,N}
    print(
        io,
        "ControlVolume$(N)D2F{$A,$B,$C,$D}\n",
        "conservative vars: $(ctr.w)\n",
        "primitive vars: $(ctr.prim)\n",
        "conservative slopes: $(ctr.sw)\n",
        "pdf vars: f\n",
        "pdf slopes: sf\n",
    )
end

function ControlVolume(W, PRIM, ND)
    w = deepcopy(W)
    prim = deepcopy(PRIM)
    sw = begin
        if ND == 1
            zero(w)
        elseif ND == 2
            slope_array(w)
        end
    end

    return ControlVolume{typeof(w),typeof(sw),ND}(w, prim, sw)
end

function ControlVolume(W, PRIM, F, ND)
    w = deepcopy(W)
    prim = deepcopy(PRIM)
    f = deepcopy(F)
    sw, sf = begin
        if ND == 1
            zero(w), zero(f)
        elseif ND == 2
            slope_array(w), slope_array(f)
        end
    end

    return ControlVolume1F{typeof(w),typeof(sw),typeof(f),typeof(sf),ND}(w, prim, sw, f, sf)
end

function ControlVolume(W, PRIM, H, B, ND)
    w = deepcopy(W)
    prim = deepcopy(PRIM)
    h = deepcopy(H)
    b = deepcopy(B)
    sw, sh, sb = begin
        if ND == 1
            zero(w), zero(h), zero(b)
        elseif ND == 2
            slope_array(w), slope_array(h), slope_array(b)
        end
    end

    return ControlVolume2F{typeof(w),typeof(sw),typeof(h),typeof(sh),ND}(w, prim, sw, h, b, sh, sb)
end

# ------------------------------------------------------------
# 1D
# ------------------------------------------------------------

"""
$(TYPEDEF)

1D control volume with no distribution function

# Fields

$(FIELDS)
"""
mutable struct ControlVolume1D{A,B} <: AbstractControlVolume1D
    w::A
    prim::B
    sw::A
end

function ControlVolume1D(W, PRIM)
    w = deepcopy(W)
    prim = deepcopy(PRIM)
    sw = zero(W)

    return ControlVolume1D{typeof(w),typeof(prim)}(w, prim, sw)
end

function Base.show(io::IO, ctr::ControlVolume1D{A,B}) where {A,B}
    print(
        io,
        "ControlVolume1D{$A,$B}\n",
        "conservative vars: $(ctr.w)\n",
        "primitive vars: $(ctr.prim)\n",
        "conservative slopes: $(ctr.sw)\n",
    )
end


"""
$(TYPEDEF)

1D control volume with 1 distribution function

# Fields

$(FIELDS)
"""
mutable struct ControlVolume1D1F{A,B} <: AbstractControlVolume1D
    w::A
    prim::A
    sw::A

    f::B
    sf::B
end

function ControlVolume1D1F(W::T, PRIM::T, F) where {T}
    w = deepcopy(W)
    prim = deepcopy(PRIM)
    sw = zero(W)

    f = deepcopy(F)
    sf = zero(f)

    return ControlVolume1D1F{typeof(w),typeof(f)}(w, prim, sw, f, sf)
end

function Base.show(io::IO, ctr::ControlVolume1D1F{A,B}) where {A,B}
    print(
        io,
        "ControlVolume1D1F{$A,$B}\n",
        "conservative vars: $(ctr.w)\n",
        "primitive vars: $(ctr.prim)\n",
        "conservative slopes: $(ctr.sw)\n",
        "pdf vars: f\n",
        "pdf slopes: sf\n",
    )
end


"""
$(TYPEDEF)

1D control volume with 2 distribution functions

# Fields

$(FIELDS)
"""
struct ControlVolume1D2F{A,B} <: AbstractControlVolume1D
    w::A
    prim::A
    sw::A

    h::B
    b::B
    sh::B
    sb::B
end

function ControlVolume1D2F(W::T1, PRIM::T1, H::T2, B::T2) where {T1<:AA,T2<:AA}
    w = deepcopy(W)
    prim = deepcopy(PRIM)
    sw = zero(W)

    h = deepcopy(H)
    b = deepcopy(B)
    sh = zero(h)
    sb = zero(b)

    return ControlVolume1D2F{typeof(w),typeof(h)}(w, prim, sw, h, b, sh, sb)
end

function Base.show(io::IO, ctr::ControlVolume1D2F{A,B}) where {A,B}
    print(
        io,
        "ControlVolume1D2F{$A,$B}\n",
        "conservative vars: $(ctr.w)\n",
        "primitive vars: $(ctr.prim)\n",
        "conservative slopes: $(ctr.sw)\n",
        "pdf vars: h, b\n",
        "pdf slopes: sh, sb\n",
    )
end


"""
$(TYPEDEF)

1D control volume with 3 distribution functions

# Fields

$(FIELDS)
"""
mutable struct ControlVolume1D3F{A,B,C,D,E} <: AbstractControlVolume1D
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
    W::AA{<:Real,2},
    PRIM::AA{<:Real,2},
    H0::AA{<:FN,3},
    H1::AA{<:FN,3},
    H2::AA{<:FN,3},
    E0::AA{<:FN,1},
    B0::AA{<:FN,1},
    L::AA{<:FN,2},
)
    w = deepcopy(W)
    prim = deepcopy(PRIM)
    sw = zero(W)

    h0 = deepcopy(H0)
    h1 = deepcopy(H1)
    h2 = deepcopy(H2)
    sh0 = zero(H0)
    sh1 = zero(H1)
    sh2 = zero(H2)

    E = deepcopy(E0)
    B = deepcopy(B0)
    ϕ = 0.0
    ψ = 0.0
    lorenz = deepcopy(L)

    return ControlVolume1D3F{typeof(w),typeof(h0),typeof(E),typeof(ϕ),typeof(lorenz)}(
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
    W::AA{<:Real,1},
    PRIM::AA{<:Real,1},
    H0::AA{<:FN,1},
    H1::AA{<:FN,1},
    H2::AA{<:FN,1},
)
    w = deepcopy(W)
    prim = deepcopy(PRIM)
    sw = zero(W)

    h0 = deepcopy(H0)
    h1 = deepcopy(H1)
    h2 = deepcopy(H2)
    sh0 = zero(H0)
    sh1 = zero(H1)
    sh2 = zero(H2)

    E = nothing
    B = nothing
    ϕ = nothing
    ψ = nothing
    lorenz = nothing

    return ControlVolume1D3F{typeof(w),typeof(h0),typeof(E),typeof(ϕ),typeof(lorenz)}(
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
    W::AA{<:Real,3},
    PRIM::AA{<:Real,3},
    H0::AA{<:FN,4},
    H1::AA{<:FN,4},
    H2::AA{<:FN,4},
    E0::AA{<:FN,2},
    B0::AA{<:FN,2},
    L::AA{<:FN,3},
)
    w = deepcopy(W)
    prim = deepcopy(PRIM)
    sw = zero(W)

    h0 = deepcopy(H0)
    h1 = deepcopy(H1)
    h2 = deepcopy(H2)
    sh0 = zero(H0)
    sh1 = zero(H1)
    sh2 = zero(H2)

    E = deepcopy(E0)
    B = deepcopy(B0)
    ϕ = zero(E[1, :]) # here is difference
    ψ = zero(B[1, :])
    lorenz = deepcopy(L)

    return ControlVolume1D3F{typeof(w),typeof(h0),typeof(E),typeof(ϕ),typeof(lorenz)}(
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

function Base.show(io::IO, ctr::ControlVolume1D3F{A,B,C,D,E}) where {A,B,C,D,E}
    print(
        io,
        "ControlVolume1D3F{$A,$B,$C,$D,$E}\n",
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
$(TYPEDEF)

1D control volume with 4 distribution functions

# Fields

$(FIELDS)
"""
mutable struct ControlVolume1D4F{A,B,C,D,E} <: AbstractControlVolume1D
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
    W::AA{<:Real,2},
    PRIM::AA{<:Real,2},
    H0::AA{<:FN,2},
    H1::AA{Float64,2},
    H2::AA{Float64,2},
    H3::AA{Float64,2},
    E0::AA{Float64,1},
    B0::AA{Float64,1},
    L::AA{Float64,2},
)

    w = deepcopy(W)
    prim = deepcopy(PRIM)
    sw = zero(W)

    h0 = deepcopy(H0)
    h1 = deepcopy(H1)
    h2 = deepcopy(H2)
    h3 = deepcopy(H3)
    sh0 = deepcopy(H0)
    sh1 = deepcopy(H1)
    sh2 = deepcopy(H2)
    sh3 = deepcopy(H3)

    E = deepcopy(E0)
    B = deepcopy(B0)
    ϕ = 0.0
    ψ = 0.0
    lorenz = deepcopy(L)

    return ControlVolume1D4F{typeof(w),typeof(h0),typeof(E),typeof(ϕ),typeof(lorenz)}(
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
    W::AA{<:Real,3},
    PRIM::AA{<:Real,3},
    H0::AA{<:FN,3},
    H1::AA{Float64,3},
    H2::AA{Float64,3},
    H3::AA{Float64,3},
    E0::AA{Float64,2},
    B0::AA{Float64,2},
    L::AA{Float64,3},
)

    w = deepcopy(W)
    prim = deepcopy(PRIM)
    sw = zero(W)

    h0 = deepcopy(H0)
    h1 = deepcopy(H1)
    h2 = deepcopy(H2)
    h3 = deepcopy(H3)
    sh0 = zero(H0)
    sh1 = zero(H1)
    sh2 = zero(H2)
    sh3 = zero(H3)

    E = deepcopy(E0)
    B = deepcopy(B0)
    ϕ = zero(E[1, :])
    ψ = zero(B[1, :])
    lorenz = deepcopy(L)

    return ControlVolume1D4F{typeof(w),typeof(h0),typeof(E),typeof(ϕ),typeof(lorenz)}(
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

function Base.show(io::IO, ctr::ControlVolume1D4F{A,B,C,D,E}) where {A,B,C,D,E}
    print(
        io,
        "ControlVolume1D4F{$A,$B,$C,$D,$E}\n",
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
$(TYPEDEF)

2D control volume with no distribution function

# Fields

$(FIELDS)
"""
mutable struct ControlVolume2D{A,B} <: AbstractControlVolume2D
    w::A
    prim::A
    sw::B
end

function ControlVolume2D(W, PRIM)
    w = deepcopy(W)
    prim = deepcopy(PRIM)
    #sw = zeros(eltype(W), (axes(W)..., Base.OneTo(2)))
    sw = slope_array(W)

    return ControlVolume2D{typeof(w),typeof(sw)}(w, prim, sw)
end

function Base.show(io::IO, ctr::ControlVolume2D{A,B}) where {A,B}
    print(
        io,
        "ControlVolume2D{$A,$B}\n",
        "conservative vars: $(ctr.w)\n",
        "primitive vars: $(ctr.prim)\n",
        "conservative slopes: $(ctr.sw)\n",
    )
end


"""
$(TYPEDEF)

2D control volume with 1 distribution function

# Fields

$(FIELDS)
"""
mutable struct ControlVolume2D1F{A,B,C,D} <: AbstractControlVolume2D
    w::A
    prim::A
    sw::B

    f::C
    sf::D
end

function ControlVolume2D1F(W, PRIM, F::AA)
    w = deepcopy(W)
    prim = deepcopy(PRIM)
    #sw = zeros(eltype(W), (axes(W)..., Base.OneTo(2)))
    sw = slope_array(W)

    f = deepcopy(F)
    #sf = zeros(eltype(F), (axes(F)..., Base.OneTo(2)))
    sf = slope_array(F)

    return ControlVolume2D1F{typeof(w),typeof(sw),typeof(f),typeof(sf)}(w, prim, sw, f, sf)
end

function Base.show(io::IO, ctr::ControlVolume2D1F{A,B,C,D}) where {A,B,C,D}
    print(
        io,
        "ControlVolume2D1F{$A,$B,$C,$D}\n",
        "conservative vars: $(ctr.w)\n",
        "primitive vars: $(ctr.prim)\n",
        "conservative slopes: $(ctr.sw)\n",
        "pdf vars: f\n",
        "pdf slopes: sf\n",
    )
end


"""
$(TYPEDEF)

2D control volume with 2 distribution functions

# Fields

$(FIELDS)
"""
struct ControlVolume2D2F{A,B,C,D} <: AbstractControlVolume2D
    w::A
    prim::A
    sw::B

    h::C
    b::C
    sh::D
    sb::D
end

function ControlVolume2D2F(W::AA, PRIM::AA, H::AA, B::AA)
    w = deepcopy(W)
    prim = deepcopy(PRIM)
    #sw = zeros(eltype(W), (axes(W)..., Base.OneTo(2)))
    sw = slope_array(W)

    h = deepcopy(H)
    b = deepcopy(B)
    #sh = zeros(eltype(H), (axes(H)..., Base.OneTo(2)))
    #sb = zeros(eltype(B), (axes(B)..., Base.OneTo(2)))
    sh = slope_array(H)
    sb = slope_array(B)

    return ControlVolume2D2F{typeof(w),typeof(sw),typeof(h),typeof(sh)}(
        w,
        prim,
        sw,
        h,
        b,
        sh,
        sb,
    )
end

function Base.show(io::IO, ctr::ControlVolume2D2F{A,B,C,D}) where {A,B,C,D}
    print(
        io,
        "ControlVolume2D2F{$A,$B,$C,$D}\n",
        "conservative vars: $(ctr.w)\n",
        "primitive vars: $(ctr.prim)\n",
        "conservative slopes: $(ctr.sw)\n",
        "pdf vars: h, b\n",
        "pdf slopes: sh, sb\n",
    )
end


"""
$(TYPEDEF)

2D control volume with 3 distribution functions

# Fields

$(FIELDS)
"""
mutable struct ControlVolume2D3F{T2,T3,T4,T5,T6,T7,T8} <: AbstractControlVolume2D
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
function ControlVolume2D3F(W::AA, PRIM::AA, H0::AA, H1::AA, H2::AA, E0::AA, B0::AA, L::AA)
    w = deepcopy(W)
    prim = deepcopy(PRIM)
    #sw = zeros(eltype(W), (axes(W)..., Base.OneTo(2))) # 2D
    sw = slope_array(W)

    h0 = deepcopy(H0)
    h1 = deepcopy(H1)
    h2 = deepcopy(H2)
    #sh0 = zeros(eltype(H0), (axes(H0)..., Base.OneTo(2)))
    #sh1 = zeros(eltype(H1), (axes(H1)..., Base.OneTo(2)))
    #sh2 = zeros(eltype(H2), (axes(H2)..., Base.OneTo(2)))
    sh0 = slope_array(H0)
    sh1 = slope_array(H1)
    sh2 = slope_array(H2)

    E = deepcopy(E0)
    B = deepcopy(B0)
    ϕ = 0.0
    ψ = 0.0
    lorenz = deepcopy(L)

    return ControlVolume2D3F{
        typeof(w),
        typeof(sw),
        typeof(h0),
        typeof(sh2),
        typeof(E),
        typeof(ϕ),
        typeof(lorenz),
    }(
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
    W::AA{<:Real,1},
    PRIM::AA{<:Real,1},
    H0::AA{<:FN,2},
    H1::AA{<:FN,2},
    H2::AA{<:FN,2},
)

    w = deepcopy(W)
    prim = deepcopy(PRIM)
    #sw = zeros(eltype(W), (axes(W)..., Base.OneTo(2))) # 2D
    sw = slope_array(W)

    h0 = deepcopy(H0)
    h1 = deepcopy(H1)
    h2 = deepcopy(H2)
    #sh0 = zeros(eltype(H0), (axes(H0)..., Base.OneTo(2)))
    #sh1 = zeros(eltype(H1), (axes(H1)..., Base.OneTo(2)))
    #sh2 = zeros(eltype(H2), (axes(H2)..., Base.OneTo(2)))
    sh0 = slope_array(H0)
    sh1 = slope_array(H1)
    sh2 = slope_array(H2)

    E = nothing
    B = nothing
    ϕ = nothing
    ψ = nothing
    lorenz = nothing

    return ControlVolume2D3F{
        typeof(w),
        typeof(sw),
        typeof(h0),
        typeof(sh2),
        typeof(E),
        typeof(ϕ),
        typeof(lorenz),
    }(
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
    ctr::ControlVolume2D3F{T1,T2,T3,T4,T5,T6,T7},
) where {T1,T2,T3,T4,T5,T6,T7}
    print(
        io,
        "ControlVolume2D3F{$T1,$T2,$T3,$T4,$T5,$T6,$T7}\n",
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
$(TYPEDEF)

Unstructured control volume with no distribution function

# Fields

$(FIELDS)
"""
mutable struct ControlVolumeUS{E,F,A,B,C} <: AbstractUnstructControlVolume
    n::E
    x::F
    dx::F

    w::A
    prim::B
    sw::C
end

function ControlVolumeUS(N, X::T, DX::T, W, PRIM) where {T<:Union{Real,AV}}
    n = deepcopy(N)
    x = deepcopy(X)
    dx = deepcopy(DX)

    w = deepcopy(W)
    prim = deepcopy(PRIM)
    #sw = zeros(eltype(W), axes(W)..., length(N[1]))
    sw = ifelse(
        length(N[1]) == 2,
        slope_array(W; reduction = true),
        slope_array(W; reduction = false),
    )

    return ControlVolumeUS{typeof(n),typeof(x),typeof(w),typeof(prim),typeof(sw)}(
        n,
        x,
        dx,
        w,
        prim,
        sw,
    )
end


"""
$(TYPEDEF)

Unstructured control volume with 1 distribution function

# Fields

$(FIELDS)
"""
struct ControlVolumeUS1F{E,F,A,B,C,D} <: AbstractUnstructControlVolume
    n::E
    x::F
    dx::F

    w::A
    prim::A
    sw::B

    f::C
    sf::D
end

function ControlVolumeUS1F(N, X, DX, W, PRIM, F::AA{T}) where {T}
    n = deepcopy(N)
    x = deepcopy(X)
    dx = deepcopy(DX)

    w = deepcopy(W)
    prim = deepcopy(PRIM)
    #sw = zeros(eltype(W), axes(W)..., length(N[1]))
    sw = ifelse(
        length(N[1]) == 2,
        slope_array(W; reduction = true),
        slope_array(W; reduction = false),
    )

    f = deepcopy(F)
    #sf = zeros(eltype(F), axes(F)..., length(N[1]))
    sf = ifelse(
        length(N[1]) == 2,
        slope_array(F; reduction = true),
        slope_array(F; reduction = false),
    )

    return ControlVolumeUS1F{typeof(n),typeof(x),typeof(w),typeof(sw),typeof(f),typeof(sf)}(
        n,
        x,
        dx,
        w,
        prim,
        sw,
        f,
        sf,
    )
end


"""
$(TYPEDEF)

Unstructured control volume with 2 distribution functions

# Fields

$(FIELDS)
"""
struct ControlVolumeUS2F{E,F,A,B,C,D} <: AbstractUnstructControlVolume
    n::E
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

function ControlVolumeUS2F(N, X, DX, W::T1, PRIM::T1, H::T2, B::T2) where {T1<:AA,T2<:AA}
    n = deepcopy(N)
    x = deepcopy(X)
    dx = deepcopy(DX)

    w = deepcopy(W)
    prim = deepcopy(PRIM)
    #sw = zeros(eltype(W), axes(W)..., length(N[1]))
    sw = ifelse(
        length(N[1]) == 2,
        slope_array(W; reduction = true),
        slope_array(W; reduction = false),
    )

    h = deepcopy(H)
    b = deepcopy(B)
    #sh = zeros(eltype(H), axes(H)..., length(N[1]))
    #sb = zeros(eltype(B), axes(B)..., length(N[1]))
    sh = ifelse(
        length(N[1]) == 2,
        slope_array(H; reduction = true),
        slope_array(H; reduction = false),
    )
    sb = ifelse(
        length(N[1]) == 2,
        slope_array(B; reduction = true),
        slope_array(B; reduction = false),
    )

    return ControlVolumeUS2F{typeof(n),typeof(x),typeof(w),typeof(sw),typeof(h),typeof(sh)}(
        n,
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
