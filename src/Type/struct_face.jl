# ============================================================
# Structs of Cell Interface
# compatible with control volume structs
# ============================================================

# ------------------------------------------------------------
# 1D
# ------------------------------------------------------------

"""
    mutable struct Interface1D{A<:Union{Real,AA}} <: AbstractInterface1D
        fw::A
    end

1D cell interface with no distribution function
Note that deepcopy constructor is needed to work with StructArrays

"""
mutable struct Interface1D{A<:Union{Real,AA}} <: AbstractInterface1D
    fw::A
end

function Base.show(io::IO, ctr::Interface1D{A}) where {A}
    print(io, "Interface1D{$A}\n", "conservative fluxes: $(ctr.fw)\n")
end


"""
    mutable struct Interface1D1F{A,B} <: AbstractInterface1D
        fw::A
        ff::B
    end

1D cell interface with 1 distribution function
Note that deepcopy constructor is needed to work with StructArrays

"""
mutable struct Interface1D1F{A,B<:AA} <: AbstractInterface1D
    fw::A
    ff::B
end

function Base.show(io::IO, ctr::Interface1D1F{A,B}) where {A,B}
    print(
        io,
        "Interface1D1F{$A,$B}\n",
        "conservative fluxes: $(ctr.fw)\n",
        "pdf fluxes: ff\n",
    )
end


"""
    mutable struct Interface1D2F{A,B} <: AbstractInterface1D
        fw::A
        fh::B
        fb::B

        function Interface1D2F(w::AA, f::AA)
            fw = zero(w)
            fh = zero(f)
            fb = zero(f)

            new{typeof(fw),typeof(fh)}(fw, fh, fb)
        end
    end

1D cell interface with 2 distribution functions

"""
mutable struct Interface1D2F{A,B<:AA} <: AbstractInterface1D
    fw::A
    fh::B
    fb::B
end

function Interface1D2F(w, f::AA)
    fw = zero(w)
    fh = zero(f)
    fb = zero(f)

    return Interface1D2F(fw, fh, fb)
end

function Base.show(io::IO, ctr::Interface1D2F{A,B}) where {A,B}
    print(
        io,
        "Interface1D2F{$A,$B}\n",
        "conservative fluxes: $(ctr.fw)\n",
        "pdf fluxes: fh, fb\n",
    )
end


"""
    mutable struct Interface1D3F{A,B,C} <: AbstractInterface1D
        fw::A
        fh0::B
        fh1::B
        fh2::B
        femL::C
        femR::C
    end

1D cell interface with 3 distribution functions

"""
mutable struct Interface1D3F{A,B,C} <: AbstractInterface1D
    fw::A
    fh0::B
    fh1::B
    fh2::B
    femL::C
    femR::C
end

# deterministic
function Interface1D3F(w::AA, f::AA, E::AA{<:Real,1})
    fw = zeros(eltype(w), axes(w))
    fh0 = zeros(eltype(f), axes(f))
    fh1 = zeros(eltype(f), axes(f))
    fh2 = zeros(eltype(f), axes(f))
    femL = zeros(eltype(E), 8)
    femR = zeros(eltype(E), 8)

    return Interface1D3F(fw, fh0, fh1, fh2, femL, femR)
end

# stochastic
function Interface1D3F(w::AA, f::AA, E::AA{<:Real,2})
    fw = zeros(eltype(w), axes(w))
    fh0 = zeros(eltype(f), axes(f))
    fh1 = zeros(eltype(f), axes(f))
    fh2 = zeros(eltype(f), axes(f))
    femL = zeros(eltype(E), 8, axes(E, 2))
    femR = zeros(eltype(E), 8, axes(E, 2))

    return Interface1D3F(fw, fh0, fh1, fh2, femL, femR)
end

# Rykov
function Interface1D3F(w::AA, f::AA)
    fw = zero(w)
    fh0 = zero(f)
    fh1 = zero(f)
    fh2 = zero(f)
    femL = nothing
    femR = nothing

    return Interface1D3F(fw, fh0, fh1, fh2, femL, femR)
end

function Base.show(io::IO, ctr::Interface1D3F{A,B,C}) where {A,B,C}
    print(
        io,
        "Interface1D3F{$A,$B,$C}\n",
        "conservative fluxes: $(ctr.fw)\n",
        "pdf fluxes: fh0, fh1, fh2\n",
        "electromagnetic fluxes: femL, femR\n",
    )
end


"""
    mutable struct Interface1D4F{A,B,C} <: AbstractInterface1D
        fw::A
        fh0::B
        fh1::B
        fh2::B
        fh3::B
        femL::C
        femR::C
    end

1D cell interface with 4 distribution functions

"""
mutable struct Interface1D4F{A,B,C} <: AbstractInterface1D
    fw::A
    fh0::B
    fh1::B
    fh2::B
    fh3::B
    femL::C
    femR::C
end

# deterministic
function Interface1D4F(w::AA, f::AA, E::AA{<:Real,1})
    fw = zeros(eltype(w), axes(w))
    fh0 = zeros(eltype(f), axes(f))
    fh1 = zeros(eltype(f), axes(f))
    fh2 = zeros(eltype(f), axes(f))
    fh3 = zeros(eltype(f), axes(f))
    femL = zeros(eltype(E), 8)
    femR = zeros(eltype(E), 8)

    return Interface1D4F(fw, fh0, fh1, fh2, fh3, femL, femR)
end

# stochastic
function Interface1D4F(w::AA, f::AA, E::AA{<:Real,2})
    fw = zeros(eltype(w), axes(w))
    fh0 = zeros(eltype(f), axes(f))
    fh1 = zeros(eltype(f), axes(f))
    fh2 = zeros(eltype(f), axes(f))
    fh3 = zeros(eltype(f), axes(f))
    femL = zeros(eltype(E), 8, axes(E, 2))
    femR = zeros(eltype(E), 8, axes(E, 2))

    return Interface1D4F(fw, fh0, fh1, fh2, fh3, femL, femR)
end

function Base.show(io::IO, ctr::Interface1D4F{A,B,C}) where {A,B,C}
    print(
        io,
        "Interface1D4F{$A,$B,$C}\n",
        "conservative fluxes: $(ctr.fw)\n",
        "pdf fluxes: fh0, fh1, fh2, fh3\n",
        "electromagnetic fluxes: femL, femR\n",
    )
end

# ------------------------------------------------------------
# 2D
# ------------------------------------------------------------

"""
    mutable struct Interface2D{A,B,C} <: AbstractInterface2D
        len::A
        n::B
        fw::C

        function Interface2D(L::Real, C::Real, S::Real, w::AA)
            len = L
            n = [C, S]

            fw = zero(w)

            new{typeof(len),typeof(n),typeof(fw)}(len, n, fw)
        end
    end

2D cell interface with no distribution function

"""
mutable struct Interface2D{A,B,C} <: AbstractInterface2D
    len::A
    n::B
    fw::C
end

function Interface2D(L::Real, C::Real, S::Real, w::AA)
    len = L
    n = [C, S]

    fw = zero(w)

    return Interface2D(len, n, fw)
end

function Base.show(io::IO, ctr::Interface2D{A,B,C}) where {A,B,C}
    print(
        io,
        "Interface2D{$A,$B,$C}\n",
        "length: $(ctr.len)\n",
        "normal vector: ($(ctr.n[1]),$(ctr.n[2]))\n",
        "conservative fluxes: $(ctr.fw)\n",
    )
end


"""
    mutable struct Interface2D1F{A,B,C,D} <: AbstractInterface2D
        len::A
        n::B
        fw::C
        ff::D

        function Interface2D1F(L::Real, C::Real, S::Real, w::AA, f::AA)
            len = L
            n = [C, S]

            fw = zero(w)
            ff = zero(f)

            new{typeof(len),typeof(n),typeof(fw),typeof(ff)}(len, n, fw, ff)
        end
    end

2D cell interface with 1 distribution function

"""
mutable struct Interface2D1F{A,B,C,D} <: AbstractInterface2D
    len::A
    n::B
    fw::C
    ff::D
end

function Interface2D1F(L::Real, C::Real, S::Real, w, f::AA)
    len = L
    n = [C, S]

    fw = zero(w)
    ff = zero(f)

    return Interface2D1F(len, n, fw, ff)
end

function Base.show(io::IO, ctr::Interface2D1F{A,B,C,D}) where {A,B,C,D}
    print(
        io,
        "Interface2D1F{$A,$B,$C}\n",
        "length: $(ctr.len)\n",
        "normal vector: ($(ctr.n[1]),$(ctr.n[2]))\n",
        "conservative fluxes: $(ctr.fw)\n",
        "pdf fluxes: ff",
    )
end


"""
    mutable struct Interface2D2F{A,B,C,D} <: AbstractInterface2D
        len::A
        n::B
        fw::C
        fh::D
        fb::D

        function Interface2D2F(L::Real, C::Real, S::Real, w::AA, f::AA)
            len = L
            n = [C, S]

            fw = zeros(eltype(w), axes(w))
            fh = zeros(eltype(f), axes(f))
            fb = zeros(eltype(f), axes(f))

            new{typeof(len),typeof(n),typeof(fw),typeof(fh)}(len, n, fw, fh, fb)
        end
    end

2D cell interface with 2 distribution functions

"""
mutable struct Interface2D2F{A,B<:AV,C,D} <: AbstractInterface2D
    len::A
    n::B
    fw::C
    fh::D
    fb::D
end

function Interface2D2F(L::Real, C::Real, S::Real, w::AA, f::AA)
    len = L
    n = @SArray [C, S]

    fw = zero(w)
    fh = zero(f)
    fb = zero(f)

    return Interface2D2F(len, n, fw, fh, fb)
end

function Base.show(io::IO, ctr::Interface2D2F{A,B,C,D}) where {A,B,C,D}
    print(
        io,
        "Interface2D2F{$A,$B,$C}\n",
        "length: $(ctr.len)\n",
        "normal vector: ($(ctr.n[1]),$(ctr.n[2]))\n",
        "conservative fluxes: $(ctr.fw)\n",
        "pdf fluxes: fh, fb",
    )
end


"""
    mutable struct Interface2D3F{A,B,C,D,E} <: AbstractInterface2D
        len::A
        n::B
        fw::C
        fh0::D
        fh1::D
        fh2::D
        femL::E
        femR::E
        femLU::E
        femLD::E
        femRU::E
        femRD::E
    end

2D cell interface with 3 distribution functions

"""
mutable struct Interface2D3F{A,B,C,D,E} <: AbstractInterface2D
    len::A
    n::B
    fw::C
    fh0::D
    fh1::D
    fh2::D
    femL::E
    femR::E
    femLU::E
    femLD::E
    femRU::E
    femRD::E
end

function Interface2D3F(L::Real, C::Real, S::Real, w::AA, f::AA, E::AA{<:Real,1})
    len = L
    n = @SArray [C, S]

    fw = zero(w)
    fh0 = zero(f)
    fh1 = zero(f)
    fh2 = zero(f)

    femL = zeros(eltype(E), 8)
    femR = zeros(eltype(E), 8)
    femLU = zeros(eltype(E), 8)
    femLD = zeros(eltype(E), 8)
    femRU = zeros(eltype(E), 8)
    femRD = zeros(eltype(E), 8)

    return Interface2D3F(len, n, fw, fh0, fh1, fh2, femL, femR, femLU, femLD, femRU, femRD)
end

function Interface2D3F(L::Real, C::Real, S::Real, w::AA, f::AA, E::AA{<:Real,2})
    len = L
    n = [C, S]

    fw = zero(w)
    fh0 = zero(f)
    fh1 = zero(f)
    fh2 = zero(f)

    femL = zeros(eltype(E), 8, axes(E, 2))
    femR = zeros(eltype(E), 8, axes(E, 2))
    femLU = zeros(eltype(E), 8, axes(E, 2))
    femLD = zeros(eltype(E), 8, axes(E, 2))
    femRU = zeros(eltype(E), 8, axes(E, 2))
    femRD = zeros(eltype(E), 8, axes(E, 2))

    return Interface2D3F(len, n, fw, fh0, fh1, fh2, femL, femR, femLU, femLD, femRU, femRD)
end

# Rykov
function Interface2D3F(L::Real, C::Real, S::Real, w::AA, f::AA)
    len = L
    n = [C, S]

    fw = zero(w)
    fh0 = zero(f)
    fh1 = zero(f)
    fh2 = zero(f)

    femL = nothing
    femR = nothing
    femLU = nothing
    femLD = nothing
    femRU = nothing
    femRD = nothing

    return Interface2D3F(len, n, fw, fh0, fh1, fh2, femL, femR, femLU, femLD, femRU, femRD)
end

function Base.show(io::IO, ctr::Interface2D3F{A,B,C,D,E}) where {A,B,C,D,E}
    print(
        io,
        "Interface2D3F{$A,$B,$C}\n",
        "length: $(ctr.len)\n",
        "normal vector: ($(ctr.n[1]),$(ctr.n[2]))\n",
        "conservative fluxes: $(ctr.fw)\n",
        "pdf fluxes: fh0, fh1, fh2\n",
        "electromagnetic fluxes: femL, femR, femLU, femLD, femRU, femRD\n",
    )
end
