"""
    struct PSpace1D{TR<:Real,TI<:Integer,TA<:AbstractArray{<:Real,1}} <: AbstractPhysicalSpace
        x0::TR
        x1::TR
        nx::TI
        x::TA
        dx::TA
    end

1D physical space with structured mesh

- @consts: x0, x1, nx, x, dx
"""
struct PSpace1D{TR<:Real,TI<:Integer,TA<:AbstractArray{<:Real,1}} <: AbstractPhysicalSpace
    x0::TR
    x1::TR
    nx::TI
    x::TA
    dx::TA
end

function PSpace1D(X0::TR, X1::TR, NX::TI, NG = 0::Integer) where {TR,TI}
    δ = (X1 - X0) / NX
    x = OffsetArray{Float64}(undef, 1-NG:NX+NG)
    dx = similar(x)

    # uniform mesh
    for i in eachindex(x)
        x[i] = X0 + (i - 0.5) * δ
        dx[i] = δ
    end

    return PSpace1D{TR,TI,typeof(x)}(X0, X1, NX, x, dx)
end

PSpace1D() = PSpace1D(0, 1, 100)
PSpace1D(X0::T, X1::T) where {T} = PSpace1D(X0, X1, 100)


"""
    struct PSpace2D{TR<:Real,TI<:Integer,TA<:AbstractArray{<:Real,2}} <: AbstractPhysicalSpace
        x0::TR
        x1::TR
        nx::TI
        y0::TR
        y1::TR
        ny::TI
        x::TA
        y::TA
        dx::TA
        dy::TA
    end

2D Physical space with structured mesh

- @consts: x0, x1, nx, y0, y1, ny, x, y, dx, dy
"""
struct PSpace2D{TR<:Real,TI<:Integer,TA<:AbstractArray{<:Real,2}} <: AbstractPhysicalSpace
    x0::TR
    x1::TR
    nx::TI
    y0::TR
    y1::TR
    ny::TI
    x::TA
    y::TA
    dx::TA
    dy::TA
end

function PSpace2D(
    X0::TR,
    X1::TR,
    NX::TI,
    Y0::TR,
    Y1::TR,
    NY::TI,
    NGX = 0::Integer,
    NGY = 0::Integer,
) where {TR,TI}
    δx = (X1 - X0) / NX
    δy = (Y1 - Y0) / NY
    x = OffsetArray{Float64}(undef, 1-NGX:NX+NGX, 1-NGY:NY+NGY)
    y = OffsetArray{Float64}(undef, 1-NGX:NX+NGX, 1-NGY:NY+NGY)
    dx = OffsetArray{Float64}(undef, 1-NGX:NX+NGX, 1-NGY:NY+NGY)
    dy = OffsetArray{Float64}(undef, 1-NGX:NX+NGX, 1-NGY:NY+NGY)

    for j in axes(x, 2)
        for i in axes(x, 1)
            x[i, j] = X0 + (i - 0.5) * δx
            y[i, j] = Y0 + (j - 0.5) * δy
            dx[i, j] = δx
            dy[i, j] = δy
        end
    end

    return PSpace2D{TR,TI,typeof(x)}(X0, X1, NX, Y0, Y1, NY, x, y, dx, dy)
end

PSpace2D() = PSpace2D(0, 1, 45, 0, 1, 45)
PSpace2D(X0::T, X1::T, Y0::T, Y1::T) where {T} = PSpace2D(X0, X1, 45, Y0, Y1, 45)


"""
    uniform_mesh(x0::Real, xnum::Int, dx::Real)

Generate uniform mesh
"""
function uniform_mesh(x0, xnum::T, dx) where {T<:Int}

    points = zeros(xnum)
    for i = 1:xnum
        points[i] = x0 + (i - 0.5) * dx
    end

    return points

end


"""
    2D: meshgrid(x::AbstractArray{<:Real,1}, y::AbstractArray{<:Real,1})
    3D: meshgrid(x::AbstractArray{<:Real,1}, y::AbstractArray{<:Real,1}, z::AbstractArray{<:Real,1})

Equivalent structured mesh generator as matlab
"""
function meshgrid(x::T, y::T) where {T<:AbstractArray{<:Real,1}}
    X = [i for j in y, i in x]
    Y = [j for j in y, i in x]

    return X, Y
end

function meshgrid(x::T, y::T, z::T) where {T<:AbstractArray{<:Real,1}}
    X = [i for k in z, j in y, i in x]
    Y = [j for k in z, j in y, i in x]
    Z = [k for k in z, j in y, i in x]

    return X, Y, Z
end


"""
    find_idx(x::AbstractArray{<:Real,1}, p::Real; mode = :nonuniform::Symbol)

Find the location index of a point in mesh

- @args x: center locations of 1D mesh
- @args p: point location
- @args mode: choose uniform / nonuniform formulations

"""
function find_idx(
    x::T,
    p::Real;
    mode = :nonuniform::Symbol,
) where {T<:AbstractArray{<:Real,1}}

    if mode == :uniform
        dx = x[2] - x[1]
        return Int(ceil((p - x[1] + 0.5 * dx) / dx)) # point location
    else
        return argmin(abs.(x .- p)) # center location
    end

end


# ------------------------------------------------------------
# Extended Base.show()
# ------------------------------------------------------------
function Base.show(io::IO, ps::PSpace1D{TR,TI,TA}) where {TR,TI,TA}
    print(io, "PhysicalSpace1D{$TR,$TI,$TA}\n",
              "domain: ($(ps.x0),$(ps.x1))\n",
              "resolution: $(ps.nx)\n",
              "ghost: $(1-firstindex(ps.x))\n")
end

function Base.show(io::IO, ps::PSpace2D{TR,TI,TA}) where {TR,TI,TA}
    print(io, "PhysicalSpace2D{$TR,$TI,$TA}\n",
              "domain: ($(ps.x0),$(ps.x1)) × ($(ps.y0),$(ps.y1))\n",
              "resolution: $(ps.nx) × $(ps.nx)\n",
              "ghost in x: $(1-firstindex(ps.x[:, 1]))\n",
              "ghost in y: $(1-firstindex(ps.y[1, :]))\n")
end
