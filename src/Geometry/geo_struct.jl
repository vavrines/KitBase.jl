"""
1D physical space with structured mesh

    @consts: x0, x1, nx, x, dx

"""
struct PSpace1D{T<:AbstractArray{Float64,1}} <: AbstractPhysicalSpace

    x0::Float64
    x1::Float64
    nx::Int64
    x::T
    dx::T

    function PSpace1D(
        X0::Float64,
        X1::Float64,
        XNUM::Int64,
        X::AbstractArray{Float64,1},
        DX::AbstractArray{Float64,1},
    )
        new{typeof(X)}(X0, X1, XNUM, X, DX)
    end

    PSpace1D() = PSpace1D(0, 1, 100)
    PSpace1D(X0::Real, X1::Real) = PSpace1D(X0, X1, 100)

    function PSpace1D(
        X0::Real,
        X1::Real,
        XNUM::Int,
        TYPE = "uniform"::AbstractString,
        NG = 0::Int,
    )

        x0 = Float64(X0)
        x1 = Float64(X1)
        nx = XNUM
        δ = (x1 - x0) / nx
        x = OffsetArray{Float64}(undef, 1-NG:nx+NG)
        dx = similar(x)

        if TYPE == "uniform" # uniform mesh
            for i in eachindex(x)
                x[i] = x0 + (i - 0.5) * δ
                dx[i] = δ
            end
        end

        # inner constructor method
        new{typeof(x)}(x0, x1, nx, x, dx)

    end

end # struct


"""
2D Physical space with structured mesh

    @consts: x0, x1, nx, y0, y1, ny, x, y, dx, dy

"""
struct PSpace2D{T<:AbstractArray{Float64,2}} <: AbstractPhysicalSpace

    x0::Float64
    x1::Float64
    nx::Int64
    y0::Float64
    y1::Float64
    ny::Int64
    x::T
    y::T
    dx::T
    dy::T

    function PSpace1D(
        X0::Float64,
        X1::Float64,
        XNUM::Int64,
        Y0::Float64,
        Y1::Float64,
        YNUM::Int64,
        X::AbstractArray{Float64,2},
        Y::AbstractArray{Float64,2},
        DX::AbstractArray{Float64,2},
        DY::AbstractArray{Float64,2},
    )
        new{typeof(X)}(X0, X1, XNUM, Y0, Y1, YNUM, X, Y, DX, DY)
    end

    PSpace2D() = PSpace2D(0, 1, 45, 0, 1, 45)
    PSpace2D(X0::Real, X1::Real, Y0::Real, Y1::Real) = PSpace2D(X0, X1, 45, Y0, Y1, 45)

    function PSpace2D(
        X0::Real,
        X1::Real,
        XNUM::Int,
        Y0::Real,
        Y1::Real,
        YNUM::Int,
        TYPE = "uniform"::String,
        NGX = 0::Int,
        NGY = 0::Int,
    )

        x0 = Float64(X0)
        x1 = Float64(X1)
        nx = XNUM
        δx = (x1 - x0) / nx
        y0 = Float64(Y0)
        y1 = Float64(Y1)
        ny = YNUM
        δy = (y1 - y0) / ny
        x = OffsetArray{Float64}(undef, 1-NGX:nx+NGX, 1-NGY:ny+NGY)
        y = OffsetArray{Float64}(undef, 1-NGX:nx+NGX, 1-NGY:ny+NGY)
        dx = OffsetArray{Float64}(undef, 1-NGX:nx+NGX, 1-NGY:ny+NGY)
        dy = OffsetArray{Float64}(undef, 1-NGX:nx+NGX, 1-NGY:ny+NGY)

        if TYPE == "uniform" # rectangular formula
            for j in axes(x, 2)
                for i in axes(x, 1)
                    x[i, j] = x0 + (i - 0.5) * δx
                    y[i, j] = y0 + (j - 0.5) * δy
                    dx[i, j] = δx
                    dy[i, j] = δy
                end
            end
        end

        # inner constructor method
        new{typeof(x)}(x0, x1, nx, y0, y1, ny, x, y, dx, dy)

    end

end # struct


"""
Generate uniform mesh

`uniform_mesh(x0::Real, xnum::Int, dx::Real)`

"""
function uniform_mesh(x0, xnum::T, dx) where {T<:Int}

    points = zeros(xnum)
    for i = 1:xnum
        points[i] = x0 + (i - 0.5) * dx
    end

    return points

end


"""
Equivalent structured mesh generator as matlab

* 2D: `meshgrid(x::AbstractArray{<:Real,1}, y::AbstractArray{<:Real,1})`
* 3D: `meshgrid(x::AbstractArray{<:Real,1}, y::AbstractArray{<:Real,1}, z::AbstractArray{<:Real,1})`

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

    @args x: center locations of 1D mesh
    @args p: point location
    @args mode: choose uniform / nonuniform formulations

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
