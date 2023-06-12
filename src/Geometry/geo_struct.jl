"""
$(TYPEDEF)

1D physical space with structured mesh

## Fields

$(FIELDS)
"""
struct PSpace1D{TR<:Real,TI<:Integer,TA<:AA} <: AbstractPhysicalSpace1D
    x0::TR
    x1::TR
    nx::TI
    x::TA
    dx::TA
end

function PSpace1D(X0::TR, X1::TR, NX::TI, NG = 0::Integer) where {TR,TI}
    TX = ifelse(TR == Float32, Float32, Float64)

    δ = (X1 - X0) / NX
    x = OffsetArray{TX}(undef, 1-NG:NX+NG)
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
$(TYPEDEF)

2D physical space with structured mesh

## Fields

$(FIELDS)
"""
struct PSpace2D{TR<:Real,TI<:Integer,TA<:AM{<:Real},TB<:AA{<:Real,4},TC,TD} <:
       AbstractPhysicalSpace2D
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
    vertices::TB
    areas::TC
    n::TD
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
    TX = ifelse(TR == Float32, Float32, Float64)

    δx = (X1 - X0) / NX
    δy = (Y1 - Y0) / NY
    x = OffsetArray{TX}(undef, 1-NGX:NX+NGX, 1-NGY:NY+NGY)
    y = similar(x)
    dx = similar(x)
    dy = similar(x)
    for j in axes(x, 2)
        for i in axes(x, 1)
            x[i, j] = X0 + (i - 0.5) * δx
            y[i, j] = Y0 + (j - 0.5) * δy
            dx[i, j] = δx
            dy[i, j] = δy
        end
    end

    vertices = similar(x, axes(x)..., 4, 2)
    for j in axes(vertices, 2), i in axes(vertices, 1)
        vertices[i, j, 1, 1] = x[i, j] - 0.5 * dx[i, j]
        vertices[i, j, 2, 1] = x[i, j] + 0.5 * dx[i, j]
        vertices[i, j, 3, 1] = x[i, j] + 0.5 * dx[i, j]
        vertices[i, j, 4, 1] = x[i, j] - 0.5 * dx[i, j]

        vertices[i, j, 1, 2] = y[i, j] - 0.5 * dy[i, j]
        vertices[i, j, 2, 2] = y[i, j] - 0.5 * dy[i, j]
        vertices[i, j, 3, 2] = y[i, j] + 0.5 * dy[i, j]
        vertices[i, j, 4, 2] = y[i, j] + 0.5 * dy[i, j]
    end

    areas = [0.0 for i in axes(x, 1), j in axes(x, 2), k = 1:4]
    for j in axes(x, 2), i in axes(x, 1)
        areas[i, j, 1] = point_distance(vertices[i, j, 1, :], vertices[i, j, 2, :])
        areas[i, j, 2] = point_distance(vertices[i, j, 2, :], vertices[i, j, 3, :])
        areas[i, j, 3] = point_distance(vertices[i, j, 3, :], vertices[i, j, 4, :])
        areas[i, j, 4] = point_distance(vertices[i, j, 4, :], vertices[i, j, 1, :])
    end

    n = [zeros(2) for i in axes(x, 1), j in axes(x, 2), k = 1:4]
    for j in axes(x, 2), i in axes(x, 1)
        n1 = unit_normal(vertices[i, j, 1, :], vertices[i, j, 2, :])
        n1 .= ifelse(
            dot(
                n1,
                [x[i, j], y[i, j]] .- (vertices[i, j, 1, :] .+ vertices[i, j, 2, :]) ./ 2,
            ) <= 0,
            n1,
            -n1,
        )

        n2 = unit_normal(vertices[i, j, 2, :], vertices[i, j, 3, :])
        n2 .= ifelse(
            dot(
                n2,
                [x[i, j], y[i, j]] .- (vertices[i, j, 2, :] .+ vertices[i, j, 3, :]) ./ 2,
            ) <= 0,
            n2,
            -n2,
        )

        n3 = unit_normal(vertices[i, j, 3, :], vertices[i, j, 4, :])
        n3 .= ifelse(
            dot(
                n3,
                [x[i, j], y[i, j]] .- (vertices[i, j, 3, :] .+ vertices[i, j, 4, :]) ./ 2,
            ) <= 0,
            n3,
            -n3,
        )

        n4 = unit_normal(vertices[i, j, 4, :], vertices[i, j, 1, :])
        n4 .= ifelse(
            dot(
                n4,
                [x[i, j], y[i, j]] .- (vertices[i, j, 4, :] .+ vertices[i, j, 1, :]) ./ 2,
            ) <= 0,
            n4,
            -n4,
        )

        n[i, j, 1] .= n1
        n[i, j, 2] .= n2
        n[i, j, 3] .= n3
        n[i, j, 4] .= n4
    end

    return PSpace2D{TR,TI,typeof(x),typeof(vertices),typeof(areas),typeof(n)}(
        X0,
        X1,
        NX,
        Y0,
        Y1,
        NY,
        x,
        y,
        dx,
        dy,
        vertices,
        areas,
        n,
    )
end

PSpace2D() = PSpace2D(0, 1, 45, 0, 1, 45)

PSpace2D(X0::T, X1::T, Y0::T, Y1::T) where {T} = PSpace2D(X0, X1, 45, Y0, Y1, 45)


"""
$(TYPEDEF)

2D Circular space in polar coordinates

## Fields

$(FIELDS)
"""
struct CSpace2D{TR<:Real,TI<:Integer,TA<:AM{<:Real},TB<:AA{<:Real,4},TC,TD} <:
       AbstractPhysicalSpace2D
    r0::TR
    r1::TR
    nr::TI
    θ0::TR
    θ1::TR
    nθ::TI
    r::TA
    θ::TA
    x::TA
    y::TA
    dr::TA
    dθ::TA
    darc::TA
    vertices::TB
    areas::TC
    n::TD
end

function CSpace2D(
    R0::TR,
    R1::TR,
    NR::TI,
    θ0,
    θ1,
    Nθ::TI,
    NGR = 0::TI,
    NGθ = 0::TI,
) where {TR,TI<:Integer}
    TX = ifelse(TR == Float32, Float32, Float64)

    δr = (R1 - R0) / NR
    δθ = (θ1 - θ0) / Nθ
    r = OffsetArray{TX}(undef, 1-NGR:NR+NGR, 1-NGθ:Nθ+NGθ)
    θ = similar(r)
    x = similar(r)
    y = similar(r)
    dr = similar(r)
    dθ = similar(r)
    darc = similar(r)

    for j in axes(r, 2)
        for i in axes(r, 1)
            r[i, j] = R0 + (i - 0.5) * δr
            θ[i, j] = θ0 + (j - 0.5) * δθ
            x[i, j] = r[i, j] * cos(θ[i, j])
            y[i, j] = r[i, j] * sin(θ[i, j])
            dr[i, j] = δr
            dθ[i, j] = δθ
            darc[i, j] = r[i, j] * dθ[i, j]
        end
    end

    vertices = similar(r, axes(r)..., 4, 2)
    for j in axes(vertices, 2), i in axes(vertices, 1)
        vertices[i, j, 1, 1] = (r[i, j] - 0.5 * dr[i, j]) * cos(θ[i, j] - 0.5 * dθ[i, j])
        vertices[i, j, 2, 1] = (r[i, j] + 0.5 * dr[i, j]) * cos(θ[i, j] - 0.5 * dθ[i, j])
        vertices[i, j, 3, 1] = (r[i, j] + 0.5 * dr[i, j]) * cos(θ[i, j] + 0.5 * dθ[i, j])
        vertices[i, j, 4, 1] = (r[i, j] - 0.5 * dr[i, j]) * cos(θ[i, j] + 0.5 * dθ[i, j])

        vertices[i, j, 1, 2] = (r[i, j] - 0.5 * dr[i, j]) * sin(θ[i, j] - 0.5 * dθ[i, j])
        vertices[i, j, 2, 2] = (r[i, j] + 0.5 * dr[i, j]) * sin(θ[i, j] - 0.5 * dθ[i, j])
        vertices[i, j, 3, 2] = (r[i, j] + 0.5 * dr[i, j]) * sin(θ[i, j] + 0.5 * dθ[i, j])
        vertices[i, j, 4, 2] = (r[i, j] - 0.5 * dr[i, j]) * sin(θ[i, j] + 0.5 * dθ[i, j])
    end

    areas = [0.0 for i in axes(x, 1), j in axes(x, 2), k = 1:4]
    for j in axes(x, 2), i in axes(x, 1)
        areas[i, j, 1] = point_distance(vertices[i, j, 1, :], vertices[i, j, 2, :])
        areas[i, j, 2] = point_distance(vertices[i, j, 2, :], vertices[i, j, 3, :])
        areas[i, j, 3] = point_distance(vertices[i, j, 3, :], vertices[i, j, 4, :])
        areas[i, j, 4] = point_distance(vertices[i, j, 4, :], vertices[i, j, 1, :])
    end

    n = [zeros(2) for i in axes(x, 1), j in axes(x, 2), k = 1:4]
    for j in axes(x, 2), i in axes(x, 1)
        n1 = unit_normal(vertices[i, j, 1, :], vertices[i, j, 2, :])
        n1 .= ifelse(
            dot(
                n1,
                [x[i, j], y[i, j]] .- (vertices[i, j, 1, :] .+ vertices[i, j, 2, :]) ./ 2,
            ) <= 0,
            n1,
            -n1,
        )

        n2 = unit_normal(vertices[i, j, 2, :], vertices[i, j, 3, :])
        n2 .= ifelse(
            dot(
                n2,
                [x[i, j], y[i, j]] .- (vertices[i, j, 2, :] .+ vertices[i, j, 3, :]) ./ 2,
            ) <= 0,
            n2,
            -n2,
        )

        n3 = unit_normal(vertices[i, j, 3, :], vertices[i, j, 4, :])
        n3 .= ifelse(
            dot(
                n3,
                [x[i, j], y[i, j]] .- (vertices[i, j, 3, :] .+ vertices[i, j, 4, :]) ./ 2,
            ) <= 0,
            n3,
            -n3,
        )

        n4 = unit_normal(vertices[i, j, 4, :], vertices[i, j, 1, :])
        n4 .= ifelse(
            dot(
                n4,
                [x[i, j], y[i, j]] .- (vertices[i, j, 4, :] .+ vertices[i, j, 1, :]) ./ 2,
            ) <= 0,
            n4,
            -n4,
        )

        n[i, j, 1] .= n1
        n[i, j, 2] .= n2
        n[i, j, 3] .= n3
        n[i, j, 4] .= n4
    end

    return CSpace2D{TR,TI,typeof(x),typeof(vertices),typeof(areas),typeof(n)}(
        R0,
        R1,
        NR,
        TX(θ0),
        TX(θ1),
        Nθ,
        r,
        θ,
        x,
        y,
        dr,
        dθ,
        darc,
        vertices,
        areas,
        n,
    )
end


"""
$(TYPEDEF)

3D physical space with structured mesh

  z  
   |  
   --- y  
  /  
x  
     8 ----- 7  
     /     /|  
    /     / |  
  5 ----- 6  
   |     |  |  
   |     | / 4  
   |     |/  
  2 ----- 3

## Fields

$(FIELDS)
"""
struct PSpace3D{TR<:Real,TI<:Integer,TA<:AA{<:Real,3},TB<:AA{<:Real,5},TC,TD} <:
       AbstractPhysicalSpace3D
    x0::TR
    x1::TR
    nx::TI
    y0::TR
    y1::TR
    ny::TI
    z0::TR
    z1::TR
    nz::TI
    x::TA
    y::TA
    z::TA
    dx::TA
    dy::TA
    dz::TA
    vertices::TB
    areas::TC
    n::TD
end


"""
$(SIGNATURES)

Generate uniform mesh
"""
function uniform_mesh(x0, xnum::Integer, dx)
    points = zeros(xnum)
    for i = 1:xnum
        points[i] = x0 + (i - 0.5) * dx
    end

    return points
end


"""
$(SIGNATURES)

Equivalent N-dimensional mesh generator as matlab
"""
ndgrid(v::AV) = copy(v)

"""
$(SIGNATURES)
"""
function ndgrid(v1::AV, v2::AV)
    m, n = length(v1), length(v2)
    v1 = reshape(v1, m, 1)
    v2 = reshape(v2, 1, n)

    return (repeat(v1, 1, n), repeat(v2, m, 1))
end

"""
$(SIGNATURES)
"""
function ndgrid(vs::AV{T}...) where {T}
    ndgrid_fill(a, v, s, snext) = begin
        for j = 1:length(a)
            a[j] = v[div(rem(j - 1, snext), s)+1]
        end
    end

    n = length(vs)
    sz = map(length, vs)
    out = ntuple(i -> Array{T}(undef, sz), n)
    s = 1
    for i = 1:n
        a = out[i]::Array
        v = vs[i]
        snext = s * size(a, i)
        ndgrid_fill(a, v, s, snext)
        s = snext
    end

    return out
end


"""
$(SIGNATURES)

Equivalent structured mesh generator as matlab...
"""
meshgrid(v::AV) = meshgrid(v, v)

"""
$(SIGNATURES)
"""
function meshgrid(x::AV, y::AV)
    X = [i for j in y, i in x]
    Y = [j for j in y, i in x]

    return X, Y
end

"""
$(SIGNATURES)
"""
function meshgrid(x::AV, y::AV, z::AV)
    X = [i for k in z, j in y, i in x]
    Y = [j for k in z, j in y, i in x]
    Z = [k for k in z, j in y, i in x]

    return X, Y, Z
end


"""
$(SIGNATURES)

Find the location index of a point in mesh
"""
function find_idx(x::AV, p; mode = :nonuniform::Symbol)
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
    print(
        io,
        "PhysicalSpace1D{$TR,$TI,$TA}\n",
        "domain: ($(ps.x0),$(ps.x1))\n",
        "resolution: $(ps.nx)\n",
        "ghost: $(1-firstindex(ps.x))\n",
    )
end

function Base.show(io::IO, ps::PSpace2D{TR,TI,TA}) where {TR,TI,TA}
    print(
        io,
        "PhysicalSpace2D{$TR,$TI,$TA}\n",
        "domain: ($(ps.x0),$(ps.x1)) × ($(ps.y0),$(ps.y1))\n",
        "resolution: $(ps.nx) × $(ps.ny)\n",
        "ghost in x: $(1-firstindex(ps.x[:, 1]))\n",
        "ghost in y: $(1-firstindex(ps.y[1, :]))\n",
    )
end

function Base.show(io::IO, ps::CSpace2D{TR,TI,TA}) where {TR,TI,TA}
    print(
        io,
        "CircularSpace2D{$TR,$TI,$TA}\n",
        "domain: ($(ps.r0),$(ps.r1)) × ($(ps.θ0),$(ps.θ1))\n",
        "resolution: $(ps.nr) × $(ps.nθ)\n",
        "ghost in r: $(1-firstindex(ps.x[:, 1]))\n",
        "ghost in θ: $(1-firstindex(ps.y[1, :]))\n",
    )
end
