"""
    linspace(start, stop, n::T) where {T<:Integer}

Python linspace function
"""
linspace(start, stop, n::T) where {T<:Integer} =
    collect(range(start, stop = stop, length = n))


"""
    heaviside(x::Real)

Heaviside step function
"""
heaviside(x::T) where {T<:Real} = ifelse(x >= 0, one(x), zero(x))


"""
    fortsign(x::Real, y::Real)

Fortran sign function
"""
fortsign(x::T, y) where {T<:Real} = abs(x) * sign(y)


"""
    mat_split(m::AbstractArray)

Split matrix into row vectors
"""
function mat_split(m::T) where {T<:AbstractArray{<:Number,2}}

    if length(m[1, :]) == 2
        nx = eltype(m).([1.0 0.0])
        ny = eltype(m).([0.0 1.0])

        return nx * m, ny * m
    elseif length(m[1, :]) == 3
        nx = eltype(m).([1.0 0.0 0.0])
        ny = eltype(m).([0.0 1.0 0.0])
        nz = eltype(m).([0.0 0.0 1.0])

        return nx * m, ny * m, nz * m
    end

end


"""
    central_diff(y::AbstractArray{<:Any,1}, x::AbstractArray{<:Any,1})
    central_diff(y::AbstractArray{<:Any,1}, dx::Any)

Central difference
"""
function central_diff(y::T, x::T) where {T<:AbstractVector}

    dy = zeros(eltype(y), axes(y))

    idx = eachindex(y) |> collect
    i0 = idx[1]
    i1 = idx[end]

    dy[i0] = (y[i0+1] - y[i0]) / (x[i0+1] - x[i0] + 1e-7)
    dy[i1] = (y[i1] - y[i1-1]) / (x[i1] - x[i1-1] + 1e-7)
    for i = i0+1:i1-1
        dy[i] = (y[i+1] - y[i-1]) / (x[i+1] - x[i-1] + 1e-7)
    end

    return dy

end

function central_diff(y::T, dx) where {T<:AbstractVector}
    x = ones(eltype(y), axes(y)) .* dx
    dy = central_diff(y, x)

    return dy
end


"""
    central_diff2(y::AbstractArray{<:Any,1}, x::AbstractArray{<:Any,1})

Central difference for second-order derivatives
"""
function central_diff2(y::T, x::T) where {T<:AbstractVector}
    dy = zeros(eltype(y), axes(y))

    i0 = eachindex(y) |> first
    i1 = eachindex(y) |> last

    for i = i0+1:i1-1
        dy[i] = (y[i+1] - 2.0 * y[i] + y[i-1]) / (x[i+1] - x[i-1])^2
    end

    return dy
end

function central_diff2(y::T, dx) where {T<:AbstractVector}
    x = ones(eltype(y), axes(y)) .* dx
    dy = central_diff2(y, x)

    return dy
end


"""
    central_diff!(dy::AbstractArray{<:Any,1}, y::AbstractArray{<:Any,1}, x::AbstractArray{<:Any,1})
    central_diff!(dy::AbstractArray{<:Any,1}, y::AbstractArray{<:Any,1}, dx::Any)

Central difference
"""
function central_diff!(dy::T1, y::T2, x::T2) where {T1<:AbstractVector,T2<:AbstractVector}

    @assert axes(dy) == axes(y) == axes(x)

    idx = eachindex(y) |> collect
    i0 = idx[1]
    i1 = idx[end]

    dy[i0] = (y[i0+1] - y[i0]) / (x[i0+1] - x[i0] + 1e-7)
    dy[i1] = (y[i1] - y[i1-1]) / (x[i1] - x[i1-1] + 1e-7)
    for i = i0+1:i1-1
        dy[i] = (y[i+1] - y[i-1]) / (x[i+1] - x[i-1] + 1e-7)
    end

end

function central_diff!(dy::T1, y::T2, dx) where {T1<:AbstractVector,T2<:AbstractVector}
    x = ones(eltype(y), axes(y)) .* dx
    central_diff!(dy, y, x)
end


"""
    central_diff2!(dy::T1, y::T2, x::T2) where {T1<:AbstractVector,T2<:AbstractVector}
    central_diff2!(dy::T1, y::T2, dx) where {T1<:AbstractVector,T2<:AbstractVector}
    
Central difference
"""
function central_diff2!(dy::T1, y::T2, x::T2) where {T1<:AbstractVector,T2<:AbstractVector}

    @assert axes(dy) == axes(y) == axes(x)

    i0 = eachindex(y) |> first
    i1 = eachindex(y) |> last

    dy[i0] = 0.0
    dy[i1] = 0.0
    for i = i0+1:i1-1
        dy[i] = (y[i+1] - 2.0 * y[i] + y[i-1]) / (x[i+1] - x[i-1])^2
    end

end

function central_diff2!(dy::T1, y::T2, dx) where {T1<:AbstractVector,T2<:AbstractVector}
    x = ones(eltype(y), axes(y)) .* dx
    central_diff2!(dy, y, x)
end


"""
    upwind_diff(
        y::AbstractArray{<:Any,1},
        x::AbstractArray{<:Any,1};
        stream = :right::Symbol,
    )

    upwind_diff(y::AbstractArray{<:Any,1}, dx::Any; stream = :right::Symbol)

Upwind difference
"""
function upwind_diff(
    y::AbstractArray{<:Any,1},
    x::AbstractArray{<:Any,1};
    stream = :right::Symbol,
)

    dy = zeros(eltype(y), axes(y))

    idx = eachindex(y) |> collect
    i0 = idx[1]
    i1 = idx[end]

    if stream == :right
        dy[i0] = 0.0
        for i = i0+1:i1
            dy[i] = (y[i] - y[i-1]) / (x[i] - x[i-1] + 1e-7)
        end
    elseif stream == :left
        dy[i1] = 0.0
        for i = i0:i1-1
            dy[i] = (y[i+1] - y[i]) / (x[i+1] - x[i] + 1e-7)
        end
    else
        throw("streaming direction should be :left or :right")
    end

    return dy

end

function upwind_diff(y::AbstractArray{<:Any,1}, dx::Any; stream = :right::Symbol)
    x = ones(eltype(y), axes(y)) .* dx
    dy = upwind_diff(y, x, stream = stream)

    return dy
end


"""
    upwind_diff!(
        dy::AbstractArray{<:Any,1},
        y::AbstractArray{<:Any,1},
        x::AbstractArray{<:Any,1};
        stream = :right::Symbol,
    )

    upwind_diff!(
        dy::AbstractArray{<:Any,1},
        y::AbstractArray{<:Any,1},
        dx::Any;
        stream = :right::Symbol,
    )

Upwind difference
"""
function upwind_diff!(
    dy::AbstractArray{<:Any,1},
    y::AbstractArray{<:Any,1},
    x::AbstractArray{<:Any,1};
    stream = :right::Symbol,
)

    @assert axes(dy) == axes(y) == axes(x)

    idx = eachindex(y) |> collect
    i0 = idx[1]
    i1 = idx[end]

    if stream == :right
        dy[i0] = 0.0
        for i = i0+1:i1
            dy[i] = (y[i] - y[i-1]) / (x[i] - x[i-1] + 1e-7)
        end
    elseif stream == :left
        dy[i1] = 0.0
        for i = i0:i1-1
            dy[i] = (y[i+1] - y[i]) / (x[i+1] - x[i] + 1e-7)
        end
    else
        throw("streaming direction should be :left or :right")
    end

    return dy

end

function upwind_diff!(
    dy::AbstractArray{<:Any,1},
    y::AbstractArray{<:Any,1},
    dx::Any;
    stream = :right::Symbol,
)
    x = ones(eltype(y), axes(y)) .* dx
    upwind_diff!(dy, y, x, stream = stream)
end


"""
    unstruct_diff(u::AbstractArray{<:Any,1}, x::AbstractArray{<:Any,1}, nx::Integer; mode = :central::Symbol)
    unstruct_diff(u::Function, x::AbstractArray{<:Any,2}, nx::Integer, dim::Integer; mode = :central::Symbol)

Finite difference for pseudo-unstructured mesh
"""
function unstruct_diff(
    u::AbstractArray{<:Any,1},
    x::AbstractArray{<:Any,1},
    nx::Integer;
    mode = :central::Symbol,
)
    uu = reshape(u, (nx, :))
    xx = reshape(x, (nx, :))

    dux = similar(xx)
    for i = 1:nx
        if mode == :central
            dux[i, :] .= central_diff(uu[i, :], xx[i, :])
        elseif mode == :upwind
            dux[i, :] .= upwind_diff(uu[i, :], xx[i, :])
        else
            throw("difference mode should be central or upwind")
        end
    end

    return reshape(dux, (1, :))
end

function unstruct_diff(
    u::Function,
    x::AbstractArray{<:Any,1},
    nx::Integer,
    dim::Integer;
    mode = :central::Symbol,
)
    uu = reshape(u.(x), (nx, :))
    xx = reshape(x, (nx, :))
    dux = zeros(eltype(x), axes(xx))

    if dim == 1
        for j in axes(xx, 2)
            if mode == :central
                dux[:, j] .= central_diff(uu[:, j], xx[:, j])
            elseif mode == :upwind
                dux[:, j] .= upwind_diff(uu[:, j], xx[:, j])
            end
        end
    elseif dim == 2
        for i in axes(xx, 1)
            if mode == :central
                dux[i, :] .= central_diff(uu[i, :], xx[i, :])
            elseif mode == :upwind
                dux[i, :] .= upwind_diff(uu[i, :], xx[i, :])
            end
        end
    end

    return dux[:]
end


"""
    lgwt(N::Integer, a::Real, b::Real)

Gauss Legendre integral for fast spectral method

* @args: number of quadrature points N, integral range [a, b]
* @args: quadrature points x & weights w

"""
function lgwt(N::Integer, a::Real, b::Real)
    x = zeros(N)
    w = zeros(N)

    N1 = N
    N2 = N + 1

    y = zeros(N1)
    y0 = zeros(N1)
    Lp = zeros(N1)
    L = zeros(N1, N2)

    # initial guess
    for i = 1:N1
        y[i] =
            cos((2.0 * (i - 1.0) + 1.0) * 4.0 * atan(1.0) / (2.0 * (N - 1.0) + 2.0)) +
            0.27 / N1 *
            sin(4.0 * atan(1.0) * (-1.0 + i * 2.0 / (N1 - 1.0)) * (N - 1.0) / N2)
        y0[i] = 2.0
    end

    # compute the zeros of the N+1 legendre Polynomial
    # using the recursion relation and the Newton method
    while maximum(abs.(y .- y0)) > 0.0000000000001
        L[:, 1] .= 1.0
        L[:, 2] .= y
        for k = 2:N1
            @. L[:, k+1] = ((2.0 * k - 1.0) * y * L[:, k] - (k - 1) * L[:, k-1]) / k
        end
        @. Lp = N2 * (L[:, N1] - y * L[:, N2]) / (1.0 - y^2)
        @. y0 = y
        @. y = y0 - L[:, N2] / Lp
    end

    # linear map from [-1 1] to [a,b]
    @. x = (a * (1.0 - y) + b * (1.0 + y)) / 2.0
    @. w = N2^2 * (b - a) / ((1.0 - y^2) * Lp^2) / N1^2

    return x, w
end


"""
    extract_last(a::AbstractArray, idx::Integer; mode=:view::Symbol)

Extract subarray except the last column
"""
function extract_last(a::AbstractArray, idx::Integer; mode = :view::Symbol)
    if mode == :copy

        if ndims(a) == 2
            sw = a[:, idx]
        elseif ndims(a) == 3
            sw = a[:, :, idx]
        elseif ndims(a) == 4
            sw = a[:, :, :, idx]
        elseif ndims(a) == 5
            sw = a[:, :, :, :, idx]
        end

    elseif mode == :view

        if ndims(a) == 2
            sw = @view a[:, idx]
        elseif ndims(a) == 3
            sw = @view a[:, :, idx]
        elseif ndims(a) == 4
            sw = @view a[:, :, :, idx]
        elseif ndims(a) == 5
            sw = @view a[:, :, :, :, idx]
        end

    else

        throw("Error in extraction mode setup")

    end

    return sw
end
