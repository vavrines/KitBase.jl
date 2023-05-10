"""
$(SIGNATURES)

Python linspace function
"""
linspace(start, stop, n::T) where {T<:Integer} =
    collect(range(start, stop = stop, length = n))


"""
$(SIGNATURES)

Heaviside step function
"""
heaviside(x::T) where {T<:Real} = ifelse(x >= 0, one(x), zero(x))


"""
$(SIGNATURES)

Fortran sign function
"""
fortsign(x, y) = abs(x) * sign(y)


"""
$(SIGNATURES)

Split matrix into row vectors

This function can be used for building physics-informed neural networks.
"""
function mat_split(m::AM{T}) where {T}
    if length(m[:, 1]) == 2
        return msx_2d(m), msy_2d(m)
    elseif length(m[:, 2]) == 3
        return msx_3d(m), msy_3d(m), msz_3d(m)
    else
        throw("arrays are not aligned")
    end
end

msx_2d(m) = eltype(m).([1.0 0.0]) * m

msx_3d(m) = eltype(m).([1.0 0.0 0.0]) * m

msy_2d(m) = eltype(m).([0.0 1.0]) * m

msy_3d(m) = eltype(m).([0.0 1.0 0.0]) * m

msz_3d(m) = eltype(m).([0.0 0.0 1.0]) * m


"""
$(SIGNATURES)

Gauss Legendre integral for fast spectral method

## Arguments
* `N`: number of quadrature points, 
* `a, b`: integral range
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
$(SIGNATURES)

Extract subarray except the last column
"""
function extract_last(a::AA{T}, idx::Integer; mode = :view::Symbol) where {T}
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


"""
$(SIGNATURES)

Calculate rate of convergence
"""
convergence_order(e1, e2, r = 2) = log(e1 / e2) / log(r)


"""
$(SIGNATURES)

Calculate L1 error
"""
L1_error(u::T, ue::T, Δx) where {T<:AA} = sum(abs.(u .- ue) .* Δx)


"""
$(SIGNATURES)

Calculate L2 error
"""
L2_error(u::T, ue::T, Δx) where {T<:AA} = sqrt(sum((abs.(u .- ue) .* Δx) .^ 2))


"""
$(SIGNATURES)

Calculate L∞ error
"""
L∞_error(u::T, ue::T, Δx) where {T<:AA} = maximum(abs.(u .- ue) .* Δx)


"""
$(SIGNATURES)

Discrete Dirac delta function
"""
dirac_delta(r) = dirac_delta(r, Class{1})

function dirac_delta(r, ::Type{Class{1}})
    if abs(r) <= 1
        return 0.125 * (3 - 2 * abs(r) + sqrt(1 + 4 * abs(r) - 4 * r^2))
    elseif 1 < abs(r) <= 2
        return 0.125 * (5 - 2 * abs(r) - sqrt(-7 + 12 * abs(r) - 4 * r^2))
    else
        return 0.0
    end
end

function dirac_delta(r, ::Type{Class{2}})
    if abs(r) <= 2
        return 0.25 * (1 + cos(0.5 * π * abs(r)))
    else
        return 0.0
    end
end

function dirac_delta(r, ::Type{Class{3}})
    if abs(r) <= 0.5
        return 3 / 8 + π / 32 - r^2 / 4
    elseif 0.5 < abs(r) <= 1.5
        return 1 / 4 + (1 - abs(r)) / 8 * sqrt(-2 + 8 * abs(r) - 4 * r^2) -
               asin(sqrt(2) * (abs(r) - 1)) / 8
    elseif 1.5 < abs(r) <= 2.5
        return 17 / 16 - π / 64 - 3 * abs(r) / 4 +
               r^2 / 8 +
               (abs(r) - 2) / 16 * sqrt(-14 + 16 * abs(r) - 4 * r^2) +
               asin(sqrt(2) * (abs(r) - 2)) / 16
    else
        return 0.0
    end
end

dirac_delta(x, y, Δx, Δy) = 1 / (Δx * Δy) * dirac_delta(x / Δx) * dirac_delta(y / Δy)

dirac_delta(x, y, Δx, Δy, T) =
    1 / (Δx * Δy) * dirac_delta(x / Δx, T) * dirac_delta(y / Δy, T)
