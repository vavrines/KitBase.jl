function ylmKCoefficient(l::Integer, m::Integer)
    k = BigInt(4) # BigInt avoids overflow
    for i = l-m+1:l+m
        k *= i
    end

    return sqrt((2.0 * l + 1) / k / π)
end


function ylmCosSinPolynomial(m::Integer, x::Variable, y::Variable)
    sum = 0.0 * (x + y)
    for j = 0:m÷2
        sum += (-1)^j * binomial(m, 2 * j) * (y^(2 * j)) * (x^(m - 2 * j))
    end

    return sum
end


function ylmSinSinPolynomial(m::Integer, x::Variable, y::Variable)
    sum = 0.0 * (x + y)
    for j = 0:(m-1)÷2
        sum += ((-1)^j) * binomial(m, 2 * j + 1) * (y^(2 * j + 1)) * (x^(m - 2 * j - 1))
    end

    return sum
end


"""
$(SIGNATURES)

Calculation of the spherical harmonic for a given order (l,m) in Cartesian coordinates

## Arguments
- `l`: degree of spherical harmonics
- `m`: order of spherical harmonics
- `x, y, z`: Cartesian coordinates

# Outputs
- spherical harmonic polynomial
"""
function ylm(l::Integer, m::Integer, x::Variable, y::Variable, z::Variable)
    if abs(m) > l
        throw("-l <= m <= l expected, but m = $m and l = $l.")
    end

    p = 1.0 * (z^2 - 1)^l + 0.0 * (x + y)

    for i = 1:l+abs(m)
        c = i <= l ? 1 / (2 * i) : 1.0
        p = c * differentiate(p, z)
    end

    if m > 0
        return sqrt(2) * ylmKCoefficient(l, m) * ylmCosSinPolynomial(m, x, y) * p
    elseif m < 0
        return sqrt(2) * ylmKCoefficient(l, abs(m)) * ylmSinSinPolynomial(abs(m), x, y) * p
    else
        return ylmKCoefficient(l, 0) * p
    end
end


"""
$(SIGNATURES)

Calculate `r^l * Ylm(x,y,z)``
"""
function rlylm(l::Integer, m::Integer, x::Variable, y::Variable, z::Variable)
    p = ylm(l, m, x, y, z)
    tout = []

    for t in terms(p)
        deg = degree(monomial(t))
        degR = l - deg
        push!(tout, (x^2 + y^2 + z^2)^Int(degR / 2) * t)
    end

    return polynomial(tout)
end


"""
$(SIGNATURES)

Evaluate spherical harmonics basis at given quadrature points
"""
function eval_spherharmonic(points::AM{T}, L) where {T<:Real}
    ne = (L + 1)^2
    nq = size(points, 1)
    m = zeros(ne, nq)

    @polyvar x y z
    if L == 0
        spe = [ylm(0, 0, x, y, z)]
    elseif L == 1
        spe = [
            ylm(0, 0, x, y, z),
            ylm(1, -1, x, y, z),
            ylm(1, 0, x, y, z),
            ylm(1, 1, x, y, z),
        ]
    elseif L == 2
        spe = [
            ylm(0, 0, x, y, z),
            ylm(1, -1, x, y, z),
            ylm(1, 0, x, y, z),
            ylm(1, 1, x, y, z),
            ylm(2, -2, x, y, z),
            ylm(2, -1, x, y, z),
            ylm(2, 0, x, y, z),
            ylm(2, 1, x, y, z),
            ylm(2, 2, x, y, z),
        ]
    end

    for j in axes(m, 2), i in axes(m, 1)
        m[i, j] = spe[i](x => points[j, 1], y => points[j, 2], z => points[j, 3])
    end

    return m
end
