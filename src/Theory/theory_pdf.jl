"""
$(SIGNATURES)

Calculate slope of particle distribution function `a = a1 + u * a2 + 0.5 * u^2 * a3`
"""
pdf_slope(u, Δ) = Δ / (u + 1e-7)

"""
$(SIGNATURES)
"""
function pdf_slope(prim::AV, sw::AV, inK)
    sl = similar(sw, axes(prim))

    if length(prim) == 3

        sl[3] =
            4.0 * prim[3]^2 / (inK + 1.0) / prim[1] * (
                2.0 * sw[3] - 2.0 * prim[2] * sw[2] +
                sw[1] * (prim[2]^2 - 0.5 * (inK + 1.0) / prim[3])
            )
        sl[2] = 2.0 * prim[3] / prim[1] * (sw[2] - prim[2] * sw[1]) - prim[2] * sl[3]
        sl[1] =
            sw[1] / prim[1] - prim[2] * sl[2] -
            0.5 * (prim[2]^2 + 0.5 * (inK + 1.0) / prim[3]) * sl[3]

    elseif length(prim) == 4

        sl[4] =
            4.0 * prim[4]^2 / (inK + 2.0) / prim[1] * (
                2.0 * sw[4] - 2.0 * prim[2] * sw[2] - 2.0 * prim[3] * sw[3] +
                sw[1] * (prim[2]^2 + prim[3]^2 - 0.5 * (inK + 2.0) / prim[4])
            )
        sl[3] = 2.0 * prim[4] / prim[1] * (sw[3] - prim[3] * sw[1]) - prim[3] * sl[4]
        sl[2] = 2.0 * prim[4] / prim[1] * (sw[2] - prim[2] * sw[1]) - prim[2] * sl[4]
        sl[1] =
            sw[1] / prim[1] - prim[2] * sl[2] - prim[3] * sl[3] -
            0.5 * (prim[2]^2 + prim[3]^2 + 0.5 * (inK + 2.0) / prim[4]) * sl[4]

    elseif length(prim) == 5

        sl[5] =
            4.0 * prim[5]^2 / (inK + 3.0) / prim[1] * (
                2.0 * sw[5] - 2.0 * prim[2] * sw[2] - 2.0 * prim[3] * sw[3] -
                2.0 * prim[4] * sw[4] +
                sw[1] * (prim[2]^2 + prim[3]^2 + prim[4]^2 - 0.5 * (inK + 3.0) / prim[5])
            )
        sl[4] = 2.0 * prim[5] / prim[1] * (sw[4] - prim[4] * sw[1]) - prim[4] * sl[5]
        sl[3] = 2.0 * prim[5] / prim[1] * (sw[3] - prim[3] * sw[1]) - prim[3] * sl[5]
        sl[2] = 2.0 * prim[5] / prim[1] * (sw[2] - prim[2] * sw[1]) - prim[2] * sl[5]
        sl[1] =
            sw[1] / prim[1] - prim[2] * sl[2] - prim[3] * sl[3] - prim[4] * sl[4] -
            0.5 * (prim[2]^2 + prim[3]^2 + prim[4]^2 + 0.5 * (inK + 3.0) / prim[5]) * sl[5]

    end

    return sl
end


"""
$(SIGNATURES)

Calculate slope of multi-component particle distribution function `a = a1 + u * a2 + 0.5 * u^2 * a3`

"""
function mixture_pdf_slope(prim::AM, sw::AM, inK)
    sl = similar(sw, axes(prim))

    for j in axes(sl, 2)
        sl[:, j] .= pdf_slope(prim[:, j], sw[:, j], inK)
    end

    return sl
end


"""
$(SIGNATURES)

Reduced distribution function

## Arguments
* `f`: particle distribution function with full velocity space
* `weights`: quadrature weights with reduced velocity setting
* `dim`: dimension of the reduced distribution function (1 by default)

For 3D -> 1D, the quadrature weights can be obtained from `VSpace2D`.
"""
function reduce_distribution(f::AM, weights::AV, dim = 1)
    @assert dim <= 3

    if dim == 1
        h = similar(f, axes(f, 1))
        for i in eachindex(h)
            h[i] = sum(@. weights * f[i, :])
        end
    elseif dim == 2
        h = similar(f, axes(f, 2))
        for j in eachindex(h)
            h[j] = sum(@. weights * f[:, j])
        end
    end

    return h
end

"""
$(SIGNATURES)
"""
function reduce_distribution(f::AA{T,3}, weights::AM, dim = 1) where {T}
    @assert dim <= 3

    if dim == 1
        h = similar(f, axes(f, 1))
        for i in eachindex(h)
            h[i] = sum(@. weights * f[i, :, :])
        end
    elseif dim == 2
        h = similar(f, axes(f, 2))
        for j in eachindex(h)
            h[j] = sum(@. weights * f[:, j, :])
        end
    elseif dim == 3
        h = similar(f, axes(f, 3))
        for k in eachindex(h)
            h[k] = sum(@. weights * f[:, :, k])
        end
    else
    end

    return h
end

"""
$(SIGNATURES)
"""
function reduce_distribution(
    f::AA{X,3},
    v::Y,
    w::Y,
    weights::AM,
    dim = 1,
) where {X,Y<:AA{T,3}} where {T}

    @assert dim <= 3

    if dim == 1
        h = similar(f, axes(f, 1))
        b = similar(h)
        for i in eachindex(h)
            h[i] = sum(@. weights * f[i, :, :])
            b[i] = sum(@. weights * (v[i, :, :]^2 + w[i, :, :]^2) * f[i, :, :])
        end
    elseif dim == 2
        h = similar(f, axes(f, 2))
        b = similar(h)
        for j in eachindex(h)
            h[j] = sum(@. weights * f[:, j, :])
            b[j] = sum(@. weights * (v[:, j, :]^2 + w[:, j, :]^2) * f[:, j, :])
        end
    elseif dim == 3
        h = similar(f, axes(f, 3))
        b = similar(h)
        for k in eachindex(h)
            h[k] = sum(@. weights * f[:, :, k])
            b[k] = sum(@. weights * (v[:, :, k]^2 + w[:, :, k]^2) * f[:, :, k])
        end
    end

    return h, b

end


"""
$(SIGNATURES)

Recover full distribution function from reduced ones

## Arguments
* `h, b`: reduced particle distribution function with 1D velocity space
* `u`: quadrature nodes in 1D velocity space
* `weights`: quadrature weights in 1D velocity space
* `v, w`: quadrature nodes in the rest velocity space (with 3D setting)

## Outputs
* `f`: particle distribution function with 3D velocity space
"""
function full_distribution(
    h::X,
    b::X,
    u::Y,
    weights::Y,
    v::Z,
    w::Z,
    ρ,
    γ = 5 / 3,
) where {X<:AV,Y<:AV,Z<:AA3}

    @assert length(h) == size(v, 1) throw(
        DimensionMismatch("reduced and full distribution function mismatch"),
    )

    Ei = 0.5 * discrete_moments(b, u, weights, 0)
    λi = 0.5 * ρ / (γ - 1.0) / Ei / 3.0 * 2.0

    f = similar(v)
    for k in axes(f, 3), j in axes(f, 2), i in axes(f, 1)
        f[i, j, k] = h[i] * (λi / π) * exp(-λi * v[i, j, k]^2) * exp(-λi * w[i, j, k]^2)
    end

    return f

end

"""
$(SIGNATURES)
"""
full_distribution(
    h::X,
    b::X,
    u::Y,
    weights::Y,
    v::Z,
    w::Z,
    prim::A,
    γ = 5 / 3,
) where {X<:AV,Y<:AV,Z<:AA3,A<:AV} = full_distribution(h, b, u, weights, v, w, prim[1], γ)


"""
$(SIGNATURES)

Shift distribution function by external force
"""
function shift_pdf!(f::AV, a, du, dt)
    q0 = eachindex(f) |> first # for OffsetArray
    q1 = eachindex(f) |> last

    if a > 0
        shift = Int(floor(a * dt / du)) # only for uniform velocity grid
        for k = q1:-1:q0+shift
            f[k] = f[k-shift]
        end
        for k = q0:shift+q0-1
            f[k] = 0.0
        end

        for k = q0+1:q1
            f[k] += (dt * a - du * shift) * (f[k-1] - f[k]) / du
        end
    else
        shift = Int(floor(-a * dt / du))
        for k = q0:q1-shift
            f[k] = f[k+shift]
        end
        for k = q1-shift+1:q1
            f[k] = 0.0
        end

        for k = q0:q1-1
            f[k] += (dt * a + du * shift) * (f[k] - f[k+1]) / du
        end
    end

    f[q0] = f[q0+1]
    f[q1] = f[q1-1]

    return nothing
end

"""
$(SIGNATURES)

Multi-component gas
"""
function shift_pdf!(f::AM, a::AV, du::AV, dt)
    for j in axes(f, 2)
        _f = @view f[:, j]
        shift_pdf!(_f, a[j], du[j], dt)
    end

    return nothing
end


"""
$(SIGNATURES)

Recover discrete Chapman-Enskog expansion
"""
function chapman_enskog(u::AV, prim::AV, a::AV, A::AV, τ)

    M = maxwellian(u, prim)
    f = @. M * (
        1 -
        τ * (a[1] * u + a[2] * u^2 + 0.5 * a[3] * u^3 + A[1] + A[2] * u + 0.5 * A[3] * u^2)
    )

    return f

end

"""
$(SIGNATURES)
"""
function chapman_enskog(u::AV, prim::AV, sw::AV, K, τ)

    Mu, Mxi, _, _1 = gauss_moments(prim, K)
    a = pdf_slope(prim, sw, K)
    swt = -prim[1] .* moments_conserve_slope(a, Mu, Mxi, 1)
    A = pdf_slope(prim, swt, K)

    return chapman_enskog(u, prim, a, A, τ)

end

"""
$(SIGNATURES)
"""
function chapman_enskog(u::AA{T1}, v::AA{T1}, prim::AV, a::AV, b::AV, A::AV, τ) where {T1}

    M = maxwellian(u, v, prim)
    f = @. M * (
        1.0 - τ * (a[1] * u + a[2] * u^2 + a[3] * u * v + 0.5 * a[4] * u * (u^2 + v^2)) -
        τ * (b[1] * v + b[2] * u * v + b[3] * v^2 + 0.5 * b[4] * v * (u^2 + v^2)) -
        τ * (A[1] + A[2] * u + A[3] * v + 0.5 * A[4] * (u^2 + v^2))
    )

    return f

end

"""
$(SIGNATURES)
"""
function chapman_enskog(u::AA{T1}, v::AA{T1}, prim::AV, swx::AV, swy::AV, K, τ) where {T1}

    Mu, Mv, Mxi, _, _1 = gauss_moments(prim, K)
    a = pdf_slope(prim, swx, K)
    b = pdf_slope(prim, swy, K)
    sw =
        -prim[1] .* (
            moments_conserve_slope(a, Mu, Mv, Mxi, 1, 0) .+
            moments_conserve_slope(b, Mu, Mv, Mxi, 0, 1)
        )
    A = pdf_slope(prim, sw, K)

    return chapman_enskog(u, v, prim, a, b, A, τ)

end

"""
$(SIGNATURES)
"""
function chapman_enskog(
    u::AA{T1},
    v::AA{T1},
    w::AA{T1},
    prim::AV,
    a::AV,
    b::AV,
    c::AV,
    A::AV,
    τ,
) where {T1}

    M = maxwellian(u, v, w, prim)
    f = @. M * (
        1.0 -
        τ * (
            a[1] * u +
            a[2] * u^2 +
            a[3] * u * v +
            a[4] * u * w +
            0.5 * a[5] * u * (u^2 + v^2 + w^2)
        ) -
        τ * (
            b[1] * v +
            b[2] * u * v +
            b[3] * v^2 +
            b[4] * w * v +
            0.5 * b[5] * v * (u^2 + v^2 + w^2)
        ) -
        τ * (
            c[1] * w +
            c[2] * u * w +
            c[3] * v * w +
            c[4] * w^2 +
            0.5 * c[5] * w * (u^2 + v^2 + w^2)
        ) -
        τ * (A[1] + A[2] * u + A[3] * v + A[4] * w + 0.5 * A[5] * (u^2 + v^2 + w^2))
    )

    return f

end

"""
$(SIGNATURES)
"""
function chapman_enskog(
    u::AA{T1},
    v::AA{T1},
    w::AA{T1},
    prim::AV,
    swx::AV,
    swy::AV,
    swz::AV,
    K,
    τ,
) where {T1}

    Mu, Mv, Mw, _, _1 = gauss_moments(prim, K)
    a = pdf_slope(prim, swx, K)
    b = pdf_slope(prim, swy, K)
    c = pdf_slope(prim, swz, K)
    sw =
        -prim[1] .* (
            moments_conserve_slope(a, Mu, Mv, Mw, 1, 0, 0) .+
            moments_conserve_slope(b, Mu, Mv, Mw, 0, 1, 0) .+
            moments_conserve_slope(c, Mu, Mv, Mw, 0, 0, 1)
        )
    A = pdf_slope(prim, sw, K)

    return chapman_enskog(u, v, w, prim, a, b, c, A, τ)

end


"""
$(SIGNATURES)

Construct a collision invariant of the Boltzmann equation
"""
function collision_invariant(α::AV, u)
    return @. exp(α[1] + α[2] * u + α[3] * u^2 / 2)
end

function collision_invariant(α::AV, vs::AbstractVelocitySpace1D)
    return collision_invariant(α, vs.u)
end

function collision_invariant(α::AV, u, v)
    return @. exp(α[1] + α[2] * u + α[3] * v + α[4] * (u^2 + v^2) / 2)
end

function collision_invariant(α::AV, vs::AbstractVelocitySpace2D)
    return collision_invariant(α, vs.u, vs.v)
end

function collision_invariant(α::AV, u, v, w)
    return @. exp(α[1] + α[2] * u + α[3] * v + α[4] * w + α[5] * (u^2 + v^2 + w^2) / 2)
end

function collision_invariant(α::AV, vs::AbstractVelocitySpace3D)
    return collision_invariant(α, vs.u, vs.v, vs.w)
end

function collision_invariant(α::AM, vs::AbstractVelocitySpace)
    M = [collision_invariant(α[:, j], vs) for j in axes(α, 2)]
    return hcat(M...)
end


"""
$(SIGNATURES)

Calculate the derivatives of distribution function using Hermite polynomials

## Arguments
* `f`: particle distribution function
* `u`: quadrature nodes in velocity space
* `weights`: quadrature weights in velocity space
* `prim`: primitive variables
* `norder`: order of Hermite polynomials

## Outputs
* `df`: derivatives of distribution function in velocity space
"""
function hermite_derivative(f::AV, u::AV, weights::AV, prim, norder)
    # factorial
    factor = OffsetArray{eltype(f)}(undef, 0:norder)
    factor[0] = 1.0
    for i = 1:norder
        factor[i] = 1.0
        for j = 1:i
            factor[i] *= j
        end
    end

    # temperature
    t = 1.0 / prim[3]

    # scaled velocity
    v = @. (u - prim[2]) / sqrt(t) * sqrt(2.0)
    w = @. exp(-v^2 / 2.0) / sqrt(2.0 * π)

    # Hermite polynomials
    ncx = length(f)
    hp = OffsetArray{eltype(f)}(undef, 0:norder+1, ncx)
    hp[0, :] .= 1.0
    hp[1, :] .= v
    hp[2, :] .= v .^ 2 .- 1.0
    for i = 3:norder+1
        hp[i, :] .= @. v * hp[i-1, :] - (i - 1.0) * hp[i-2, :] # recursive formula
    end

    # moments of f
    fmoments = OffsetArray{eltype(f)}(undef, 0:norder)
    for i = 0:norder
        fmoments[i] = sum(f .* weights .* hp[i, :]) * sqrt(2.0 / t)
    end

    # derivatives of f
    df = zero(f)
    for i = 0:norder
        df .+= fmoments[i] * hp[i+1, :] / factor[i]
    end
    @. df *= -sqrt(2.0 / t) * w

    return df
end

"""
$(SIGNATURES)

Calculate the forcing term of distribution function in the Boltzmann equation
"""
function hermite_force(f::AV, u::AV, weights::AV, prim, norder, a)
    return hermite_derivative(f, u, weights, prim, norder) .* a
end
