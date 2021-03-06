"""
    pdf_slope(u, Δ)

    pdf_slope(prim::A, sw::B, inK) where {A<:AbstractArray{<:Real,1},B<:AbstractArray{<:Real,1}}

Calculate slope of particle distribution function `a = a1 + u * a2 + 0.5 * u^2 * a3`
"""
pdf_slope(u, Δ::T) where {T<:Real} = Δ / (u + 1e-7)

function pdf_slope(
    prim::X,
    sw::Y,
    inK,
) where {X<:AbstractArray{<:Real,1},Y<:AbstractArray{<:Real,1}}

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
    mixture_pdf_slope(prim::X, sw::Y, inK) where {X<:AbstractArray{<:Real,2},Y<:AbstractArray{<:Real,2}}

Calculate slope of multi-component particle distribution function `a = a1 + u * a2 + 0.5 * u^2 * a3`

"""
function mixture_pdf_slope(
    prim::X,
    sw::Y,
    inK,
) where {X<:AbstractArray{<:Real,2},Y<:AbstractArray{<:Real,2}}

    sl = similar(sw, axes(prim))

    for j in axes(sl, 2)
        sl[:, j] .= pdf_slope(prim[:, j], sw[:, j], inK)
    end

    return sl

end


"""
    reduce_distribution(
        f::X,
        weights::Y,
        dim = 1,
    ) where {X<:AbstractArray{<:AbstractFloat,2},Y<:AbstractArray{<:AbstractFloat,1}}

    reduce_distribution(
        f::X,
        weights::Y,
        dim = 1,
    ) where {X<:AbstractArray{<:AbstractFloat,3},Y<:AbstractArray{<:AbstractFloat,2}}

Reduced distribution function

* @arg : particle distribution function with full velocity space
* @arg : quadrature weights with reduced velocity setting (v & w by default)

"""
function reduce_distribution(
    f::X,
    weights::Y,
    dim = 1,
) where {X<:AbstractArray{<:AbstractFloat,2},Y<:AbstractArray{<:AbstractFloat,1}}

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
    else
        throw("dimension dismatch")
    end

    return h

end

function reduce_distribution(
    f::X,
    weights::Y,
    dim = 1,
) where {X<:AbstractArray{<:AbstractFloat,3},Y<:AbstractArray{<:AbstractFloat,2}}

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

function reduce_distribution(
    f::X,
    v::Y,
    w::Y,
    weights::Z,
    dim = 1,
) where {
    X<:AbstractArray{<:AbstractFloat,3},
    Y<:AbstractArray{<:AbstractFloat,3},
    Z<:AbstractArray{<:AbstractFloat,2},
}

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
    else
    end

    return h, b

end


"""
    full_distribution(
        h::X,
        b::X,
        u::Y,
        weights::Y,
        v::Z,
        w::Z,
        ρ,
        γ = 5 / 3,
    ) where {
        X<:AbstractArray{<:AbstractFloat,1},
        Y<:AbstractArray{<:AbstractFloat,1},
        Z<:AbstractArray{<:AbstractFloat,3},
    }

    full_distribution(
        h::X,
        b::X,
        u::Y,
        weights::Y,
        v::Z,
        w::Z,
        prim::A,
        γ = 5 / 3,
    ) where {
        X<:AbstractArray{<:AbstractFloat,1},
        Y<:AbstractArray{<:AbstractFloat,1},
        Z<:AbstractArray{<:AbstractFloat,3},
        A<:AbstractArray{<:Real,1},
    }

Recover full distribution function from reduced ones

* @arg h & b : reduced particle distribution function with 1D velocity space
* @arg u : quadrature nodes in 1D velocity space
* @arg weights : quadrature weights in 1D velocity space
* @arg v & w : quadrature nodes in the rest velocity space (with 3D setting)
* @return f : particle distribution function with 3D velocity space

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
) where {
    X<:AbstractArray{<:AbstractFloat,1},
    Y<:AbstractArray{<:AbstractFloat,1},
    Z<:AbstractArray{<:AbstractFloat,3},
}

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

full_distribution(
    h::X,
    b::X,
    u::Y,
    weights::Y,
    v::Z,
    w::Z,
    prim::A,
    γ = 5 / 3,
) where {
    X<:AbstractArray{<:AbstractFloat,1},
    Y<:AbstractArray{<:AbstractFloat,1},
    Z<:AbstractArray{<:AbstractFloat,3},
    A<:AbstractArray{<:Real,1},
} = full_distribution(h, b, u, weights, v, w, prim[1], γ)


"""
    shift_pdf!(
        f::T,
        a,
        du,
        dt,
    ) where {T<:AbstractArray{<:AbstractFloat,1}}

    shift_pdf!(
        f::X,
        a::Y,
        du::Z,
        dt,
    ) where {X<:AbstractArray{<:AbstractFloat,2},Y<:AbstractArray{<:Real,1},Z<:AbstractArray{<:AbstractFloat,1}}

Shift distribution function by external force
"""
function shift_pdf!(f::T, a, du, dt) where {T<:AbstractArray{<:AbstractFloat,1}}

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

#--- multi-component gas ---#
function shift_pdf!(
    f::X,
    a::Y,
    du::Z,
    dt,
) where {
    X<:AbstractArray{<:AbstractFloat,2},
    Y<:AbstractArray{<:Real,1},
    Z<:AbstractArray{<:AbstractFloat,1},
}
    for j in axes(f, 2)
        _f = @view f[:, j]
        shift_pdf!(_f, a[j], du[j], dt)
    end

    return nothing
end


"""
    function chapman_enskog(
        u::AbstractVector{T1},
        prim::AbstractVector{T2},
        a::AbstractVector{T3},
        A::AbstractVector{T3},
        τ::Real,
    ) where {T1<:AbstractFloat,T2<:Real,T3<:Real}

    function chapman_enskog(
        u::AbstractMatrix{T1},
        v::AbstractMatrix{T1},
        prim::AbstractVector{T2},
        a::AbstractVector{T3},
        b::AbstractVector{T3},
        A::AbstractVector{T3},
        τ::Real,
    ) where {T1<:AbstractFloat,T2<:Real,T3<:Real}

Recover discrete Chapman-Enskog expansion
"""
function chapman_enskog(
    u::AbstractVector{T1},
    prim::AbstractVector{T2},
    a::AbstractVector{T3},
    A::AbstractVector{T3},
    τ::Real,
) where {T1<:AbstractFloat,T2<:Real,T3<:Real}
    M = maxwellian(u, prim)
    f = @. M * (
        1 -
        τ * (a[1] * u + a[2] * u^2 + 0.5 * a[3] * u^3 + A[1] + A[2] * u + 0.5 * A[3] * u^2)
    )

    return f
end

function chapman_enskog(
    u::AbstractMatrix{T1},
    v::AbstractMatrix{T1},
    prim::AbstractVector{T2},
    a::AbstractVector{T3},
    b::AbstractVector{T3},
    A::AbstractVector{T3},
    τ::Real,
) where {T1<:AbstractFloat,T2<:Real,T3<:Real}
    M = maxwellian(u, v, prim)
    f = @. M * (
        1.0 - τ * (a[1] * u + a[2] * u^2 + a[3] * u * v + 0.5 * a[4] * u * (u^2 + v^2)) -
        τ * (b[1] * v + b[2] * u * v + b[3] * v^2 + 0.5 * b[4] * v * (u^2 + v^2)) -
        τ * (A[1] + A[2] * u + A[3] * v + 0.5 * A[4] * (u^2 + v^2))
    )

    return f
end
