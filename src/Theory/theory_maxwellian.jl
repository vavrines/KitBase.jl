"""
    maxwellian(u, ρ, U, λ)
    maxwellian(u, prim)
    maxwellian(u, v, ρ, U, V, λ)
    maxwellian(u, v, prim)
    maxwellian(u, v, w, ρ, U, V, W, λ)
    maxwellian(u, v, w, prim)

Maxwellian in discrete form

* @args: particle velocity quadrature points
* @args: density, velocity and inverse of temperature
* @return: Maxwellian distribution function

"""
maxwellian(u, ρ, U, λ) = @. ρ * sqrt(λ / π) * exp(-λ * (u - U)^2) # 1V

maxwellian(u, prim::AbstractVector{T}) where {T<:Real} =
    maxwellian(u, prim[1], prim[2], prim[end]) # in case of input with length 4/5

#--- 2V ---#
maxwellian(u, v, ρ, U, V, λ) = @. ρ * (λ / π) * exp(-λ * ((u - U)^2 + (v - V)^2))

maxwellian(u, v, prim::AbstractVector{T}) where {T<:Real} =
    maxwellian(u, v, prim[1], prim[2], prim[3], prim[end]) # in case of input with length 5

#--- 3V ---#
maxwellian(u, v, w, ρ, U, V, W, λ) =
    @. ρ * sqrt((λ / π)^3) * exp(-λ * ((u - U)^2 + (v - V)^2 + (w - W)^2))

maxwellian(u, v, w, prim::AbstractVector{T}) where {T<:Real} =
    maxwellian(u, v, w, prim[1], prim[2], prim[3], prim[4], prim[5])


"""
    maxwellian!(M, u, ρ, U, λ)
    maxwellian!(M, u, prim)

    # Rykov
    maxwellian!(Ht, Bt, Rt, Hr, Br, Rr, u, prim, K, Kr)

    maxwellian!(M, u, v, ρ, U, V, λ)
    maxwellian!(M, u, v, prim)

    # Rykov
    maxwellian!(Ht, Bt, Rt, Hr, Br, Rr, u, v, prim, K, Kr)

    maxwellian!(M, u, v, ρ, U, V, W, λ)
    maxwellian!(M, u, v, w, prim)

In-place Maxwellian

* @args: particle velocity quadrature points
* @args: density, velocity and inverse of temperature
* @return: Maxwellian distribution function

"""
function maxwellian!(
    M::AbstractVector{T1},
    u::AbstractVector{T2},
    ρ,
    U,
    λ,
) where {T1<:AbstractFloat,T2<:AbstractFloat}
    @. M = ρ * sqrt(λ / π) * exp(-λ * (u - U)^2) # 1V
    return nothing
end

maxwellian!(
    M::AbstractVector{T1},
    u::AbstractVector{T2},
    prim::AbstractVector{T3},
) where {T1<:AbstractFloat,T2<:AbstractFloat,T3<:Real} =
    maxwellian!(M, u, prim[1], prim[2], prim[end])

# Rykov
function maxwellian!(
    Ht::T1,
    Bt::T1,
    Rt::T1,
    Hr::T1,
    Br::T1,
    Rr::T1,
    u::T2,
    prim::T3,
    K,
    Kr,
) where {
    T1<:AbstractVector{<:AbstractFloat},
    T2<:AbstractVector{<:AbstractFloat},
    T3<:AbstractVector{<:Real},
}
    @. Ht = prim[1] * sqrt(prim[4] / π) * exp(-prim[4] * (u - prim[2])^2)
    @. Bt = Ht * K / (2.0 * prim[4])
    @. Rt = Ht * Kr / (2.0 * prim[5])

    @. Hr = prim[1] * sqrt(prim[3] / π) * exp(-prim[3] * (u - prim[2])^2)
    @. Br = Hr * K / (2.0 * prim[3])
    @. Rr = Hr * Kr / (2.0 * prim[3])

    return nothing
end

#--- 2V ---#
function maxwellian!(
    M::T1,
    u::T2,
    v::T2,
    ρ,
    U,
    V,
    λ,
) where {T1<:AbstractArray{<:AbstractFloat,2},T2<:AbstractArray{<:AbstractFloat}}
    @. M = ρ * (λ / π) * exp(-λ * ((u - U)^2 + (v - V)^2))
    return nothing
end

maxwellian!(
    M::T1,
    u::T2,
    v::T2,
    prim::T3,
) where {
    T1<:AbstractArray{<:AbstractFloat,2},
    T2<:AbstractArray{<:AbstractFloat},
    T3<:AbstractArray{<:Real,1},
} = maxwellian!(M, u, v, prim[1], prim[2], prim[3], prim[end])

# Rykov
function maxwellian!(
    Ht::T1,
    Bt::T1,
    Rt::T1,
    Hr::T1,
    Br::T1,
    Rr::T1,
    u::T2,
    v::T2,
    prim::T3,
    K,
    Kr,
) where {
    T1<:AbstractArray{<:AbstractFloat,2},
    T2<:AbstractArray{<:AbstractFloat},
    T3<:AbstractArray{<:Real,1},
}
    @. Ht = prim[1] * (prim[5] / π) * exp(-prim[5] * ((u - prim[2])^2 + (v - prim[3])^2))
    @. Bt = Ht * K / (2.0 * prim[5])
    @. Rt = Ht * Kr / (2.0 * prim[6])

    @. Hr = prim[1] * (prim[4] / π) * exp(-prim[4] * ((u - prim[2])^2 + (v - prim[3])^2))
    @. Br = Hr * K / (2.0 * prim[4])
    @. Rr = Hr * Kr / (2.0 * prim[4])

    return nothing
end

#--- 3V ---#
function maxwellian!(
    M::T1,
    u::T2,
    v::T2,
    w::T2,
    ρ,
    U,
    V,
    W,
    λ,
) where {T1<:AbstractArray{<:AbstractFloat,3},T2<:AbstractArray{<:AbstractFloat}}
    @. M = ρ * sqrt((λ / π)^3) * exp(-λ * ((u - U)^2 + (v - V)^2 + (w - W)^2))
    return nothing
end

maxwellian!(
    M::T1,
    u::T2,
    v::T2,
    w::T2,
    prim::T3,
) where {
    T1<:AbstractArray{<:AbstractFloat,3},
    T2<:AbstractArray{<:AbstractFloat},
    T3<:AbstractArray{<:Real,1},
} = maxwellian!(M, u, v, w, prim[1], prim[2], prim[3], prim[4], prim[5])


"""
    mixture_maxwellian(u::X, prim::Y) where {X<:AbstractArray{<:AbstractFloat,2},Y<:AbstractArray{<:Real,2}}

    mixture_maxwellian(
        u::X,
        v::X,
        prim::Y,
    ) where {X<:AbstractArray{<:AbstractFloat,3},Y<:AbstractArray{<:Real,2}}

    mixture_maxwellian(
        u::X,
        v::X,
        w::X,
        prim::Y,
    ) where {X<:AbstractArray{<:AbstractFloat,4},Y<:AbstractArray{<:Real,2}}

Multi-component Maxwellian in discrete form
"""
function mixture_maxwellian(
    u::X,
    prim::Y,
) where {X<:AbstractArray{<:AbstractFloat,2},Y<:AbstractArray{<:Real,2}}

    mixM = similar(u)
    for j in axes(mixM, 2)
        mixM[:, j] .= maxwellian(u[:, j], prim[:, j])
    end

    return mixM

end

function mixture_maxwellian(
    u::X,
    v::X,
    prim::Y,
) where {X<:AbstractArray{<:AbstractFloat,3},Y<:AbstractArray{<:Real,2}}

    mixM = similar(u)
    for k in axes(mixM, 3)
        mixM[:, :, k] .= maxwellian(u[:, :, k], v[:, :, k], prim[:, k])
    end

    return mixM

end

function mixture_maxwellian(
    u::X,
    v::X,
    w::X,
    prim::Y,
) where {X<:AbstractArray{<:AbstractFloat,4},Y<:AbstractArray{<:Real,2}}

    mixM = similar(u)
    for l in axes(mixM, 4)
        mixM[:, :, :, l] .=
            maxwellian(u[:, :, :, l], v[:, :, :, l], w[:, :, :, l], prim[:, l])
    end

    return mixM

end


"""
    mixture_maxwellian!(
        M::T1,
        u::T2,
        prim::T3,
    ) where {
        T1<:AbstractArray{<:AbstractFloat,2},
        T2<:AbstractArray{<:AbstractFloat,2},
        T3<:AbstractArray{<:Real,2},
    }

    mixture_maxwellian!(
        M::T1,
        u::T2,
        v::T2,
        prim::T3,
    ) where {
        T1<:AbstractArray{<:AbstractFloat,3},
        T2<:AbstractArray{<:AbstractFloat,3},
        T3<:AbstractArray{<:Real,2},
    }

    mixture_maxwellian!(
        M::T1,
        u::T2,
        v::T2,
        w::T2,
        prim::T3,
    ) where {
        T1<:AbstractArray{<:AbstractFloat,4},
        T2<:AbstractArray{<:AbstractFloat,4},
        T3<:AbstractArray{<:Real,2},
    }

In-place multi-component Maxwellian

"""
function mixture_maxwellian!(
    M::T1,
    u::T2,
    prim::T3,
) where {
    T1<:AbstractArray{<:AbstractFloat,2},
    T2<:AbstractArray{<:AbstractFloat,2},
    T3<:AbstractArray{<:Real,2},
}

    for j in axes(M, 2)
        _M = @view M[:, j]
        maxwellian!(_M, u[:, j], prim[:, j])
    end

    return nothing

end

function mixture_maxwellian!(
    M::T1,
    u::T2,
    v::T2,
    prim::T3,
) where {
    T1<:AbstractArray{<:AbstractFloat,3},
    T2<:AbstractArray{<:AbstractFloat,3},
    T3<:AbstractArray{<:Real,2},
}

    for k in axes(M, 3)
        _M = @view M[:, :, k]
        maxwellian!(_M, u[:, :, k], v[:, :, k], prim[:, k])
    end

    return nothing

end

function mixture_maxwellian!(
    M::T1,
    u::T2,
    v::T2,
    w::T2,
    prim::T3,
) where {
    T1<:AbstractArray{<:AbstractFloat,4},
    T2<:AbstractArray{<:AbstractFloat,4},
    T3<:AbstractArray{<:Real,2},
}

    for l in axes(M, 4)
        _M = @view M[:, :, :, l]
        maxwellian!(_M, u[:, :, :, l], v[:, :, :, l], w[:, :, :, l], prim[:, l])
    end

    return nothing

end
