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

maxwellian(u, prim::AbstractVector) = maxwellian(u, prim[1], prim[2], prim[end]) # in case of input with length 4/5

#--- 2V ---#
maxwellian(u, v, ρ, U, V, λ) = @. ρ * (λ / π) * exp(-λ * ((u - U)^2 + (v - V)^2))

maxwellian(u, v, prim::AbstractVector) =
    maxwellian(u, v, prim[1], prim[2], prim[3], prim[end]) # in case of input with length 5

#--- 3V ---#
maxwellian(u, v, w, ρ, U, V, W, λ) =
    @. ρ * sqrt((λ / π)^3) * exp(-λ * ((u - U)^2 + (v - V)^2 + (w - W)^2))

maxwellian(u, v, w, prim::AbstractVector) =
    maxwellian(u, v, w, prim[1], prim[2], prim[3], prim[4], prim[5])


"""
    * Maxwellian
    maxwellian!(M, u, ρ, U, λ)
    maxwellian!(M, u, ρ, U, V, λ)
    maxwellian!(M, u, ρ, U, V, W, λ)
    maxwellian!(M, u, prim)
    maxwellian!(M, u, v, prim)
    maxwellian!(M, u, v, w, prim)

    * Rykov
    maxwellian!(Ht, Bt, Rt, Hr, Br, Rr, u, prim, K, Kr)
    maxwellian!(Ht, Bt, Rt, Hr, Br, Rr, u, v, prim, K, Kr)

In-place Maxwellian

* @args: particle velocity quadrature points
* @args: density, velocity and inverse of temperature
* @return: Maxwellian distribution function

"""
function maxwellian!(
    M::AbstractVector,
    u::AbstractVector,
    ρ,
    U,
    λ,
)

    @. M = ρ * sqrt(λ / π) * exp(-λ * (u - U)^2) # 1V
    return nothing

end

maxwellian!(
    M::AbstractVector,
    u::AbstractVector,
    prim::AbstractVector,
) = maxwellian!(M, u, prim[1], prim[2], prim[end])

# Rykov
function maxwellian!(
    Ht::T,
    Bt::T,
    Rt::T,
    Hr::T,
    Br::T,
    Rr::T,
    u::AbstractVector,
    prim::AbstractVector,
    K,
    Kr,
) where {T<:AbstractVector}

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
    M::AbstractArray,
    u::T,
    v::T,
    ρ,
    U,
    V,
    λ,
) where {T<:AbstractArray}

    @. M = ρ * (λ / π) * exp(-λ * ((u - U)^2 + (v - V)^2))
    return nothing

end

maxwellian!(
    M::AbstractArray,
    u::T,
    v::T,
    prim::AbstractVector,
) where {T<:AbstractArray} = maxwellian!(M, u, v, prim[1], prim[2], prim[3], prim[end])

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
    prim::AbstractVector,
    K,
    Kr,
) where {
    T1<:AbstractArray,
    T2<:AbstractArray,
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
    M::AbstractArray,
    u::T,
    v::T,
    w::T,
    ρ,
    U,
    V,
    W,
    λ,
) where {T<:AbstractArray}
    @. M = ρ * sqrt((λ / π)^3) * exp(-λ * ((u - U)^2 + (v - V)^2 + (w - W)^2))
    return nothing
end

maxwellian!(
    M::AbstractArray,
    u::T,
    v::T,
    w::T,
    prim::AbstractVector,
) where {
    T<:AbstractArray,
} = maxwellian!(M, u, v, w, prim[1], prim[2], prim[3], prim[4], prim[5])


"""
    mixture_maxwellian(u, prim)
    mixture_maxwellian(u, v, prim)
    mixture_maxwellian(u, v, w, prim)

Multi-component Maxwellian in discrete form
"""
function mixture_maxwellian(
    u::AbstractMatrix,
    prim::AbstractMatrix,
)

    mixM = similar(u)
    mixture_maxwellian!(mixM, u, prim)

    return mixM

end

#--- 2V ---#
function mixture_maxwellian(
    u::X,
    v::X,
    prim::AbstractMatrix,
) where {X<:AbstractArray}

    mixM = similar(u)
    mixture_maxwellian!(mixM, u, v, prim)

    return mixM

end

#--- 3V ---#
function mixture_maxwellian(
    u::X,
    v::X,
    w::X,
    prim::AbstractMatrix,
) where {X<:AbstractArray}

    mixM = similar(u)
    mixture_maxwellian!(mixM, u, v, w, prim)

    return mixM

end


"""
    mixture_maxwellian!(M, u, prim)
    mixture_maxwellian!(M, u, v, prim)
    mixture_maxwellian!(M, u, v, w, prim)

In-place multi-component Maxwellian

"""
function mixture_maxwellian!(
    M::AbstractMatrix,
    u::AbstractArray{T,2},
    prim::AbstractMatrix,
) where {T}

    for j in axes(M, 2)
        _M = @view M[:, j]
        maxwellian!(_M, u[:, j], prim[:, j])
    end

    return nothing

end

#--- 2V ---#
function mixture_maxwellian!(
    M::AbstractArray{T1,3},
    u::T2,
    v::T2,
    prim::AbstractMatrix,
) where {T1,T2<:AbstractArray{T3,3}} where T3

    for k in axes(M, 3)
        _M = @view M[:, :, k]
        maxwellian!(_M, u[:, :, k], v[:, :, k], prim[:, k])
    end

    return nothing

end

function mixture_maxwellian!(
    M::AbstractMatrix,
    u::T,
    v::T,
    prim::AbstractMatrix,
) where {T<:AbstractMatrix}

    for k in axes(M, 2)
        _M = @view M[:, k]
        maxwellian!(_M, u[:, k], v[:, k], prim[:, k])
    end

    return nothing

end

#--- 3V ---#
function mixture_maxwellian!(
    M::AbstractArray{T1,4},
    u::T2,
    v::T2,
    w::T2,
    prim::AbstractMatrix,
) where {T1,T2<:AbstractArray{T3,4}} where T3

    for l in axes(M, 4)
        _M = @view M[:, :, :, l]
        maxwellian!(_M, u[:, :, :, l], v[:, :, :, l], w[:, :, :, l], prim[:, l])
    end

    return nothing

end

function mixture_maxwellian!(
    M::AbstractMatrix,
    u::T,
    v::T,
    w::T,
    prim::AbstractMatrix,
) where {T<:AbstractMatrix}

    for l in axes(M, 2)
        _M = @view M[:, l]
        maxwellian!(_M, u[:, l], v[:, l], w[:, l], prim[:, l])
    end

    return nothing

end
