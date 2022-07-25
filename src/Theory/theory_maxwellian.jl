"""
$(SIGNATURES)

Maxwellian in discrete form

1V
"""
maxwellian(u, ρ, U, λ) = @. ρ * sqrt(λ / π) * exp(-λ * (u - U)^2)

"""
$(SIGNATURES)
"""
maxwellian(u, prim::AV) = maxwellian(u, prim[1], prim[2], prim[end]) # in case of input with length 4/5

"""
$(SIGNATURES)

2V
"""
maxwellian(u, v, ρ, U, V, λ) = @. ρ * (λ / π) * exp(-λ * ((u - U)^2 + (v - V)^2))

"""
$(SIGNATURES)
"""
maxwellian(u, v, prim::AV) = maxwellian(u, v, prim[1], prim[2], prim[3], prim[end]) # in case of input with length 5

"""
$(SIGNATURES)

3V
"""
maxwellian(u, v, w, ρ, U, V, W, λ) =
    @. ρ * sqrt((λ / π)^3) * exp(-λ * ((u - U)^2 + (v - V)^2 + (w - W)^2))

"""
$(SIGNATURES)
"""
maxwellian(u, v, w, prim::AV) =
    maxwellian(u, v, w, prim[1], prim[2], prim[3], prim[4], prim[5])


"""
$(SIGNATURES)

In-place Maxwellian

1V
"""
function maxwellian!(M::AV, u::AV, ρ, U, λ)
    @. M = ρ * sqrt(λ / π) * exp(-λ * (u - U)^2)
    return nothing
end

"""
$(SIGNATURES)
"""
maxwellian!(M::AV, u::AV, prim::AV) = maxwellian!(M, u, prim[1], prim[2], prim[end])

"""
$(SIGNATURES)

Rykov model
"""
function maxwellian!(
    Ht::T,
    Bt::T,
    Rt::T,
    Hr::T,
    Br::T,
    Rr::T,
    u::AV,
    prim::AV,
    K,
    Kr,
) where {T<:AV}

    @. Ht = prim[1] * sqrt(prim[4] / π) * exp(-prim[4] * (u - prim[2])^2)
    @. Bt = Ht * K / (2.0 * prim[4])
    @. Rt = Ht * Kr / (2.0 * prim[5])

    @. Hr = prim[1] * sqrt(prim[3] / π) * exp(-prim[3] * (u - prim[2])^2)
    @. Br = Hr * K / (2.0 * prim[3])
    @. Rr = Hr * Kr / (2.0 * prim[3])

    return nothing

end

"""
$(SIGNATURES)

2V
"""
function maxwellian!(M::AA, u::T, v::T, ρ, U, V, λ) where {T<:AA}
    @. M = ρ * (λ / π) * exp(-λ * ((u - U)^2 + (v - V)^2))
    return nothing
end

"""
$(SIGNATURES)
"""
maxwellian!(M::AA, u::T, v::T, prim::AV) where {T<:AA} =
    maxwellian!(M, u, v, prim[1], prim[2], prim[3], prim[end])

"""
$(SIGNATURES)

Rykov model
"""
function maxwellian!(
    Ht::T1,
    Bt::T1,
    Rt::T1,
    Hr::T1,
    Br::T1,
    Rr::T1,
    u::T2,
    v::T2,
    prim::AV,
    K,
    Kr,
) where {T1<:AA,T2<:AA}

    @. Ht = prim[1] * (prim[5] / π) * exp(-prim[5] * ((u - prim[2])^2 + (v - prim[3])^2))
    @. Bt = Ht * K / (2.0 * prim[5])
    @. Rt = Ht * Kr / (2.0 * prim[6])

    @. Hr = prim[1] * (prim[4] / π) * exp(-prim[4] * ((u - prim[2])^2 + (v - prim[3])^2))
    @. Br = Hr * K / (2.0 * prim[4])
    @. Rr = Hr * Kr / (2.0 * prim[4])

    return nothing

end

"""
$(SIGNATURES)

3V
"""
function maxwellian!(M::AA, u::T, v::T, w::T, ρ, U, V, W, λ) where {T<:AA}
    @. M = ρ * sqrt((λ / π)^3) * exp(-λ * ((u - U)^2 + (v - V)^2 + (w - W)^2))
    return nothing
end

"""
$(SIGNATURES)
"""
maxwellian!(M::AA, u::T, v::T, w::T, prim::AV) where {T<:AA} =
    maxwellian!(M, u, v, w, prim[1], prim[2], prim[3], prim[4], prim[5])


"""
$(SIGNATURES)

Reduced Maxwellian distribution related to energy
"""
energy_maxwellian(h, λ, K) = @. h * K / (2.0 * λ)

"""
$(SIGNATURES)
"""
energy_maxwellian(h, prim::AV, K) = energy_maxwellian(h, prim[end], K)


"""
$(SIGNATURES)

Multi-component Maxwellian in discrete form

1V
"""
function mixture_maxwellian(u::AM, prim::AM)
    mixM = similar(u)
    mixture_maxwellian!(mixM, u, prim)

    return mixM
end

"""
$(SIGNATURES)

2V
"""
function mixture_maxwellian(u::X, v::X, prim::AM) where {X<:AA}
    mixM = similar(u)
    mixture_maxwellian!(mixM, u, v, prim)

    return mixM
end

"""
$(SIGNATURES)

3V
"""
function mixture_maxwellian(u::X, v::X, w::X, prim::AM) where {X<:AA}
    mixM = similar(u)
    mixture_maxwellian!(mixM, u, v, w, prim)

    return mixM
end


"""
$(SIGNATURES)

In-place multi-component Maxwellian

1V
"""
function mixture_maxwellian!(M::AM, u::AA{T,2}, prim::AM) where {T}
    @inbounds for j in axes(M, 2)
        _M = @view M[:, j]
        maxwellian!(_M, u[:, j], prim[:, j])
    end

    return nothing
end

"""
$(SIGNATURES)

2V
"""
function mixture_maxwellian!(
    M::AA{T1,3},
    u::T2,
    v::T2,
    prim::AM,
) where {T1,T2<:AA{T3,3}} where {T3}

    @inbounds for k in axes(M, 3)
        _M = @view M[:, :, k]
        maxwellian!(_M, u[:, :, k], v[:, :, k], prim[:, k])
    end

    return nothing

end

"""
$(SIGNATURES)
"""
function mixture_maxwellian!(M::AM, u::T, v::T, prim::AM) where {T<:AM}

    @inbounds for k in axes(M, 2)
        _M = @view M[:, k]
        maxwellian!(_M, u[:, k], v[:, k], prim[:, k])
    end

    return nothing

end

"""
$(SIGNATURES)

3V
"""
function mixture_maxwellian!(
    M::AA{T1,4},
    u::T2,
    v::T2,
    w::T2,
    prim::AM,
) where {T1,T2<:AA{T3,4}} where {T3}

    @inbounds for l in axes(M, 4)
        _M = @view M[:, :, :, l]
        maxwellian!(_M, u[:, :, :, l], v[:, :, :, l], w[:, :, :, l], prim[:, l])
    end

    return nothing

end

function mixture_maxwellian!(M::AM, u::T, v::T, w::T, prim::AM) where {T<:AM}

    @inbounds for l in axes(M, 2)
        _M = @view M[:, l]
        maxwellian!(_M, u[:, l], v[:, l], w[:, l], prim[:, l])
    end

    return nothing

end

"""
$(SIGNATURES)

Reduced Maxwellian distribution related to energy for mixture
"""
function mixture_energy_maxwellian(h::AM, prim::AM, K)
    b = zero(h)
    for j in axes(prim, 2)
        b[:, j] .= energy_maxwellian(h[:, j], prim[:, j], K)
    end

    return b
end

function mixture_energy_maxwellian(h::AA{T,3}, prim::AM, K) where {T}
    b = zero(h)
    for j in axes(prim, 2)
        b[:, :, j] .= energy_maxwellian(h[:, :, j], prim[:, j], K)
    end

    return b
end
