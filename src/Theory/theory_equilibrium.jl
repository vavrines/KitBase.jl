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

Compute Maxwellian directly from distribution function
"""
function f_maxwellian(f::AV, u::AV, weights::AV, γ = 3)
    w = moments_conserve(f, u, weights)
    prim = conserve_prim(w, γ)

    return maxwellian(u, prim)
end

"""
$(SIGNATURES)
"""
f_maxwellian(
    f::AV,
    vs = VSpace1D(-6, 6, size(f, 1); precision = Float32)::VSpace1D,
    γ = 3,
) = f_maxwellian(f, vs.u, vs.weights, γ)

"""
$(SIGNATURES)
"""
function f_maxwellian(
    f::AM,
    vs = VSpace1D(-6, 6, size(f, 1); precision = Float32)::VSpace1D,
    γ = 3,
)
    M = [f_maxwellian(f[:, i], vs, γ) for i in axes(f, 2)]

    return hcat(M...)
end

"""
$(SIGNATURES)
"""
function f_maxwellian(h::AV, b::AV, u::AV, weights::AV, γ = 5 / 3)
    w = moments_conserve(h, b, u, weights)
    prim = conserve_prim(w, γ)
    H = maxwellian(u, prim)
    B = energy_maxwellian(H, prim[end], internal_dof(γ, 1))

    return H, B
end

f_maxwellian(
    h::AV,
    b::AV,
    vs = VSpace1D(-6, 6, size(h, 1); precision = Float32)::VSpace1D,
    γ = 5 / 3,
) = f_maxwellian(h, b, vs.u, vs.weights, γ)

function f_maxwellian(
    h::AM,
    b::AM,
    vs = VSpace1D(-6, 6, size(h, 1); precision = Float32)::VSpace1D,
    γ = 5 / 3,
)
    H = zero(h)
    B = zero(b)
    for j in axes(h, 2)
        H[:, j], B[:, j] = f_maxwellian(h[:, j], b[:, j], vs, γ)
    end

    return H, B
end

# ------------------------------------------------------------
# Shakhov model
# ------------------------------------------------------------

"""
$(SIGNATURES)

Shakhov non-equilibrium part

1F1V
"""
function shakhov(u::AV, M::T, q, prim::AV, Pr)::T where {T<:AV}
    M_plus = @. 0.8 * (1.0 - Pr) * prim[end]^2 / prim[1] *
       (u - prim[2]) *
       q *
       (2.0 * prim[end] * (u - prim[2])^2 - 5.0) *
       M

    return M_plus
end

"""
$(SIGNATURES)

2F1V
"""
function shakhov(u::AV, H::T, B::T, q, prim::AV, Pr, K) where {T<:AV}

    H_plus = @. 0.8 * (1.0 - Pr) * prim[end]^2 / prim[1] *
       (u - prim[2]) *
       q *
       (2.0 * prim[end] * (u - prim[2])^2 + K - 5.0) *
       H
    B_plus = @. 0.8 * (1.0 - Pr) * prim[end]^2 / prim[1] *
       (u - prim[2]) *
       q *
       (2.0 * prim[end] * (u - prim[2])^2 + K - 3.0) *
       B

    return H_plus, B_plus

end

"""
$(SIGNATURES)

1F2V
"""
function shakhov(u::T, v::T, M::T1, q::AV, prim::AV, Pr)::T1 where {T<:AM,T1<:AM}

    M_plus = @. 0.8 * (1.0 - Pr) * prim[end]^2 / prim[1] *
       ((u - prim[2]) * q[1] + (v - prim[3]) * q[2]) *
       (2.0 * prim[end] * ((u - prim[2])^2 + (v - prim[3])^2) - 5.0) *
       M

    return M_plus

end

"""
$(SIGNATURES)

2F2V
"""
function shakhov(u::T, v::T, H::X, B::X, q::AV, prim::AV, Pr, K) where {T<:AM,X<:AM}

    H_plus = @. 0.8 * (1.0 - Pr) * prim[end]^2 / prim[1] *
       ((u - prim[2]) * q[1] + (v - prim[3]) * q[2]) *
       (2.0 * prim[end] * ((u - prim[2])^2 + (v - prim[3])^2) + K - 5.0) *
       H
    B_plus = @. 0.8 * (1.0 - Pr) * prim[end]^2 / prim[1] *
       ((u - prim[2]) * q[1] + (v - prim[3]) * q[2]) *
       (2.0 * prim[end] * ((u - prim[2])^2 + (v - prim[3])^2) + K - 3.0) *
       B

    return H_plus, B_plus

end

"""
$(SIGNATURES)

1F3V
"""
function shakhov(
    u::T,
    v::T,
    w::T,
    M::T1,
    q::AV,
    prim::AV,
    Pr,
)::T1 where {T<:AA{T2,3},T1<:AA{T3,3}} where {T2,T3}

    M_plus = @. 0.8 * (1.0 - Pr) * prim[end]^2 / prim[1] *
       ((u - prim[2]) * q[1] + (v - prim[3]) * q[2] + (w - prim[4]) * q[3]) *
       (2.0 * prim[end] * ((u - prim[2])^2 + (v - prim[3])^2 + (w - prim[4])^2) - 5.0) *
       M

    return M_plus

end


"""
$(SIGNATURES)

In-place Shakhov non-equilibrium part

1F1V
"""
function shakhov!(S::T1, u::AV, M::T1, q, prim, Pr) where {T1<:AV}
    @. S = @. 0.8 * (1.0 - Pr) * prim[end]^2 / prim[1] *
       (u - prim[2]) *
       q *
       (2.0 * prim[end] * (u - prim[2])^2 - 5.0) *
       M

    return nothing
end

"""
$(SIGNATURES)

2F1V
"""
function shakhov!(SH::T1, SB::T1, u::AV, H::T1, B::T1, q, prim, Pr, K) where {T1<:AV}

    @. SH =
        0.8 * (1.0 - Pr) * prim[end]^2 / prim[1] *
        (u - prim[2]) *
        q *
        (2.0 * prim[end] * (u - prim[2])^2 + K - 5.0) *
        H
    @. SB =
        0.8 * (1.0 - Pr) * prim[end]^2 / prim[1] *
        (u - prim[2]) *
        q *
        (2.0 * prim[end] * (u - prim[2])^2 + K - 3.0) *
        B

    return nothing

end

"""
$(SIGNATURES)

1F2V
"""
function shakhov!(S::T1, u::T2, v::T2, M::T1, q::AV, prim, Pr) where {T1<:AM,T2<:AM}

    @. S =
        0.8 * (1.0 - Pr) * prim[end]^2 / prim[1] *
        ((u - prim[2]) * q[1] + (v - prim[3]) * q[2]) *
        (2.0 * prim[end] * ((u - prim[2])^2 + (v - prim[3])^2) - 5.0) *
        M

    return nothing

end

"""
$(SIGNATURES)

2F2V
"""
function shakhov!(
    SH::T1,
    SB::T1,
    u::T2,
    v::T2,
    H::T1,
    B::T1,
    q::AV,
    prim,
    Pr,
    K,
) where {T1<:AM,T2<:AM}

    @. SH =
        0.8 * (1.0 - Pr) * prim[end]^2 / prim[1] *
        ((u - prim[2]) * q[1] + (v - prim[3]) * q[2]) *
        (2.0 * prim[end] * ((u - prim[2])^2 + (v - prim[3])^2) + K - 5.0) *
        H
    @. SB =
        0.8 * (1.0 - Pr) * prim[end]^2 / prim[1] *
        ((u - prim[2]) * q[1] + (v - prim[3]) * q[2]) *
        (2.0 * prim[end] * ((u - prim[2])^2 + (v - prim[3])^2) + K - 3.0) *
        B

    return nothing

end

"""
$(SIGNATURES)

1F3V
"""
function shakhov!(
    S::T1,
    u::T2,
    v::T2,
    w::T2,
    M::T1,
    q::AV,
    prim,
    Pr,
) where {T1<:AA{T3,3},T2<:AA{T4,3}} where {T3,T4}

    @. S =
        0.8 * (1.0 - Pr) * prim[end]^2 / prim[1] *
        ((u - prim[2]) * q[1] + (v - prim[3]) * q[2] + (w - prim[4]) * q[3]) *
        (2.0 * prim[end] * ((u - prim[2])^2 + (v - prim[3])^2 + (w - prim[4])^2) - 5.0) *
        M

    return nothing

end

# ------------------------------------------------------------
# ES-BGK model
# ------------------------------------------------------------

"""
$(SIGNATURES)

Compute Gaussian distribution in ES-BGK model

1F1V
"""
function esbgk(f, u, weights, prim, Pr)
    P = stress(f, prim, u, weights) / prim[1]
    T = 0.5 / prim[end]
    ν = 1 - 1 / Pr
    M = (1 - ν) * T + ν * P

    prod = zero(u)
    for i in eachindex(u)
        c = u[i] - prim[2]
        prod[i] = c^2 / M
    end

    return prim[1] / (sqrt(2π * M)) .* exp.(-prod ./ 2)
end

"""
$(SIGNATURES)

1F2V
"""
function esbgk(f, u, v, weights, prim, Pr)
    P = stress(f, prim, u, v, weights) / prim[1]
    t = 0.5 / prim[end]
    T = diagm([t, t])
    ν = 1 - 1 / Pr
    M = (1 - ν) * T + ν * P
    MI = inv(M)

    c = zeros(2)
    prod = zero(u)
    for i in eachindex(u)
        c[1] = u[i] - prim[2]
        c[2] = v[i] - prim[3]
        prod[i] = c' * MI * c
    end

    return prim[1] / det(sqrt(2π * M)) .* exp.(-prod ./ 2)
end

"""
$(SIGNATURES)

1F3V
"""
function esbgk(f, u, v, w, weights, prim, Pr)
    P = stress(f, prim, u, v, w, weights) / prim[1]
    t = 0.5 / prim[end]
    T = diagm([t, t, t])
    ν = 1 - 1 / Pr
    M = (1 - ν) * T + ν * P
    MI = inv(M)

    c = zeros(3)
    prod = zero(u)
    for i in eachindex(u)
        c[1] = u[i] - prim[2]
        c[2] = v[i] - prim[3]
        c[3] = w[i] - prim[4]
        prod[i] = c' * MI * c
    end

    return prim[1] / det(sqrt(2π * M)) .* exp.(-prod ./ 2)
end

# ------------------------------------------------------------
# Polyatomic gas
# ------------------------------------------------------------

"""
$(SIGNATURES)

Rykov model
"""
function polyatomic_maxwellian!(
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
) where {T}

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

Rykov model
"""
function polyatomic_maxwellian!(
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
) where {T1,T2}

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

Rykov non-equilibrium part 

- t: translation
- r: rotation

## Arguments
* `q`: [qt, qr]@1D, [qtx, qty, qrx, qry]@2D
* `prim`: [ρ, U, λ, λₜ, λᵣ]@1D, [ρ, U, V, λ, λₜ, λᵣ]@2D
* `σ`, ω0, ω1: rotation parameters
"""
function rykov!(
    Ht_plus::T1,
    Bt_plus::T1,
    Rt_plus::T1,
    Hr_plus::T1,
    Br_plus::T1,
    Rr_plus::T1,
    u::T2,
    Ht::T1,
    Bt::T1,
    Rt::T1,
    Hr::T1,
    Br::T1,
    Rr::T1,
    q::T3,
    prim::AV,
    Pr,
    K,
    σ,
    ω0,
    ω1,
) where {T1,T2,T3}

    @. Ht_plus =
        (
            0.8 * (1.0 - Pr) * prim[4]^2 / prim[1] *
            (u - prim[2]) *
            q[1] *
            (2.0 * prim[4] * (u - prim[2])^2 + K - 5.0)
        ) * Ht
    @. Bt_plus =
        (
            0.8 * (1.0 - Pr) * prim[4]^2 / prim[1] *
            (u - prim[2]) *
            q[1] *
            (2.0 * prim[4] * (u - prim[2])^2 + K - 3.0)
        ) * Bt
    @. Rt_plus =
        (
            0.8 * (1.0 - Pr) * prim[4]^2 / prim[1] *
            (u - prim[2]) *
            q[1] *
            (2.0 * prim[4] * (u - prim[2])^2 + K - 5.0) +
            4.0 * (1 - σ) * (u - prim[2]) * q[2] / prim[1] * prim[4] * prim[5]
        ) * Rt

    @. Hr_plus =
        (
            0.8 * ω0 * (1.0 - Pr) * prim[3]^2 / prim[1] *
            (u - prim[2]) *
            q[1] *
            (2.0 * prim[3] * (u - prim[2])^2 + K - 5.0)
        ) * Hr
    @. Br_plus =
        (
            0.8 * ω0 * (1.0 - Pr) * prim[3]^2 / prim[1] *
            (u - prim[2]) *
            q[1] *
            (2.0 * prim[3] * (u - prim[2])^2 + K - 3.0)
        ) * Br
    @. Rr_plus =
        (
            0.8 * ω0 * (1.0 - Pr) * prim[3]^2 / prim[1] *
            (u - prim[2]) *
            q[1] *
            (2.0 * prim[3] * (u - prim[2])^2 + K - 5.0) +
            4.0 * ω1 * (1 - σ) * (u - prim[2]) * q[2] / prim[1] * prim[3] * prim[3]
        ) * Rr

    return nothing

end

"""
$(SIGNATURES)
"""
function rykov!(
    Ht_plus::T1,
    Bt_plus::T1,
    Rt_plus::T1,
    Hr_plus::T1,
    Br_plus::T1,
    Rr_plus::T1,
    u::T2,
    v::T2,
    Ht::T1,
    Bt::T1,
    Rt::T1,
    Hr::T1,
    Br::T1,
    Rr::T1,
    q::T3,
    prim::AV,
    Pr,
    K,
    σ,
    ω0,
    ω1,
) where {T1,T2,T3}

    @. Ht_plus =
        (
            0.8 * (1.0 - Pr) * prim[5]^2 / prim[1] *
            ((u - prim[2]) * q[1] + (v - prim[3]) * q[2]) *
            (2.0 * prim[5] * ((u - prim[2])^2 + (v - prim[3])^2) + K - 5.0)
        ) * Ht
    @. Bt_plus =
        (
            0.8 * (1.0 - Pr) * prim[5]^2 / prim[1] *
            ((u - prim[2]) * q[1] + (v - prim[3]) * q[2]) *
            (2.0 * prim[5] * ((u - prim[2])^2 + (v - prim[3])^2) + K - 3.0)
        ) * Bt
    @. Rt_plus =
        (
            0.8 * (1.0 - Pr) * prim[5]^2 / prim[1] *
            ((u - prim[2]) * q[1] + (v - prim[3]) * q[2]) *
            (2.0 * prim[5] * ((u - prim[2])^2 + (v - prim[3])^2) + K - 5) +
            4.0 * (1 - σ) * ((u - prim[2]) * q[3] + (v - prim[3]) * q[4]) / prim[1] *
            prim[5] *
            prim[6]
        ) * Rt

    @. Hr_plus =
        (
            0.8 * ω0 * (1.0 - Pr) * prim[4]^2 / prim[1] *
            ((u - prim[2]) * q[1] + (v - prim[3]) * q[2]) *
            (2.0 * prim[4] * ((u - prim[2])^2 + (v - prim[3])^2) + K - 5.0)
        ) * Hr
    @. Br_plus =
        (
            0.8 * ω0 * (1.0 - Pr) * prim[4]^2 / prim[1] *
            ((u - prim[2]) * q[1] + (v - prim[3]) * q[2]) *
            (2.0 * prim[4] * ((u - prim[2])^2 + (v - prim[3])^2) + K - 3.0)
        ) * Br
    @. Rr_plus =
        (
            0.8 * ω0 * (1.0 - Pr) * prim[4]^2 / prim[1] *
            ((u - prim[2]) * q[1] + (v - prim[3]) * q[2]) *
            (2.0 * prim[4] * ((u - prim[2])^2 + (v - prim[3])^2) + K - 5.0) +
            4.0 * ω1 * (1 - σ) * ((u - prim[2]) * q[3] + (v - prim[3]) * q[4]) / prim[1] *
            prim[4] *
            prim[4]
        ) * Rr

    return nothing

end


"""
$(SIGNATURES)

BIP model (1V)

_F. Bernard, A. Iollo, and G. Puppo. BGK polyatomic model for rarefied flows. Journal of Scientific Computing 78 (2019): 1893-1916._
"""
function polyatomic_maxwellian!(M1, M2, M3, u, prim, K, Kr)
    @. M1 = prim[1] * (prim[4] / π)^0.5 * exp(-prim[4] * (u - prim[2])^2)
    @. M2 = M1 * K / 2.0 / prim[4]
    @. M3 = M1 * Kr / 2.0 / prim[5]

    return nothing
end

"""
$(SIGNATURES)

BIP model (2V)
"""
function polyatomic_maxwellian!(M1, M2, M3, u, v, prim, K, Kr)
    @. M1 = prim[1] * (prim[5] / π) * exp(-prim[5] * ((u - prim[2])^2 + (v - prim[3])^2))
    @. M2 = M1 * K / (2.0 * prim[5])
    @. M3 = M1 * Kr / (2.0 * prim[6])

    return nothing
end

"""
$(SIGNATURES)

BIP model (3V)
"""
function polyatomic_maxwellian!(M1, M2, M3, u, v, w, prim, K, Kr)
    @. M1 = prim[1] * (prim[6] / π) * exp(-prim[6] * ((u - prim[2])^2 + (v - prim[3])^2 + (w - prim[4])^2))
    @. M2 = M1 * K / (2.0 * prim[6])
    @. M3 = M1 * Kr / (2.0 * prim[7])

    return nothing
end

# ------------------------------------------------------------
# Multi-species gas
# ------------------------------------------------------------

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

Compute Maxwellian for polyatomic gas mixture
"""
function mixture_polyatomic_maxwellian!(M1, M2, M3, u, prim, K, Kr)
    for i in axes(M1)[end]
        _M1 = extract_last(M1, i)
        _M2 = extract_last(M2, i)
        _M3 = extract_last(M3, i)
        _u = extract_last(u, i)
        polyatomic_maxwellian!(_M1, _M2, _M3, _u, prim[:, i], K, Kr)
    end

    return nothing
end

"""
$(SIGNATURES)
"""
function mixture_polyatomic_maxwellian!(M1, M2, M3, u, v, prim, K, Kr)
    for i in axes(M1)[end]
        _M1 = extract_last(M1, i)
        _M2 = extract_last(M2, i)
        _M3 = extract_last(M3, i)
        _u = extract_last(u, i)
        _v = extract_last(v, i)
        polyatomic_maxwellian!(_M1, _M2, _M3, _u, _v, prim[:, i], K, Kr)
    end

    return nothing
end

"""
$(SIGNATURES)
"""
function mixture_polyatomic_maxwellian!(M1, M2, M3, u, v, w, prim, K, Kr)
    for i in axes(M1)[end]
        _M1 = extract_last(M1, i)
        _M2 = extract_last(M2, i)
        _M3 = extract_last(M3, i)
        _u = extract_last(u, i)
        _v = extract_last(v, i)
        _w = extract_last(w, i)
        polyatomic_maxwellian!(_M1, _M2, _M3, _u, _v, _w, prim[:, i], K, Kr)
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
