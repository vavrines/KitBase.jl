"""
    rykov!(
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
        prim::T4,
        Pr,
        K,
        σ,
        ω0,
        ω1,
    ) where {
        T1<:AA{<:Real,1},
        T2<:AA{<:Real,1},
        T3<:AA{<:Real,1},
        T4<:AA{<:Real,1},
    }

    rykov!(
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
        prim::T4,
        Pr,
        K,
        σ,
        ω0,
        ω1,
    ) where {
        T1<:AA{<:Real,2},
        T2<:AA{<:Real,2},
        T3<:AA{<:Real,1},
        T4<:AA{<:Real,1},
    }

Rykov non-equilibrium part 

- t: translation
- r: rotation
* @arg q: [qt, qr]@1D, [qtx, qty, qrx, qry]@2D
* @arg prim: [ρ, U, λ, λₜ, λᵣ]@1D, [ρ, U, V, λ, λₜ, λᵣ]@2D
* @arg σ, ω0, ω1: rotation parameters

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
    prim::T4,
    Pr,
    K,
    σ,
    ω0,
    ω1,
) where {
    T1<:AbstractVector{<:Real},
    T2<:AbstractVector{<:Real},
    T3<:AbstractVector{<:Real},
    T4<:AbstractVector{<:Real},
}

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
    prim::T4,
    Pr,
    K,
    σ,
    ω0,
    ω1,
) where {
    T1<:AA{<:Real,2},
    T2<:AA{<:Real,2},
    T3<:AA{<:Real,1},
    T4<:AA{<:Real,1},
}

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
    rykov_zr(Tₜ, T₀, Z₀)

Calculate dimensionless rotationa number in Rykov model

"""
rykov_zr(Tₜ, T₀, Z₀) = Z₀ / (1.0 + π^1.5 / 2.0 * sqrt(T₀ / Tₜ) + (π + 0.25 * π^2) * T₀ / Tₜ)
