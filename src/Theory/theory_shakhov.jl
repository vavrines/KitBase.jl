"""
    shakhov(u, M, q, prim, Pr)
    shakhov(u, H, B, q, prim, Pr, K)
    shakhov(u, v, M, q, prim, Pr)
    shakhov(u, v, H, B, q, prim, Pr, K)
    shakhov(u, v, w, M, q, prim, Pr)

Shakhov non-equilibrium part

* @arg: particle velocity quadrature points
* @arg: discrete Maxwellian
* @arg: primitive variables, Prandtl number, heat flux, inner degree of freedom

"""
function shakhov(
    u::AbstractVector{X},
    M::AbstractVector{Y},
    q,
    prim::AbstractVector{Z},
    Pr,
) where {X<:AbstractFloat,Y<:AbstractFloat,Z<:Real} # 1F1V

    M_plus = @. 0.8 * (1.0 - Pr) * prim[end]^2 / prim[1] *
       (u - prim[2]) *
       q *
       (2.0 * prim[end] * (u - prim[2])^2 - 5.0) *
       M

    return M_plus

end

#--- 2F1V ---#
function shakhov(
    u::AbstractVector{X},
    H::Y,
    B::Y,
    q,
    prim::AbstractVector{Z},
    Pr,
    K,
) where {X<:AbstractFloat,Y<:AbstractVector{<:AbstractFloat},Z<:Real}

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

#--- 1F2V ---#
function shakhov(
    u::T,
    v::T,
    M::AbstractMatrix{X},
    q::AbstractVector{Y},
    prim::AbstractVector{Z},
    Pr,
) where {T<:AbstractArray{<:AbstractFloat,2},X<:AbstractFloat,Y<:AbstractFloat,Z<:Real}

    M_plus = @. 0.8 * (1.0 - Pr) * prim[end]^2 / prim[1] *
       ((u - prim[2]) * q[1] + (v - prim[3]) * q[2]) *
       (2.0 * prim[end] * ((u - prim[2])^2 + (v - prim[3])^2) - 5.0) *
       M

    return M_plus

end

#--- 2F2V ---#
function shakhov(
    u::T,
    v::T,
    H::X,
    B::X,
    q::AbstractVector{Y},
    prim::AbstractVector{Z},
    Pr,
    K,
) where {
    T<:AbstractArray{<:AbstractFloat,2},
    X<:AbstractArray{<:AbstractFloat,2},
    Y<:Real,
    Z<:Real,
}

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

#--- 1F3V ---#
function shakhov(
    u::T,
    v::T,
    w::T,
    M::AbstractArray{X,3},
    q::AbstractVector{Y},
    prim::AbstractVector{Z},
    Pr,
) where {T<:AbstractArray{<:AbstractFloat,3},X<:AbstractFloat,Y<:Real,Z<:Real}

    M_plus = @. 0.8 * (1.0 - Pr) * prim[end]^2 / prim[1] *
       ((u - prim[2]) * q[1] + (v - prim[3]) * q[2] + (w - prim[4]) * q[3]) *
       (2.0 * prim[end] * ((u - prim[2])^2 + (v - prim[3])^2 + (w - prim[4])^2) - 5.0) *
       M

    return M_plus

end


"""
    shakhov!(S, u, M, q, prim, Pr)
    shakhov!(SH, SB, u, H, B, q, prim, Pr, K)
    shakhov!(S, u, v, M, q, prim, Pr)
    shakhov!(SH, SB, u, v, H, B, q, prim, Pr, K)
    shakhov!(S, u, v, w, M, q, prim, Pr)

In-place Shakhov non-equilibrium part

* @arg: particle velocity quadrature points
* @arg: discrete Maxwellian
* @arg: primitive variables, Prandtl number, heat flux, inner degree of freedom

"""
function shakhov!(
    S::T1,
    u::AbstractVector{T2},
    M::T1,
    q,
    prim,
    Pr,
) where {T1<:AbstractArray{<:AbstractFloat,1},T2<:AbstractFloat} # 1F1V

    @. S = @. 0.8 * (1.0 - Pr) * prim[end]^2 / prim[1] *
       (u - prim[2]) *
       q *
       (2.0 * prim[end] * (u - prim[2])^2 - 5.0) *
       M

    return nothing

end

#--- 2F1V ---#
function shakhov!(
    SH::T1,
    SB::T1,
    u::AbstractVector{T2},
    H::T1,
    B::T1,
    q,
    prim,
    Pr,
    K,
) where {T1<:AbstractArray{<:AbstractFloat,1},T2<:AbstractFloat}

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

#--- 1F2V ---#
function shakhov!(
    S::T1,
    u::T2,
    v::T2,
    M::T1,
    q::AbstractVector{T3},
    prim,
    Pr,
) where {T1<:AbstractArray{<:AbstractFloat,2},T2<:AbstractArray{<:AbstractFloat,2},T3<:Real}

    @. S =
        0.8 * (1.0 - Pr) * prim[end]^2 / prim[1] *
        ((u - prim[2]) * q[1] + (v - prim[3]) * q[2]) *
        (2.0 * prim[end] * ((u - prim[2])^2 + (v - prim[3])^2) - 5.0) *
        M

    return nothing

end

#--- 2F2V ---#
function shakhov!(
    SH::T1,
    SB::T1,
    u::T2,
    v::T2,
    H::T1,
    B::T1,
    q::AbstractVector{T3},
    prim,
    Pr,
    K,
) where {T1<:AbstractArray{<:AbstractFloat,2},T2<:AbstractArray{<:AbstractFloat,2},T3<:Real}

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

#--- 1F3V ---#
function shakhov!(
    S::T1,
    u::T2,
    v::T2,
    w::T2,
    M::T1,
    q::AbstractVector{T3},
    prim,
    Pr,
) where {T1<:AbstractArray{<:AbstractFloat,3},T2<:AbstractArray{<:AbstractFloat,3},T3<:Real}

    @. S =
        0.8 * (1.0 - Pr) * prim[end]^2 / prim[1] *
        ((u - prim[2]) * q[1] + (v - prim[3]) * q[2] + (w - prim[4]) * q[3]) *
        (2.0 * prim[end] * ((u - prim[2])^2 + (v - prim[3])^2 + (w - prim[4])^2) - 5.0) *
        M

    return nothing

end
