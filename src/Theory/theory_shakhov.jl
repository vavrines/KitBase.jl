"""
$(SIGNATURES)

Shakhov non-equilibrium part
"""
function shakhov(u::AV{X}, M::AV{Y}, q, prim::AV{Z}, Pr) where {X<:FN,Y<:FN,Z<:Real} # 1F1V
    M_plus = @. 0.8 * (1.0 - Pr) * prim[end]^2 / prim[1] *
       (u - prim[2]) *
       q *
       (2.0 * prim[end] * (u - prim[2])^2 - 5.0) *
       M

    return M_plus
end

#--- 2F1V ---#
function shakhov(
    u::AV{X},
    H::Y,
    B::Y,
    q,
    prim::AV{Z},
    Pr,
    K,
) where {X<:FN,Y<:AV{<:FN},Z<:Real}

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
    M::AM{X},
    q::AV{Y},
    prim::AV{Z},
    Pr,
) where {T<:AA{<:FN,2},X<:FN,Y<:FN,Z<:Real}

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
    q::AV{Y},
    prim::AV{Z},
    Pr,
    K,
) where {T<:AA{<:FN,2},X<:AA{<:FN,2},Y<:Real,Z<:Real}

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
    M::AA{X,3},
    q::AV{Y},
    prim::AV{Z},
    Pr,
) where {T<:AA{<:FN,3},X<:FN,Y<:Real,Z<:Real}

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
function shakhov!(S::T1, u::AV{T2}, M::T1, q, prim, Pr) where {T1<:AA{<:FN,1},T2<:FN}
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
function shakhov!(
    SH::T1,
    SB::T1,
    u::AV{T2},
    H::T1,
    B::T1,
    q,
    prim,
    Pr,
    K,
) where {T1<:AA{<:FN,1},T2<:FN}

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
function shakhov!(
    S::T1,
    u::T2,
    v::T2,
    M::T1,
    q::AV{T3},
    prim,
    Pr,
) where {T1<:AA{<:FN,2},T2<:AA{<:FN,2},T3<:Real}

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
    q::AV{T3},
    prim,
    Pr,
    K,
) where {T1<:AA{<:FN,2},T2<:AA{<:FN,2},T3<:Real}

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
    q::AV{T3},
    prim,
    Pr,
) where {T1<:AA{<:FN,3},T2<:AA{<:FN,3},T3<:Real}

    @. S =
        0.8 * (1.0 - Pr) * prim[end]^2 / prim[1] *
        ((u - prim[2]) * q[1] + (v - prim[3]) * q[2] + (w - prim[4]) * q[3]) *
        (2.0 * prim[end] * ((u - prim[2])^2 + (v - prim[3])^2 + (w - prim[4])^2) - 5.0) *
        M

    return nothing

end
