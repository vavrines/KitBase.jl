"""
$(SIGNATURES)

Shakhov non-equilibrium part

1F1V
"""
function shakhov(u::AV, M::AV, q, prim::AV, Pr)
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
function shakhov(
    u::AV,
    H::T,
    B::T,
    q,
    prim::AV,
    Pr,
    K,
) where {T<:AV}

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
function shakhov(
    u::T,
    v::T,
    M::AM,
    q::AV,
    prim::AV,
    Pr,
) where {T<:AM}

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
function shakhov(
    u::T,
    v::T,
    H::X,
    B::X,
    q::AV,
    prim::AV,
    Pr,
    K,
) where {T<:AM,X<:AM}

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
    M::AA{T1,3},
    q::AV,
    prim::AV,
    Pr,
) where {T<:AA{T2,3},T1} where {T2}

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
function shakhov!(
    SH::T1,
    SB::T1,
    u::AV,
    H::T1,
    B::T1,
    q,
    prim,
    Pr,
    K,
) where {T1<:AV}

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
    q::AV,
    prim,
    Pr,
) where {T1<:AM,T2<:AM}

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
