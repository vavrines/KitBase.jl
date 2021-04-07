"""
    shakhov(
        u::X,
        M::Y,
        q,
        prim::Z,
        Pr,
    ) where {
        X<:AbstractArray{<:AbstractFloat,1},
        Y<:AbstractArray{<:AbstractFloat,1},
        Z<:AbstractArray{<:Real,1},
    }

    shakhov(
        u::T,
        H::X,
        B::X,
        q,
        prim::Y,
        Pr,
        K,
    ) where {
        T<:AbstractArray{<:AbstractFloat,1},
        X<:AbstractArray{<:AbstractFloat,1},
        Y<:AbstractArray{<:Real,1},
    }

    shakhov(
        u::T,
        v::T,
        M::X,
        q::Y,
        prim::Z,
        Pr,
    ) where {
        T<:AbstractArray{<:AbstractFloat,2},
        X<:AbstractArray{<:AbstractFloat,2},
        Y<:AbstractArray{<:AbstractFloat,1},
        Z<:AbstractArray{<:Real,1},
    }

    shakhov(
        u::T,
        v::T,
        H::X,
        B::X,
        q::Y,
        prim::Z,
        Pr,
        K,
    ) where {
        T<:AbstractArray{<:AbstractFloat,2},
        X<:AbstractArray{<:AbstractFloat,2},
        Y<:AbstractArray{<:Real,1},
        Z<:AbstractArray{<:Real,1},
    }

    shakhov(
        u::T,
        v::T,
        w::T,
        M::X,
        q::Y,
        prim::Z,
        Pr,
    ) where {
        T<:AbstractArray{<:AbstractFloat,3},
        X<:AbstractArray{<:AbstractFloat,3},
        Y<:AbstractArray{<:Real,1},
        Z<:AbstractArray{<:Real,1},
    }

Shakhov non-equilibrium part

* @arg: particle velocity quadrature points
* @arg: discrete Maxwellian
* @arg: primitive variables, Prandtl number, heat flux, inner degree of freedom

"""
function shakhov(
    u::X,
    M::Y,
    q,
    prim::Z,
    Pr,
) where {
    X<:AbstractArray{<:AbstractFloat,1},
    Y<:AbstractArray{<:AbstractFloat,1},
    Z<:AbstractArray{<:Real,1},
} # 1F1V

    M_plus = @. 0.8 * (1.0 - Pr) * prim[end]^2 / prim[1] *
       (u - prim[2]) *
       q *
       (2.0 * prim[end] * (u - prim[2])^2 - 5.0) *
       M

    return M_plus

end

#--- 2F1V ---#
function shakhov(
    u::T,
    H::X,
    B::X,
    q,
    prim::Y,
    Pr,
    K,
) where {
    T<:AbstractArray{<:AbstractFloat,1},
    X<:AbstractArray{<:AbstractFloat,1},
    Y<:AbstractArray{<:Real,1},
}

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
    M::X,
    q::Y,
    prim::Z,
    Pr,
) where {
    T<:AbstractArray{<:AbstractFloat,2},
    X<:AbstractArray{<:AbstractFloat,2},
    Y<:AbstractArray{<:AbstractFloat,1},
    Z<:AbstractArray{<:Real,1},
}

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
    q::Y,
    prim::Z,
    Pr,
    K,
) where {
    T<:AbstractArray{<:AbstractFloat,2},
    X<:AbstractArray{<:AbstractFloat,2},
    Y<:AbstractArray{<:Real,1},
    Z<:AbstractArray{<:Real,1},
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
    M::X,
    q::Y,
    prim::Z,
    Pr,
) where {
    T<:AbstractArray{<:AbstractFloat,3},
    X<:AbstractArray{<:AbstractFloat,3},
    Y<:AbstractArray{<:Real,1},
    Z<:AbstractArray{<:Real,1},
}

    M_plus = @. 0.8 * (1.0 - Pr) * prim[end]^2 / prim[1] *
       ((u - prim[2]) * q[1] + (v - prim[3]) * q[2] + (w - prim[4]) * q[3]) *
       (2.0 * prim[end] * ((u - prim[2])^2 + (v - prim[3])^2 + (w - prim[4])^2) - 5.0) *
       M

    return M_plus

end


"""
shakhov!(
    S::T1,
    u::T2,
    M::T1,
    q,
    prim,
    Pr,
) where {T1<:AbstractArray{<:AbstractFloat,1},T2<:AbstractArray{<:AbstractFloat,1}}

shakhov!(
    SH::T1,
    SB::T1,
    u::T2,
    H::T1,
    B::T1,
    q,
    prim,
    Pr,
    K,
) where {T1<:AbstractArray{<:AbstractFloat,1},T2<:AbstractArray{<:AbstractFloat,1}}

shakhov!(
    S::T1,
    u::T2,
    v::T2,
    M::T1,
    q::T3,
    prim,
    Pr,
) where {
    T1<:AbstractArray{<:AbstractFloat,2},
    T2<:AbstractArray{<:AbstractFloat,2},
    T3<:AbstractArray{<:Real,1},
}

shakhov!(
    SH::T1,
    SB::T1,
    u::T2,
    v::T2,
    H::T1,
    B::T1,
    q::T3,
    prim,
    Pr,
    K,
) where {
    T1<:AbstractArray{<:AbstractFloat,2},
    T2<:AbstractArray{<:AbstractFloat,2},
    T3<:AbstractArray{<:Real,1},
}

shakhov!(
    S::T1,
    u::T2,
    v::T2,
    w::T2,
    M::T1,
    q::T3,
    prim,
    Pr,
) where {
    T1<:AbstractArray{<:AbstractFloat,3},
    T2<:AbstractArray{<:AbstractFloat,3},
    T3<:AbstractArray{<:Real,1},
}

In-place Shakhov non-equilibrium part

* @arg: particle velocity quadrature points
* @arg: discrete Maxwellian
* @arg: primitive variables, Prandtl number, heat flux, inner degree of freedom

"""
function shakhov!(
    S::T1,
    u::T2,
    M::T1,
    q,
    prim,
    Pr,
) where {T1<:AbstractArray{<:AbstractFloat,1},T2<:AbstractArray{<:AbstractFloat,1}} # 1F1V

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
    u::T2,
    H::T1,
    B::T1,
    q,
    prim,
    Pr,
    K,
) where {T1<:AbstractArray{<:AbstractFloat,1},T2<:AbstractArray{<:AbstractFloat,1}}

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
    q::T3,
    prim,
    Pr,
) where {
    T1<:AbstractArray{<:AbstractFloat,2},
    T2<:AbstractArray{<:AbstractFloat,2},
    T3<:AbstractArray{<:Real,1},
}

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
    q::T3,
    prim,
    Pr,
    K,
) where {
    T1<:AbstractArray{<:AbstractFloat,2},
    T2<:AbstractArray{<:AbstractFloat,2},
    T3<:AbstractArray{<:Real,1},
}

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
    q::T3,
    prim,
    Pr,
) where {
    T1<:AbstractArray{<:AbstractFloat,3},
    T2<:AbstractArray{<:AbstractFloat,3},
    T3<:AbstractArray{<:Real,1},
}

    @. S =
        0.8 * (1.0 - Pr) * prim[end]^2 / prim[1] *
        ((u - prim[2]) * q[1] + (v - prim[3]) * q[2] + (w - prim[4]) * q[3]) *
        (2.0 * prim[end] * ((u - prim[2])^2 + (v - prim[3])^2 + (w - prim[4])^2) - 5.0) *
        M

    return nothing

end
