"""
Calculate moments of Gaussian distribution `G = (λ / π)^(D / 2) * exp[-λ(c^2 + ξ^2)]`

* internality: `gauss_moments(prim::T) where {T<:AA{<:Real,1}}`
* no internality: `gauss_moments(prim::T, inK) where {T<:AA{<:Real,1}}`

"""
function gauss_moments(prim::T) where {T<:AA{<:Real,1}}

    if eltype(prim) <: Int
        MuL = OffsetArray(similar(prim, Float64, 7), 0:6)
    else
        MuL = OffsetArray(similar(prim, 7), 0:6)
    end
    MuR = similar(MuL)
    Mu = similar(MuL)

    MuL[0] = 0.5 * SpecialFunctions.erfc(-sqrt(prim[end]) * prim[2])
    MuL[1] = prim[2] * MuL[0] + 0.5 * exp(-prim[end] * prim[2]^2) / sqrt(π * prim[end])
    MuR[0] = 0.5 * SpecialFunctions.erfc(sqrt(prim[end]) * prim[2])
    MuR[1] = prim[2] * MuR[0] - 0.5 * exp(-prim[end] * prim[2]^2) / sqrt(π * prim[end])
    for i = 2:6
        MuL[i] = prim[2] * MuL[i-1] + 0.5 * (i - 1) * MuL[i-2] / prim[end]
        MuR[i] = prim[2] * MuR[i-1] + 0.5 * (i - 1) * MuR[i-2] / prim[end]
    end
    @. Mu = MuL + MuR

    if length(prim) == 3

        return Mu, MuL, MuR

    elseif length(prim) == 4

        Mv = similar(MuL)
        Mv[0] = 1.0
        Mv[1] = prim[3]
        for i = 2:6
            Mv[i] = prim[3] * Mv[i-1] + 0.5 * (i - 1) * Mv[i-2] / prim[end]
        end

        return Mu, Mv, MuL, MuR

    elseif length(prim) == 5

        Mv = similar(MuL)
        Mv[0] = 1.0
        Mv[1] = prim[3]
        for i = 2:6
            Mv[i] = prim[3] * Mv[i-1] + 0.5 * (i - 1) * Mv[i-2] / prim[end]
        end

        Mw = similar(MuL)
        Mw[0] = 1.0
        Mw[1] = prim[4]
        for i = 2:6
            Mw[i] = prim[4] * Mw[i-1] + 0.5 * (i - 1) * Mw[i-2] / prim[end]
        end

        return Mu, Mv, Mw, MuL, MuR

    end

end

# ------------------------------------------------------------
# A more general function
# deal with absent internality by setting inK = 0
# ------------------------------------------------------------
function gauss_moments(prim::T, inK) where {T<:AA{<:Real,1}}

    if eltype(prim) <: Int
        MuL = OffsetArray(similar(prim, Float64, 7), 0:6)
    else
        MuL = OffsetArray(similar(prim, 7), 0:6)
    end
    MuR = similar(MuL)
    Mu = similar(MuL)

    MuL[0] = 0.5 * SpecialFunctions.erfc(-sqrt(prim[end]) * prim[2])
    MuL[1] = prim[2] * MuL[0] + 0.5 * exp(-prim[end] * prim[2]^2) / sqrt(π * prim[end])
    MuR[0] = 0.5 * SpecialFunctions.erfc(sqrt(prim[end]) * prim[2])
    MuR[1] = prim[2] * MuR[0] - 0.5 * exp(-prim[end] * prim[2]^2) / sqrt(π * prim[end])
    for i = 2:6
        MuL[i] = prim[2] * MuL[i-1] + 0.5 * (i - 1) * MuL[i-2] / prim[end]
        MuR[i] = prim[2] * MuR[i-1] + 0.5 * (i - 1) * MuR[i-2] / prim[end]
    end
    @. Mu = MuL + MuR

    if length(prim) == 3

        Mxi = similar(MuL, 0:2)
        Mxi[0] = 1.0
        Mxi[1] = 0.5 * inK / prim[end]
        Mxi[2] = (inK^2 + 2.0 * inK) / (4.0 * prim[end]^2)

        return Mu, Mxi, MuL, MuR

    elseif length(prim) == 4

        Mv = similar(MuL)
        Mv[0] = 1.0
        Mv[1] = prim[3]
        for i = 2:6
            Mv[i] = prim[3] * Mv[i-1] + 0.5 * (i - 1) * Mv[i-2] / prim[end]
        end

        Mxi = similar(MuL, 0:2)
        Mxi[0] = 1.0
        Mxi[1] = 0.5 * inK / prim[end]
        Mxi[2] = (inK^2 + 2.0 * inK) / (4.0 * prim[end]^2)

        return Mu, Mv, Mxi, MuL, MuR

    elseif length(prim) == 5

        Mv = similar(MuL)
        Mv[0] = 1.0
        Mv[1] = prim[3]
        for i = 2:6
            Mv[i] = prim[3] * Mv[i-1] + 0.5 * (i - 1) * Mv[i-2] / prim[end]
        end

        Mw = similar(MuL)
        Mw[0] = 1.0
        Mw[1] = prim[4]
        for i = 2:6
            Mw[i] = prim[4] * Mw[i-1] + 0.5 * (i - 1) * Mw[i-2] / prim[end]
        end

        return Mu, Mv, Mw, MuL, MuR

    end

end


"""
    moments_conserve(Mu::OffsetArray{<:FN,1}, alpha::Integer)

    moments_conserve(Mu::OffsetArray{<:Real,1}, Mxi::OffsetArray{<:Real,1},
        alpha::Integer, delta::Integer)

    moments_conserve(Mu::OffsetArray{<:Real,1}, Mv::OffsetArray{<:Real,1},
        Mw::OffsetArray{<:Real,1}, alpha::Integer, beta::Integer, delta::Integer)

Calculate conservative moments of particle distribution

"""
moments_conserve(Mu::T, alpha::I) where {T<:OffsetArray{<:FN,1},I<:Integer} =
    Mu[alpha]

function moments_conserve(
    Mu::T,
    Mxi::T,
    alpha::I,
    delta::I,
) where {T<:OffsetArray{<:FN,1},I<:Integer}

    uv = similar(Mu, 3)
    uv[1] = Mu[alpha] * Mxi[delta÷2]
    uv[2] = Mu[alpha+1] * Mxi[delta÷2]
    uv[3] = 0.5 * (Mu[alpha+2] * Mxi[delta÷2] + Mu[alpha] * Mxi[(delta+2)÷2])

    return uv

end

function moments_conserve(
    Mu::T,
    Mv::T,
    Mw::T,
    alpha::I,
    beta::I,
    delta::I,
) where {T<:OffsetArray{<:FN,1},I<:Integer}

    if length(Mw) == 3 # internal motion
        uv = similar(Mu, 4)
        uv[1] = Mu[alpha] * Mv[beta] * Mw[delta÷2]
        uv[2] = Mu[alpha+1] * Mv[beta] * Mw[delta÷2]
        uv[3] = Mu[alpha] * Mv[beta+1] * Mw[delta÷2]
        uv[4] =
            0.5 * (
                Mu[alpha+2] * Mv[beta] * Mw[delta÷2] +
                Mu[alpha] * Mv[beta+2] * Mw[delta÷2] +
                Mu[alpha] * Mv[beta] * Mw[(delta+2)÷2]
            )
    else
        uv = similar(Mu, 5)
        uv[1] = Mu[alpha] * Mv[beta] * Mw[delta]
        uv[2] = Mu[alpha+1] * Mv[beta] * Mw[delta]
        uv[3] = Mu[alpha] * Mv[beta+1] * Mw[delta]
        uv[4] = Mu[alpha] * Mv[beta] * Mw[delta+1]
        uv[5] =
            0.5 * (
                Mu[alpha+2] * Mv[beta] * Mw[delta] +
                Mu[alpha] * Mv[beta+2] * Mw[delta] +
                Mu[alpha] * Mv[beta] * Mw[delta+2]
            )
    end

    return uv

end

# ------------------------------------------------------------
# Discrete moments of conservative variables
# ------------------------------------------------------------
#--- 1F1V ---#
function moments_conserve(
    f::X,
    u::T,
    ω::T,
) where {X<:AA{<:FN,1},T<:AA{<:FN,1}}
    w = similar(f, 3)
    w[1] = discrete_moments(f, u, ω, 0)
    w[2] = discrete_moments(f, u, ω, 1)
    w[3] = 0.5 * discrete_moments(f, u, ω, 2)

    return w
end

#--- 2F1V ---#
function moments_conserve(
    h::X,
    b::X,
    u::T,
    ω::T,
) where {X<:AA{<:FN,1},T<:AA{<:FN,1}}
    w = similar(h, 3)
    w[1] = discrete_moments(h, u, ω, 0)
    w[2] = discrete_moments(h, u, ω, 1)
    w[3] = 0.5 * (discrete_moments(h, u, ω, 2) + discrete_moments(b, u, ω, 0))

    return w
end

#--- 1F2V ---#
function moments_conserve(
    f::X,
    u::T,
    v::T,
    ω::T,
) where {X<:AA{<:FN,2},T<:AA{<:FN,2}}
    w = similar(f, 4)
    w[1] = discrete_moments(f, u, ω, 0)
    w[2] = discrete_moments(f, u, ω, 1)
    w[3] = discrete_moments(f, v, ω, 1)
    w[4] = 0.5 * (discrete_moments(f, u, ω, 2) + discrete_moments(f, v, ω, 2))

    return w
end

#--- 2F2V ---#
function moments_conserve(
    h::X,
    b::X,
    u::T,
    v::T,
    ω::T,
) where {X<:AA{<:FN,2},T<:AA{<:FN,2}}
    w = similar(h, 4)
    w[1] = discrete_moments(h, u, ω, 0)
    w[2] = discrete_moments(h, u, ω, 1)
    w[3] = discrete_moments(h, v, ω, 1)
    w[4] =
        0.5 * (
            discrete_moments(h, u, ω, 2) +
            discrete_moments(h, v, ω, 2) +
            discrete_moments(b, u, ω, 0)
        )

    return w
end

#--- 3F2V ---#
function moments_conserve(
    h0::X,
    h1::X,
    h2::X,
    u::T,
    v::T,
    ω::T,
) where {X<:AA{<:FN,2},T<:AA{<:FN,2}}
    w = similar(h0, 5)
    w[1] = discrete_moments(h0, u, ω, 0)
    w[2] = discrete_moments(h0, u, ω, 1)
    w[3] = discrete_moments(h0, v, ω, 1)
    w[4] = discrete_moments(h1, u, ω, 0)
    w[5] =
        0.5 * (
            discrete_moments(h0, u, ω, 2) +
            discrete_moments(h0, v, ω, 2) +
            discrete_moments(h2, u, ω, 0)
        )

    return w
end

#--- 1F3V ---#
function moments_conserve(
    f::X,
    u::T,
    v::T,
    w::T,
    ω::T,
) where {X<:AA{<:FN,3},T<:AA{<:FN,3}}
    moments = similar(f, 5)

    moments[1] = discrete_moments(f, u, ω, 0)
    moments[2] = discrete_moments(f, u, ω, 1)
    moments[3] = discrete_moments(f, v, ω, 1)
    moments[4] = discrete_moments(f, w, ω, 1)
    moments[5] =
        0.5 * (
            discrete_moments(f, u, ω, 2) +
            discrete_moments(f, v, ω, 2) +
            discrete_moments(f, w, ω, 2)
        )

    return moments
end

#--- 4F1V ---#
function moments_conserve(
    h0::X,
    h1::X,
    h2::X,
    h3::X,
    u::T,
    ω::T,
) where {X<:AA{<:FN,1},T<:AA{<:FN,1}}
    moments = similar(h0, 5)

    moments[1] = discrete_moments(h0, u, ω, 0)
    moments[2] = discrete_moments(h0, u, ω, 1)
    moments[3] = discrete_moments(h1, u, ω, 0)
    moments[4] = discrete_moments(h2, u, ω, 0)
    moments[5] = 0.5 * discrete_moments(h0, u, ω, 2) + 0.5 * discrete_moments(h3, u, ω, 0)

    return moments
end


"""
Calculate conservative moments of diatomic particle distribution

- 1D: `diatomic_moments_conserve(
    h::X,
    b::X,
    r::X,
    u::T,
    ω::T,
) where {X<:AA{<:FN,1},T<:AA{<:FN,1}}`
- 2D: `diatomic_moments_conserve(
    h0::X,
    h1::X,
    h2::X,
    u::T,
    v::T,
    ω::T,
) where {X<:AA{<:FN,2},T<:AA{<:FN,2}}`

"""
function diatomic_moments_conserve(
    h::X,
    b::X,
    r::X,
    u::T,
    ω::T,
) where {X<:AA{<:FN,1},T<:AA{<:FN,1}}
    w = similar(h, 4)
    w[1] = discrete_moments(h, u, ω, 0)
    w[2] = discrete_moments(h, u, ω, 1)
    w[3] =
        0.5 * (
            discrete_moments(h, u, ω, 2) +
            discrete_moments(b, u, ω, 0) +
            discrete_moments(r, u, ω, 0)
        )
    w[4] = 0.5 * discrete_moments(r, u, ω, 0)

    return w
end

#--- 3F2V ---#
function diatomic_moments_conserve(
    h0::X,
    h1::X,
    h2::X,
    u::T,
    v::T,
    ω::T,
) where {X<:AA{<:FN,2},T<:AA{<:FN,2}}
    w = similar(h0, 5)
    w[1] = discrete_moments(h0, u, ω, 0)
    w[2] = discrete_moments(h0, u, ω, 1)
    w[3] = discrete_moments(h0, v, ω, 1)
    w[4] =
        0.5 * (
            discrete_moments(h0, u, ω, 2) +
            discrete_moments(h0, v, ω, 2) +
            discrete_moments(h1, u, ω, 0) +
            discrete_moments(h2, u, ω, 0)
        )
    w[5] = 0.5 * discrete_moments(h2, u, ω, 0)

    return w
end


"""
Calculate slope-related conservative moments
`a = a1 + u * a2 + 0.5 * u^2 * a3`

"""
moments_conserve_slope(a, Mu::T, alpha::I) where {T<:OffsetArray{<:Real,1},I<:Int} =
    a * moments_conserve(Mu, alpha)

moments_conserve_slope(
    a::X,
    Mu::Y,
    Mxi::Y,
    alpha::I,
) where {X<:AA{<:Real,1},Y<:OffsetArray{<:Real,1},I<:Int} =
    a[1] .* moments_conserve(Mu, Mxi, alpha + 0, 0) .+
    a[2] .* moments_conserve(Mu, Mxi, alpha + 1, 0) .+
    0.5 * a[3] .* moments_conserve(Mu, Mxi, alpha + 2, 0) .+
    0.5 * a[3] .* moments_conserve(Mu, Mxi, alpha + 0, 2)

function moments_conserve_slope(
    a::X,
    Mu::Y,
    Mv::Y,
    Mw::Y,
    alpha::I,
    beta::I,
    delta = 0::I,
) where {X<:AA{<:Real,1},Y<:OffsetArray{<:Real,1},I<:Integer}

    if length(a) == 4
        return a[1] .* moments_conserve(Mu, Mv, Mw, alpha + 0, beta + 0, 0) .+
            a[2] .* moments_conserve(Mu, Mv, Mw, alpha + 1, beta + 0, 0) .+
            a[3] .* moments_conserve(Mu, Mv, Mw, alpha + 0, beta + 1, 0) .+
            0.5 * a[4] .* moments_conserve(Mu, Mv, Mw, alpha + 2, beta + 0, 0) .+
            0.5 * a[4] .* moments_conserve(Mu, Mv, Mw, alpha + 0, beta + 2, 0) .+
            0.5 * a[4] .* moments_conserve(Mu, Mv, Mw, alpha + 0, beta + 0, 2)
    elseif length(a) == 5
        return a[1] .* moments_conserve(Mu, Mv, Mw, alpha + 0, beta + 0, delta + 0) .+
            a[2] .* moments_conserve(Mu, Mv, Mw, alpha + 1, beta + 0, delta + 0) .+
            a[3] .* moments_conserve(Mu, Mv, Mw, alpha + 0, beta + 1, delta + 0) .+
            a[4] .* moments_conserve(Mu, Mv, Mw, alpha + 0, beta + 0, delta + 1) .+
            0.5 * a[5] .* moments_conserve(Mu, Mv, Mw, alpha + 2, beta + 0, delta + 0) .+
            0.5 * a[5] .* moments_conserve(Mu, Mv, Mw, alpha + 0, beta + 2, delta + 0) .+
            0.5 * a[5] .* moments_conserve(Mu, Mv, Mw, alpha + 0, beta + 0, delta + 2)
    end

end


"""
    * direct quadrature
    discrete_moments(f, ω)

    * velocity moments
    discrete_moments(f, u, ω, n)

Discrete moments of particle distribution

"""
discrete_moments(f, ω) = sum(@. ω * f)

discrete_moments(f, u, ω, n) = sum(@. ω * u^n * f)


"""
    pressure(f, prim, u, ω)
    pressure(h, b, prim, u, ω, K)
    pressure(h, b, prim, u, v, ω, K)

Calculate pressure from particle distribution function

"""
pressure(f, prim, u, ω) = sum(@. ω * (u - prim[2])^2 * f)

pressure(h, b, prim, u, ω, K) =
    (sum(@. ω * (u - prim[2])^2 * h) + sum(@. ω * b)) / (K + 1.0)

pressure(h, b, prim, u, v, ω, K) =
    (sum(@. ω * ((u - prim[2])^2 + (v - prim[3])^2) * h) + sum(@. ω * b)) / (K + 2.0)


"""
Calculate stress tensor from particle distribution function

"""
stress(f, prim, u, ω) = sum(@. ω * (u - prim[2]) * (u - prim[2]) * f)

function stress(f, prim, u, v, ω)
    P = similar(prim, 2, 2)

    P[1, 1] = sum(@. ω * (u - prim[2]) * (u - prim[2]) * f)
    P[1, 2] = sum(@. ω * (u - prim[2]) * (v - prim[3]) * f)
    P[2, 1] = P[1, 2]
    P[1, 2] = sum(@. ω * (v - prim[3]) * (v - prim[3]) * f)

    return P
end


"""
    heat_flux(f, prim, u, ω)
    heat_flux(h, b, prim, u, ω)
    heat_flux(h, b, r, prim, u, ω)
    heat_flux(f, prim, u, v, ω)
    heat_flux(h, b, prim, u, v, ω)
    heat_flux(h, b, r, prim, u, v, ω)
    heat_flux(f, prim, u, v, w, ω)

Calculate heat flux from particle distribution function

Multiple dispatch doesn't consider unstructured multi-dimensional velocity space.
In that case a new method needs to be defined.

"""
heat_flux(f, prim, u, ω) = 0.5 * sum(@. ω * (u - prim[2]) * (u - prim[2])^2 * f) # 1F1V

#--- 2F1V ---#
heat_flux(
    h::X,
    b::X,
    prim::Y,
    u::Z,
    ω::Z,
) where {
    X<:AA{<:FN,1},
    Y<:AA{<:Real,1},
    Z<:AA{<:FN,1},
} = 0.5 * (sum(@. ω * (u - prim[2]) * (u - prim[2])^2 * h) + sum(@. ω * (u - prim[2]) * b))

#--- 3F1V Rykov ---#
function heat_flux(
    h::X,
    b::X,
    r::X,
    prim::Y,
    u::Z,
    ω::Z,
) where {
    X<:AA{<:FN,1},
    Y<:AA{<:Real,1},
    Z<:AA{<:FN,1},
}

    q = similar(h, 2)

    q[1] =
        0.5 *
        (sum(@. ω * (u - prim[2]) * (u - prim[2])^2 * h) + sum(@. ω * (u - prim[2]) * b))
    q[2] = 0.5 * (sum(@. ω * (u - prim[2]) * r))

    return q

end

#--- 1F2V ---#
function heat_flux(
    h::X,
    prim::Y,
    u::Z,
    v::Z,
    ω::Z,
) where {
    X<:AA{<:FN,2},
    Y<:AA{<:Real,1},
    Z<:AA{<:FN,2},
}

    q = similar(h, 2)
    q[1] = 0.5 * sum(@. ω * (u - prim[2]) * ((u - prim[2])^2 + (v - prim[3])^2) * h)
    q[2] = 0.5 * sum(@. ω * (v - prim[3]) * ((u - prim[2])^2 + (v - prim[3])^2) * h)

    return q

end

#--- 2F2V ---#
function heat_flux(
    h::X,
    b::X,
    prim::Y,
    u::Z,
    v::Z,
    ω::Z,
) where {
    X<:AA{<:FN,2},
    Y<:AA{<:Real,1},
    Z<:AA{<:FN,2},
}

    q = similar(h, 2)

    q[1] =
        0.5 * (
            sum(@. ω * (u - prim[2]) * ((u - prim[2])^2 + (v - prim[3])^2) * h) +
            sum(@. ω * (u - prim[2]) * b)
        )
    q[2] =
        0.5 * (
            sum(@. ω * (v - prim[3]) * ((u - prim[2])^2 + (v - prim[3])^2) * h) +
            sum(@. ω * (v - prim[3]) * b)
        )

    return q

end

#--- 3F2V Rykov ---#
function heat_flux(
    h::X,
    b::X,
    r::X,
    prim::Y,
    u::Z,
    v::Z,
    ω::Z,
) where {
    X<:AA{<:FN,2},
    Y<:AA{<:Real,1},
    Z<:AA{<:FN,2},
}

    q = similar(h, 4)

    q[1] =
        0.5 * (
            sum(@. ω * (u - prim[2]) * ((u - prim[2])^2 + (v - prim[3])^2) * h) +
            sum(@. ω * (u - prim[2]) * b)
        )
    q[2] =
        0.5 * (
            sum(@. ω * (v - prim[3]) * ((u - prim[2])^2 + (v - prim[3])^2) * h) +
            sum(@. ω * (v - prim[3]) * b)
        )
    q[3] = 0.5 * sum(@. ω * (u - prim[2]) * r)
    q[4] = 0.5 * sum(@. ω * (v - prim[3]) * r)

    return q

end

#--- 1F3V ---#
function heat_flux(
    f::X,
    prim::Y,
    u::Z,
    v::Z,
    w::Z,
    ω::Z,
) where {
    X<:AA{<:FN,3},
    Y<:AA{<:Real,1},
    Z<:AA{<:FN,3},
}

    q = similar(f, 3)

    q[1] =
        0.5 * sum(
            @. ω *
               (u - prim[2]) *
               ((u - prim[2])^2 + (v - prim[3])^2 + (w - prim[4])^2) *
               f
        )
    q[2] =
        0.5 * sum(
            @. ω *
               (v - prim[3]) *
               ((u - prim[2])^2 + (v - prim[3])^2 + (w - prim[4])^2) *
               f
        )
    q[3] =
        0.5 * sum(
            @. ω *
               (w - prim[4]) *
               ((u - prim[2])^2 + (v - prim[3])^2 + (w - prim[4])^2) *
               f
        )

    return q

end
