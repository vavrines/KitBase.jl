# ------------------------------------------------------------
# Single-component gas
# ------------------------------------------------------------

"""
$(SIGNATURES)

Calculate moments of Gaussian distribution `G = (λ / π)^(D / 2) * exp[-λ(c^2 + ξ^2)]`
"""
function gauss_moments(prim)
    if eltype(prim) <: Integer
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

"""
$(SIGNATURES)

Calculate moments of Gaussian distribution with internal energy
"""
function gauss_moments(prim, inK)
    if length(prim) == 3
        Mu, MuL, MuR = gauss_moments(prim)

        Mxi = similar(MuL, 0:2)
        Mxi[0] = 1.0
        Mxi[1] = 0.5 * inK / prim[end]
        Mxi[2] = (inK^2 + 2.0 * inK) / (4.0 * prim[end]^2)

        return Mu, Mxi, MuL, MuR
    elseif length(prim) == 4
        Mu, Mv, MuL, MuR = gauss_moments(prim)

        Mxi = similar(MuL, 0:2)
        Mxi[0] = 1.0
        Mxi[1] = 0.5 * inK / prim[end]
        Mxi[2] = (inK^2 + 2.0 * inK) / (4.0 * prim[end]^2)

        return Mu, Mv, Mxi, MuL, MuR
    elseif length(prim) == 5
        return gauss_moments(prim)
    end
end


"""
$(SIGNATURES)

Calculate conservative moments of particle distribution
"""
moments_conserve(Mu::OffsetVector, alpha::Integer) = Mu[alpha]

"""
$(SIGNATURES)
"""
function moments_conserve(
    Mu::T,
    Mxi::T,
    alpha::Integer,
    delta::Integer,
) where {T<:OffsetVector}

    uv = similar(Mu, 3)
    uv[1] = Mu[alpha] * Mxi[delta÷2]
    uv[2] = Mu[alpha+1] * Mxi[delta÷2]
    uv[3] = 0.5 * (Mu[alpha+2] * Mxi[delta÷2] + Mu[alpha] * Mxi[(delta+2)÷2])

    return uv

end

"""
$(SIGNATURES)
"""
function moments_conserve(
    Mu::T,
    Mv::T,
    Mw::T,
    alpha::Integer,
    beta::Integer,
    delta::Integer,
) where {T<:OffsetVector}

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

"""
$(SIGNATURES)

Discrete moments of conservative variables

1F1V
"""
function moments_conserve(f, u, ω, ::Type{VDF{1,1}})
    w = similar(f, 3)
    w[1] = discrete_moments(f, u, ω, 0)
    w[2] = discrete_moments(f, u, ω, 1)
    w[3] = 0.5 * discrete_moments(f, u, ω, 2)

    return w
end

"""
$(SIGNATURES)

1F2V
"""
function moments_conserve(f, u, v, ω, ::Type{VDF{1,2}})
    w = similar(f, 4)
    w[1] = discrete_moments(f, u, ω, 0)
    w[2] = discrete_moments(f, u, ω, 1)
    w[3] = discrete_moments(f, v, ω, 1)
    w[4] = 0.5 * (discrete_moments(f, u, ω, 2) + discrete_moments(f, v, ω, 2))

    return w
end

"""
$(SIGNATURES)

2F1V
"""
function moments_conserve(h, b, u, ω, ::Type{VDF{2,1}})
    w = similar(h, 3)
    w[1] = discrete_moments(h, u, ω, 0)
    w[2] = discrete_moments(h, u, ω, 1)
    w[3] = 0.5 * (discrete_moments(h, u, ω, 2) + discrete_moments(b, u, ω, 0))

    return w
end

"""
$(SIGNATURES)

2F2V
"""
function moments_conserve(h, b, u, v, ω, ::Type{VDF{2,2}})
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

"""
$(SIGNATURES)

3F2V
"""
function moments_conserve(h0, h1, h2, u, v, ω, ::Type{VDF{3,2}})
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

"""
$(SIGNATURES)

1F3V
"""
function moments_conserve(f, u, v, w, ω, ::Type{VDF{1,3}})
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

"""
$(SIGNATURES)

4F1V
"""
function moments_conserve(h0, h1, h2, h3, u, ω, ::Type{VDF{4,1}})
    moments = similar(h0, 5)

    moments[1] = discrete_moments(h0, u, ω, 0)
    moments[2] = discrete_moments(h0, u, ω, 1)
    moments[3] = discrete_moments(h1, u, ω, 0)
    moments[4] = discrete_moments(h2, u, ω, 0)
    moments[5] = 0.5 * discrete_moments(h0, u, ω, 2) + 0.5 * discrete_moments(h3, u, ω, 0)

    return moments
end

"""
$(SIGNATURES)

Shortcut methods

1F1V
"""
function moments_conserve(f::AV, u::T, ω::T) where {T<:AV}
    return moments_conserve(f, u, ω, VDF{1,1})
end

"""
$(SIGNATURES)

2F1V & 1F2V
"""
function moments_conserve(a1::AA, a2::AA, a3::T, a4::T) where {T<:AA}
    if minimum(a2) < -0.9
        f, u, v, ω = a1, a2, a3, a4
        return moments_conserve(f, u, v, ω, VDF{1,2})
    else
        h, b, u, ω = a1, a2, a3, a4
        return moments_conserve(h, b, u, ω, VDF{2,1})
    end
end

"""
$(SIGNATURES)

2F2V & 1F3V
"""
function moments_conserve(a1::AA, a2::AA, a3::T, a4::T, a5::T) where {T<:AA}
    if minimum(a2) < -0.9
        f, u, v, w, ω = a1, a2, a3, a4, a5
        return moments_conserve(f, u, v, w, ω, VDF{1,3})
    else
        h, b, u, v, ω = a1, a2, a3, a4, a5
        return moments_conserve(h, b, u, v, ω, VDF{2,2})
    end
end

"""
$(SIGNATURES)

4F1V & 3F2V
"""
function moments_conserve(a1::AA, a2::AA, a3::AA, a4::AA, a5::T, a6::T) where {T<:AA}
    if minimum(a4) < -0.9
        h0, h1, h2, u, v, ω = a1, a2, a3, a4, a5, a6
        return moments_conserve(h0, h1, h2, u, v, ω, VDF{3,2})
    else
        h0, h1, h2, h3, u, ω = a1, a2, a3, a4, a5, a6
        return moments_conserve(h0, h1, h2, h3, u, ω, VDF{4,1})
    end
end


"""
$(SIGNATURES)

Calculate slope-related conservative moments
`a = a1 + u * a2 + 0.5 * u^2 * a3`

"""
moments_conserve_slope(a, Mu::OffsetVector, alpha::Integer) =
    a * moments_conserve(Mu, alpha)

"""
$(SIGNATURES)
"""
moments_conserve_slope(a::AV, Mu::Y, Mxi::Y, alpha::Integer) where {Y} =
    a[1] .* moments_conserve(Mu, Mxi, alpha + 0, 0) .+
    a[2] .* moments_conserve(Mu, Mxi, alpha + 1, 0) .+
    0.5 * a[3] .* moments_conserve(Mu, Mxi, alpha + 2, 0) .+
    0.5 * a[3] .* moments_conserve(Mu, Mxi, alpha + 0, 2)

"""
$(SIGNATURES)
"""
function moments_conserve_slope(
    a::AV,
    Mu::Y,
    Mv::Y,
    Mw::Y,
    alpha::Integer,
    beta::Integer,
    delta = 0::Integer,
) where {Y}

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
$(SIGNATURES)

Calculate conservative moments from microscopic moments
"""
function flux_conserve!(fw, args...)
    if length(fw) == 3
        return flux_conserve_1d!(fw, args...)
    elseif length(fw) == 4
        return flux_conserve_2d!(fw, args...)
    elseif length(fw) == 5
        return flux_conserve_3d!(fw, args...)
    end

    return nothing
end

#--- 1f1v ---#
function flux_conserve_1d!(fw, ff, u, ω)
    fw[1] = sum(ω .* ff)
    fw[2] = sum(u .* ω .* ff)
    fw[end] = 0.5 * sum(u .^ 2 .* ω .* ff)
end

#--- 2f1v ---#
function flux_conserve_1d!(fw, fh, fb, u, ω)
    fw[1] = sum(ω .* fh)
    fw[2] = sum(u .* ω .* fh)
    fw[end] = 0.5 * (sum(u .^ 2 .* ω .* fh) + sum(ω .* fb))
end

#--- 1f2v ---#
function flux_conserve_2d!(fw, ff, u, v, ω)
    fw[1] = sum(ω .* ff)
    fw[2] = sum(u .* ω .* ff)
    fw[3] = sum(v .* ω .* ff)
    fw[end] = 0.5 * sum((u .^ 2 .+ v .^ 2) .* ω .* ff)
end

#--- 2f2v ---#
function flux_conserve_2d!(fw, fh, fb, u, v, ω)
    fw[1] = sum(ω .* fh)
    fw[2] = sum(u .* ω .* fh)
    fw[3] = sum(v .* ω .* fh)
    fw[end] = 0.5 * (sum((u .^ 2 .+ v .^ 2) .* ω .* fh) + sum(ω .* fb))
end

#--- 1f3v ---#
function flux_conserve_3d!(fw, ff, u, v, w, ω)
    fw[1] = sum(ω .* ff)
    fw[2] = sum(u .* ω .* ff)
    fw[3] = sum(v .* ω .* ff)
    fw[4] = sum(w .* ω .* ff)
    fw[end] = 0.5 * sum((u .^ 2 .+ v .^ 2 + w .^ 2) .* ω .* ff)
end


"""
$(SIGNATURES)

Discrete moments of particle distribution
"""
discrete_moments(f, ω) = sum(@. ω * f)

"""
$(SIGNATURES)
"""
discrete_moments(f, u, ω, n) = sum(@. ω * u^n * f)


"""
$(SIGNATURES)

Calculate pressure `p=nkT`
"""
(pressure(prim::AV{T})::T) where {T} = 0.5 * prim[1] / prim[end]

"""
$(SIGNATURES)

Calculate pressure from particle distribution function
"""
pressure(f, prim, u, ω) = sum(@. ω * (u - prim[2])^2 * f)

"""
$(SIGNATURES)
"""
pressure(h, b, prim, u, ω, K) =
    (sum(@. ω * (u - prim[2])^2 * h) + sum(@. ω * b)) / (K + 1.0)

"""
$(SIGNATURES)
"""
pressure(h, b, prim, u, v, ω, K, ::Type{VDF{2,2}}) =
    (sum(@. ω * ((u - prim[2])^2 + (v - prim[3])^2) * h) + sum(@. ω * b)) / (K + 2.0)

"""
$(SIGNATURES)
"""
pressure(f, prim, u, v, w, ω, K, ::Type{VDF{1,3}}) =
    sum(@. ω * ((u - prim[2])^2 + (v - prim[3])^2 + (w - prim[4])^2) * f) / (K + 3.0)


"""
$(SIGNATURES)

Calculate stress tensor from particle distribution function
"""
stress(f, prim, u, ω) = sum(@. ω * (u - prim[2]) * (u - prim[2]) * f)

"""
$(SIGNATURES)
"""
function stress(f, prim, u, v, ω)
    P = similar(prim, 2, 2)

    P[1, 1] = sum(@. ω * (u - prim[2]) * (u - prim[2]) * f)
    P[1, 2] = sum(@. ω * (u - prim[2]) * (v - prim[3]) * f)
    P[2, 1] = P[1, 2]
    P[2, 2] = sum(@. ω * (v - prim[3]) * (v - prim[3]) * f)

    return P
end

"""
$(SIGNATURES)
"""
function stress(f, prim, u, v, w, ω)
    P = similar(prim, 3, 3)

    P[1, 1] = sum(@. ω * (u - prim[2]) * (u - prim[2]) * f)
    P[1, 2] = sum(@. ω * (u - prim[2]) * (v - prim[3]) * f)
    P[1, 3] = sum(@. ω * (u - prim[2]) * (w - prim[4]) * f)

    P[2, 1] = P[1, 2]
    P[2, 2] = sum(@. ω * (v - prim[3]) * (v - prim[3]) * f)
    P[2, 3] = sum(@. ω * (v - prim[3]) * (w - prim[4]) * f)

    P[3, 1] = P[1, 3]
    P[3, 2] = P[2, 3]
    P[3, 3] = sum(@. ω * (w - prim[4]) * (w - prim[4]) * f)

    return P
end


"""
$(SIGNATURES)

Calculate heat flux from particle distribution function

Multiple dispatch doesn't consider unstructured multi-dimensional velocity space.
In that case a new method needs to be defined.

1F1V
"""
heat_flux(f, prim, u, ω) = 0.5 * sum(@. ω * (u - prim[2]) * (u - prim[2])^2 * f) # 1F1V

"""
$(SIGNATURES)

2F1V
"""
heat_flux(h::X, b::X, prim::AV, u::Z, ω::Z) where {X<:AV,Z<:AV} =
    0.5 * (sum(@. ω * (u - prim[2]) * (u - prim[2])^2 * h) + sum(@. ω * (u - prim[2]) * b))

"""
$(SIGNATURES)

3F1V Rykov model
"""
function heat_flux(h::X, b::X, r::X, prim::AV, u::Z, ω::Z) where {X<:AV,Z<:AV}

    q = similar(h, 2)

    q[1] =
        0.5 *
        (sum(@. ω * (u - prim[2]) * (u - prim[2])^2 * h) + sum(@. ω * (u - prim[2]) * b))
    q[2] = 0.5 * (sum(@. ω * (u - prim[2]) * r))

    return q

end

"""
$(SIGNATURES)

1F2V
"""
function heat_flux(h::AM, prim::AV, u::Z, v::Z, ω::Z) where {Z<:AM}

    q = similar(h, 2)
    q[1] = 0.5 * sum(@. ω * (u - prim[2]) * ((u - prim[2])^2 + (v - prim[3])^2) * h)
    q[2] = 0.5 * sum(@. ω * (v - prim[3]) * ((u - prim[2])^2 + (v - prim[3])^2) * h)

    return q

end

"""
$(SIGNATURES)

2F2V
"""
function heat_flux(h::X, b::X, prim::AV, u::Z, v::Z, ω::Z) where {X<:AM,Z<:AM}

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

"""
$(SIGNATURES)

3F2V Rykov model
"""
function heat_flux(h::X, b::X, r::X, prim::AV, u::Z, v::Z, ω::Z) where {X<:AM,Z<:AM}

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

"""
$(SIGNATURES)

1F3V
"""
function heat_flux(
    f::AA{T,3},
    prim::AV,
    u::Z,
    v::Z,
    w::Z,
    ω::Z,
) where {T,Z<:AA{T1,3}} where {T1}

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

# ------------------------------------------------------------
# Polyatomic gas
# ------------------------------------------------------------

"""
$(SIGNATURES)

Calculate conservative moments of polyatomic distribution function
"""
function polyatomic_moments_conserve!(w, h::X, b::X, r::X, u::T, ω::T) where {X<:AV,T<:AV}
    w[1] = discrete_moments(h, u, ω, 0)
    w[2] = discrete_moments(h, u, ω, 1)
    w[3] =
        0.5 * (
            discrete_moments(h, u, ω, 2) +
            discrete_moments(b, u, ω, 0) +
            discrete_moments(r, u, ω, 0)
        )
    w[4] = 0.5 * discrete_moments(r, u, ω, 0)

    return nothing
end

"""
$(SIGNATURES)
"""
function polyatomic_moments_conserve!(
    w,
    h0::X,
    h1::X,
    h2::X,
    u::T,
    v::T,
    ω::T,
    ::Type{VDF{3,2}},
) where {X<:AA,T<:AA}
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

    return nothing
end

"""
$(SIGNATURES)
"""
function polyatomic_moments_conserve!(
    cons,
    h::X,
    r::X,
    u::T,
    v::T,
    w::T,
    ω::T,
    ::Type{VDF{2,3}},
) where {X<:AA,T<:AA}
    cons[1] = discrete_moments(h, u, ω, 0)
    cons[2] = discrete_moments(h, u, ω, 1)
    cons[3] = discrete_moments(h, v, ω, 1)
    cons[4] = discrete_moments(h, w, ω, 1)
    cons[5] =
        0.5 * (
            discrete_moments(h, u, ω, 2) +
            discrete_moments(h, v, ω, 2) +
            discrete_moments(h, w, ω, 2) +
            discrete_moments(r, u, ω, 0)
        )
    cons[6] = 0.5 * discrete_moments(r, u, ω, 0)

    return nothing
end


"""
$(SIGNATURES)
"""
function polyatomic_moments_conserve(h::X, b::X, r::X, u::T, ω::T) where {X<:AV,T<:AV}
    w = similar(h, 4)
    polyatomic_moments_conserve!(w, h, b, r, u, ω)
    return w
end

"""
$(SIGNATURES)
"""
function polyatomic_moments_conserve(
    h0::X,
    h1::X,
    h2::X,
    u::T,
    v::T,
    ω::T,
    ::Type{VDF{3,2}},
) where {X<:AA,T<:AA}
    w = similar(h0, 5)
    polyatomic_moments_conserve!(w, h0, h1, h2, u, v, ω, VDF{3,2})
    return w
end

"""
$(SIGNATURES)
"""
function polyatomic_moments_conserve(
    h::X,
    r::X,
    u::T,
    v::T,
    w::T,
    ω::T,
    ::Type{VDF{2,3}},
) where {X<:AA,T<:AA}
    cons = similar(h, 6)
    polyatomic_moments_conserve!(cons, h, r, u, v, w, ω, VDF{2,3})
    return cons
end

# ------------------------------------------------------------
# Multi-component gas
# ------------------------------------------------------------

"""
$(SIGNATURES)

Calculate moments of Gaussian distribution in multi-component gas
"""
function mixture_gauss_moments(prim::AM, inK)
    if eltype(prim) <: Integer
        Mu = OffsetArray(similar(prim, Float64, 7, axes(prim, 2)), 0:6, axes(prim, 2))
    else
        Mu = OffsetArray(similar(prim, 7, axes(prim, 2)), 0:6, axes(prim, 2))
    end

    MuL = similar(Mu)
    MuR = similar(Mu)

    if size(prim, 1) == 3

        Mxi = similar(Mu, 0:2, axes(prim, 2))
        for j in axes(prim, 2)
            _tu, _txi, _tuL, _tuR = gauss_moments(prim[:, j], inK)

            Mu[:, j] .= _tu
            Mxi[:, j] .= _txi
            MuL[:, j] .= _tuL
            MuR[:, j] .= _tuR
        end

        return Mu, Mxi, MuL, MuR

    elseif size(prim, 1) == 4

        Mv = similar(Mu)
        Mxi = similar(Mu, 0:2, axes(prim, 2))
        for j in axes(prim, 2)
            _tu, _tv, _txi, _tuL, _tuR = gauss_moments(prim[:, j], inK)

            Mu[:, j] .= _tu
            Mv[:, j] .= _tv
            Mxi[:, j] .= _txi
            MuL[:, j] .= _tuL
            MuR[:, j] .= _tuR
        end

        return Mu, Mv, Mxi, MuL, MuR

    elseif size(prim, 1) == 5

        Mv = similar(Mu)
        Mw = similar(Mu)

        for j in axes(prim, 2)
            _tu, _tv, _tw, _tuL, _tuR = gauss_moments(prim[:, j], inK)

            Mu[:, j] .= _tu
            Mv[:, j] .= _tv
            Mw[:, j] .= _tw
            MuL[:, j] .= _tuL
            MuR[:, j] .= _tuR
        end

        return Mu, Mv, Mw, MuL, MuR

    end
end


"""
$(SIGNATURES)

Calculate conservative moments of particle distribution in multi-component gas
"""
function mixture_moments_conserve(
    Mu::T,
    Mxi::T,
    alpha::I,
    delta::I,
) where {T<:OffsetMatrix,I<:Integer}

    Muv = similar(Mu, 3, size(Mu, 2))
    for j in axes(Muv, 2)
        Muv[:, j] .= moments_conserve(Mu[:, j], Mxi[:, j], alpha, delta)
    end

    return Muv

end

"""
$(SIGNATURES)
"""
function mixture_moments_conserve(
    Mu::T,
    Mv::T,
    Mw::T,
    alpha::I,
    beta::I,
    delta::I,
) where {T<:OffsetMatrix,I<:Integer}

    Muv = ifelse(size(Mw, 1) == 3, similar(Mu, 4, size(Mu, 2)), similar(Mu, 5, size(Mu, 2)))
    for j in axes(Muv, 2)
        Muv[:, j] .= moments_conserve(Mu[:, j], Mv[:, j], Mw[:, j], alpha, beta, delta)
    end

    return Muv

end

"""
$(SIGNATURES)

Discrete moments of conservative variables

1F1V
"""
function mixture_moments_conserve(f::X, u::T, ω::T) where {X<:AM,T<:AM}
    w = similar(f, 3, size(f, 2))
    for j in axes(w, 2)
        w[:, j] .= moments_conserve(f[:, j], u, ω)
    end

    return w
end

"""
$(SIGNATURES)

2F1V
"""
function mixture_moments_conserve(h::X, b::X, u::T, ω::T) where {X<:AM,T<:AM}

    w = similar(h, 3, size(h, 2))
    for j in axes(w, 2)
        w[:, j] .= moments_conserve(h[:, j], b[:, j], u, ω)
    end

    return w

end

"""
$(SIGNATURES)

4F1V
"""
function mixture_moments_conserve(
    h0::X,
    h1::X,
    h2::X,
    h3::X,
    u::T,
    ω::T,
) where {X<:AM,T<:AM}

    moments = similar(h0, 5, size(h0, 2))
    for j in axes(moments, 2)
        moments[:, j] .=
            moments_conserve(h0[:, j], h1[:, j], h2[:, j], h3[:, j], u[:, j], ω[:, j])
    end

    return moments

end

"""
$(SIGNATURES)

1F2V
"""
function mixture_moments_conserve(
    f::AA{T1,3},
    u::T,
    v::T,
    ω::T,
) where {T<:AA{T2,3}} where {T1,T2}

    w = similar(f, 4, size(f, 3))
    for j in axes(w, 2)
        w[:, j] .= moments_conserve(f[:, :, j], u[:, :, j], v[:, :, j], ω[:, :, j])
    end

    return w

end

"""
$(SIGNATURES)

2F2V
"""
function mixture_moments_conserve(
    h::X,
    b::X,
    u::T,
    v::T,
    ω::T,
) where {X<:AA{T1,3},T<:AA{T2,3}} where {T1,T2}

    w = similar(h, 4, size(f, 3))
    for j in axes(w, 2)
        w[:, j] .=
            moments_conserve(h[:, :, j], b[:, :, j], u[:, :, j], v[:, :, j], ω[:, :, j])
    end

    return w

end

"""
$(SIGNATURES)

3F2V
"""
function mixture_moments_conserve(
    h0::X,
    h1::X,
    h2::X,
    u::T,
    v::T,
    ω::T,
) where {X<:AA{T1,3},T<:AA{T2,3}} where {T1,T2}

    w = similar(h0, 5, size(h0, 3))
    for j in axes(w, 2)
        w[:, j] .= moments_conserve(
            h0[:, :, j],
            h1[:, :, j],
            h2[:, :, j],
            u[:, :, j],
            v[:, :, j],
            ω[:, :, j],
        )
    end

    return w

end

"""
$(SIGNATURES)

1F3V
"""
function mixture_moments_conserve(
    f::X,
    u::T,
    v::T,
    w::T,
    ω::T,
) where {X<:AA{T1,4},T<:AA{T2,4}} where {T1,T2}

    moments = similar(f, 5, size(f, 4))
    for j in axes(w, 2)
        moments[:, j] .= moments_conserve(
            f[:, :, :, j],
            u[:, :, :, j],
            v[:, :, :, j],
            w[:, :, :, j],
            ω[:, :, :, j],
        )
    end

    return moments

end


"""
$(SIGNATURES)

Calculate slope-related conservative moments under the assumption
`a = a1 + u * a2 + 0.5 * u^2 * a3`
"""
function mixture_moments_conserve_slope(
    a::AM,
    Mu::Y,
    Mxi::Y,
    alpha::Integer,
) where {Y<:OffsetMatrix}

    au = similar(a, 3, axes(a, 2))
    for j in axes(au, 2)
        au[:, j] .= moments_conserve_slope(a[:, j], Mu[:, j], Mxi[:, j], alpha)
    end

    return au

end

"""
$(SIGNATURES)
"""
function mixture_moments_conserve_slope(
    a::AM,
    Mu::Y,
    Mv::Y,
    Mxi::Y,
    alpha::I,
    beta::I,
) where {Y<:OffsetMatrix,I<:Integer}

    au = similar(a, 4, axes(a, 2))
    for j in axes(au, 2)
        au[:, j] .=
            moments_conserve_slope(a[:, j], Mu[:, j], Mv[:, j], Mxi[:, j], alpha, beta)
    end

    return au

end

"""
$(SIGNATURES)
"""
function mixture_moments_conserve_slope(
    a::X,
    Mu::Y,
    Mv::Y,
    Mw::Y,
    alpha::I,
    beta::I,
    delta::I,
) where {X<:AM,Y<:OffsetMatrix,I<:Integer}

    au = similar(a, 5, axes(a, 2))
    for j in axes(au, 2)
        au[:, j] .= moments_conserve_slope(
            a[:, j],
            Mu[:, j],
            Mv[:, j],
            Mw[:, j],
            alpha,
            beta,
            delta,
        )
    end

    return au

end


"""
$(SIGNATURES)

Compute conservative moments of polyatomic gas mixture
"""
function mixture_polyatomic_moments_conserve(h, b, r, u, ω)
    w = similar(h, 4, size(h)[end])
    for i in axes(w, 2)
        _w = @view w[:, i]
        _h = extract_last(h, i)
        _b = extract_last(b, i)
        _r = extract_last(r, i)
        _u = extract_last(u, i)
        _ω = extract_last(ω, i)

        polyatomic_moments_conserve!(_w, _h, _b, _r, _u, _ω)
    end

    return w
end

"""
$(SIGNATURES)
"""
function mixture_polyatomic_moments_conserve(h, b, r, u, v, ω, ::Type{VDF{3,2}})
    w = similar(h, 5, size(h)[end])
    for i in axes(w, 2)
        _w = @view w[:, i]
        _h = extract_last(h, i)
        _b = extract_last(b, i)
        _r = extract_last(r, i)
        _u = extract_last(u, i)
        _v = extract_last(v, i)
        _ω = extract_last(ω, i)

        polyatomic_moments_conserve!(_w, _h, _b, _r, _u, _v, _ω, VDF{3,2})
    end

    return w
end

"""
$(SIGNATURES)
"""
function mixture_polyatomic_moments_conserve(h, r, u, v, w, ω, ::Type{VDF{2,3}})
    cons = similar(h, 6, size(h)[end])
    for i in axes(cons, 2)
        _cons = @view cons[:, i]
        _h = extract_last(h, i)
        _r = extract_last(r, i)
        _u = extract_last(u, i)
        _v = extract_last(v, i)
        _w = extract_last(w, i)
        _ω = extract_last(ω, i)

        polyatomic_moments_conserve!(_cons, _h, _r, _u, _v, _w, _ω, VDF{2,3})
    end

    return cons
end
