# ------------------------------------------------------------
# Limiter functions
# ------------------------------------------------------------

"""
$(SIGNATURES)

Linear average
"""
linear(sL, sR) = 0.5 * (sL + sR)


"""
$(SIGNATURES)

van Leer limiter
"""
vanleer(sL, sR) =
    (fortsign(1.0, sL) + fortsign(1.0, sR)) * abs(sL) * abs(sR) /
    (abs(sL) + abs(sR) + 1.e-7)

"""
$(SIGNATURES)

Triangle case
"""
function vanleer(sL, s, sR, connect = 2)
    δ = begin
        if connect == 2
            (
                (fortsign(1.0, sL) + fortsign(1.0, s)) * abs(sL) * abs(s) /
                (abs(sL) + abs(s) + 1.e-7),
                (fortsign(1.0, s) + fortsign(1.0, sR)) * abs(s) * abs(sR) /
                (abs(s) + abs(sR) + 1.e-7),
            )
        elseif connect == 3
            (
                (fortsign(1.0, sL) + fortsign(1.0, s)) * abs(sL) * abs(s) /
                (abs(sL) + abs(s) + 1.e-7),
                (fortsign(1.0, s) + fortsign(1.0, sR)) * abs(s) * abs(sR) /
                (abs(s) + abs(sR) + 1.e-7),
                (fortsign(1.0, sL) + fortsign(1.0, sR)) * abs(sL) * abs(sR) /
                (abs(sL) + abs(sR) + 1.e-7),
            )
        end
    end

    id = findmin(abs.(δ))[2]

    return δ[id]
end


"""
$(SIGNATURES)

MinMod limiter
"""
minmod(sL, sR) = 0.5 * (fortsign(1.0, sL) + fortsign(1.0, sR)) * min(abs(sR), abs(sL))

"""
$(SIGNATURES)

Triangle case
"""
function minmod(sL, s, sR, connect = 2)
    δ = begin
        if connect == 2
            (
                0.5 * (fortsign(1.0, sL) + fortsign(1.0, s)) * min(abs(s), abs(sL)),
                0.5 * (fortsign(1.0, s) + fortsign(1.0, sR)) * min(abs(sR), abs(s)),
            )
        elseif connect == 3
            (
                0.5 * (fortsign(1.0, sL) + fortsign(1.0, s)) * min(abs(s), abs(sL)),
                0.5 * (fortsign(1.0, s) + fortsign(1.0, sR)) * min(abs(sR), abs(s)),
                0.5 * (fortsign(1.0, sL) + fortsign(1.0, sR)) * min(abs(sR), abs(sL)),
            )
        end
    end

    id = findmin(abs.(δ))[2]

    return δ[id]
end


"""
$(SIGNATURES)

SuperBee limiter
"""
function superbee(sL, sR)
    if sR >= 0.5 * sL && sR <= 2.0 * sL
        return 0.5 * (fortsign(1.0, sL) + fortsign(1.0, sR)) * max(abs(sL), abs(sR))
    elseif sR < 0.5 * sL && sR > 2.0 * sL
        return (fortsign(1.0, sL) + fortsign(1.0, sR)) * min(abs(sL), abs(sR))
    else
        return 0.0
    end
end


"""
$(SIGNATURES)

van Albaba limiter
"""
vanalbaba(sL, sR) = (sL^2 * sR + sL * sR^2) / (sL^2 + sR^2 + 1.e-7)


"""
$(SIGNATURES)

5th-order WENO-JS interpolation
"""
function weno5(wL2, wL1, wN, wR1, wR2)
    ϵ = 1e-6

    β0 = 13.0 / 12.0 * (wN - 2.0 * wR1 + wR2)^2 + 1.0 / 4.0 * (3.0 * wN - 4.0 * wR1 + wR2)^2
    β1 = 13.0 / 12.0 * (wL1 - 2.0 * wN + wR1)^2 + 1.0 / 4.0 * (wL1 - wR1)^2
    β2 = 13.0 / 12.0 * (wL2 - 2.0 * wL1 + wN)^2 + 1.0 / 4.0 * (wL2 - 4.0 * wL1 + 3.0 * wN)^2

    # right interface
    dr0 = 0.3
    dr1 = 0.6
    dr2 = 0.1

    αr0 = dr0 / (ϵ + β0)^2
    αr1 = dr1 / (ϵ + β1)^2
    αr2 = dr2 / (ϵ + β2)^2

    ωr0 = αr0 / (αr0 + αr1 + αr2)
    ωr1 = αr1 / (αr0 + αr1 + αr2)
    ωr2 = αr2 / (αr0 + αr1 + αr2)

    qr0 = 1.0 / 3.0 * wN + 5.0 / 6.0 * wR1 - 1.0 / 6.0 * wR2
    qr1 = -1.0 / 6.0 * wL1 + 5.0 / 6.0 * wN + 1.0 / 3.0 * wR1
    qr2 = 1.0 / 3.0 * wL2 - 7.0 / 6.0 * wL1 + 11.0 / 6.0 * wN

    wR = ωr0 * qr0 + ωr1 * qr1 + ωr2 * qr2

    # left interface
    dl0 = 0.1
    dl1 = 0.6
    dl2 = 0.3

    αl0 = dl0 / (ϵ + β0)^2
    αl1 = dl1 / (ϵ + β1)^2
    αl2 = dl2 / (ϵ + β2)^2

    ωl0 = αl0 / (αl0 + αl1 + αl2)
    ωl1 = αl1 / (αl0 + αl1 + αl2)
    ωl2 = αl2 / (αl0 + αl1 + αl2)

    ql0 = 11.0 / 6.0 * wN - 7.0 / 6.0 * wR1 + 1.0 / 3.0 * wR2
    ql1 = 1.0 / 3.0 * wL1 + 5.0 / 6.0 * wN - 1.0 / 6.0 * wR1
    ql2 = -1.0 / 6.0 * wL2 + 5.0 / 6.0 * wL1 + 1.0 / 3.0 * wN

    wL = αl0 * ql0 + αl1 * ql1 + αl2 * ql2

    return wL, wR
end
