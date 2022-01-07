"""
$(SIGNATURES)

Initialize Rankine-Hugoniot relation
"""
function ib_rh(
    set::AbstractSetup,
    ps::AbstractPhysicalSpace,
    vs::Union{AbstractVelocitySpace,Nothing},
    gas::AbstractProperty,
)

    gam = gas.γ
    MaL = gas.Ma
    MaR = sqrt((MaL^2 * (gam - 1.0) + 2.0) / (2.0 * gam * MaL^2 - (gam - 1.0)))
    ratioT =
        (1.0 + (gam - 1.0) / 2.0 * MaL^2) * (2.0 * gam / (gam - 1.0) * MaL^2 - 1.0) /
        (MaL^2 * (2.0 * gam / (gam - 1.0) + (gam - 1.0) / 2.0))

    if set.nSpecies == 1

        primL = [1.0, MaL * sqrt(gam / 2.0), 1.0]
        primR = [
            primL[1] * (gam + 1.0) * MaL^2 / ((gam - 1.0) * MaL^2 + 2.0),
            MaR * sqrt(gam / 2.0) * sqrt(ratioT),
            primL[3] / ratioT,
        ]

        primL = ifelse(set.space[5:6] == "3v", [primL[1:2]; zeros(2); primL[end]], primL)
        primR = ifelse(set.space[5:6] == "3v", [primR[1:2]; zeros(2); primR[end]], primR)

        wL = prim_conserve(primL, gam)
        wR = prim_conserve(primR, gam)

        fw = function (x)
            if x <= (ps.x0 + ps.x1) / 2
                return wL
            else
                return wR
            end
        end

        bc = function (x)
            if x <= (ps.x0 + ps.x1) / 2
                return primL
            else
                return primR
            end
        end

        if set.space[1:4] == "1d0f"
            return fw, bc
        elseif set.space == "1d1f1v"
            ff = function (x)
                w = fw(x)
                prim = conserve_prim(w, gas.γ)
                h = maxwellian(vs.u, prim)
                return h
            end

            return fw, ff, bc
        elseif set.space == "1d2f1v"
            ff = function (x)
                w = fw(x)
                prim = conserve_prim(w, gas.γ)
                h = maxwellian(vs.u, prim)
                b = h * gas.K / 2 / prim[end]
                return h, b
            end

            return fw, ff, bc
        elseif set.space == "1d1f3v"
            ff = function (x)
                w = fw(x)
                prim = conserve_prim(w, gas.γ)
                h = maxwellian(vs.u, vs.v, vs.w, prim)
                return h
            end

            return fw, ff, bc
        end

    elseif set.nSpecies == 2

        primL = zeros(3, 2)
        primL[:, 1] .= [mi, MaL * sqrt(gam / 2.0), mi / 1.0]
        primL[:, 2] .= [me, MaL * sqrt(gam / 2.0), me / 1.0]

        primR = zeros(3, 2)
        for j = 1:2
            primR[1, j] = primL[1, j] * (gam + 1.0) * MaL^2 / ((gam - 1.0) * MaL^2 + 2.0)
            primR[2, j] = MaR * sqrt(gam / 2.0 / (mi * ni + me * ne)) * sqrt(ratioT)
            primR[3, j] = primL[3, j] / ratioT
        end

        wL = mixture_prim_conserve(primL, gam)
        wR = mixture_prim_conserve(primR, gam)

        fw = function (x)
            if x <= (ps.x0 + ps.x1) / 2
                return wL
            else
                return wR
            end
        end

        bc = function (x)
            if x <= (ps.x0 + ps.x1) / 2
                return primL
            else
                return primR
            end
        end

        if set.space == "1d2f1v"
            hL = mixture_maxwellian(vs.u, primL)
            hR = mixture_maxwellian(vs.u, primR)
            bL = similar(hL)
            bR = similar(hR)
            for j = 1:2
                bL[:, j] .= hL[:, j] .* gas.K ./ (2.0 .* primL[end, j])
                bR[:, j] .= hR[:, j] .* gas.K ./ (2.0 .* primR[end, j])
            end

            ff = function (x)
                if x <= (ps.x0 + ps.x1) / 2
                    return hL, bL
                else
                    return hR, bR
                end
            end

            return fw, ff, bc
        end

    end

    return nothing

end

"""
$(SIGNATURES)
"""
function ib_rh(Ma::Real, gam::Real)
    MaL = Ma
    MaR = sqrt((MaL^2 * (gam - 1.0) + 2.0) / (2.0 * gam * MaL^2 - (gam - 1.0)))
    ratioT =
        (1.0 + (gam - 1.0) / 2.0 * MaL^2) * (2.0 * gam / (gam - 1.0) * MaL^2 - 1.0) /
        (MaL^2 * (2.0 * gam / (gam - 1.0) + (gam - 1.0) / 2.0))

    primL = [1.0, MaL * sqrt(gam / 2.0), 1.0]
    primR = [
        primL[1] * (gam + 1.0) * MaL^2 / ((gam - 1.0) * MaL^2 + 2.0),
        MaR * sqrt(gam / 2.0) * sqrt(ratioT),
        primL[3] / ratioT,
    ]

    return primL, primR
end
