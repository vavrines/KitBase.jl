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

    p = (x0=ps.x0, x1=ps.x1, u=vs.u, γ=gas.γ, K=gas.K)

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

        p = (p..., wL=wL, wR=wR, primL=primL, primR=primR)

        fw = function (x, p)
            if x <= (p.x0 + p.x1) / 2
                return p.wL
            else
                return p.wR
            end
        end

        bc = function (x, p)
            if x <= (p.x0 + p.x1) / 2
                return p.primL
            else
                return p.primR
            end
        end

        if set.space[1:4] == "1d0f"
            return fw, bc, p
        elseif set.space == "1d1f1v"
            ff = function (x, p)
                w = ifelse(x <= (p.x0 + p.x1) / 2, p.wL, p.wR)
                prim = conserve_prim(w, p.γ)
                h = maxwellian(p.u, prim)
                return h
            end

            return fw, ff, bc, p
        elseif set.space == "1d2f1v"
            ff = function (x, p)
                w = ifelse(x <= (p.x0 + p.x1) / 2, p.wL, p.wR)
                prim = conserve_prim(w, p.γ)
                h = maxwellian(p.u, prim)
                b = @. h * p.K / 2 / prim[end]
                return h, b
            end

            return fw, ff, bc, p
        elseif set.space == "1d1f3v"
            p = (p..., v=vs.v, w=vs.w)
            ff = function (x, p)
                w = ifelse(x <= (p.x0 + p.x1) / 2, p.wL, p.wR)
                prim = conserve_prim(w, p.γ)
                h = maxwellian(p.u, p.v, p.w, prim)
                return h
            end

            return fw, ff, bc, p
        end

    elseif set.nSpecies == 2
        mi, me = gas.mi, gas.me
        ni, ne = gas.ni, gas.ne

        primL = zeros(3, 2)
        primL[:, 1] .= [mi, MaL * sqrt(gam / 2.0 / (mi * ni + me * ne)), mi / 1.0]
        primL[:, 2] .= [me, MaL * sqrt(gam / 2.0 / (mi * ni + me * ne)), me / 1.0]

        primR = zeros(3, 2)
        for j in 1:2
            primR[1, j] = primL[1, j] * (gam + 1.0) * MaL^2 / ((gam - 1.0) * MaL^2 + 2.0)
            primR[2, j] = MaR * sqrt(gam / 2.0 / (mi * ni + me * ne)) * sqrt(ratioT)
            primR[3, j] = primL[3, j] / ratioT
        end

        wL = mixture_prim_conserve(primL, gam)
        wR = mixture_prim_conserve(primR, gam)

        p = (p..., wL=wL, wR=wR, primL=primL, primR=primR)

        fw = function (x, p)
            if x <= (p.x0 + p.x1) / 2
                return p.wL
            else
                return p.wR
            end
        end

        bc = function (x, p)
            if x <= (p.x0 + p.x1) / 2
                return p.primL
            else
                return p.primR
            end
        end

        if set.space == "1d2f1v"
            hL = mixture_maxwellian(vs.u, primL)
            hR = mixture_maxwellian(vs.u, primR)
            bL = similar(hL)
            bR = similar(hR)
            for j in 1:2
                bL[:, j] .= hL[:, j] .* gas.K ./ (2.0 .* primL[end, j])
                bR[:, j] .= hR[:, j] .* gas.K ./ (2.0 .* primR[end, j])
            end

            p = (p..., hL=hL, bL=bL, hR=hR, bR=bR, K=gas.K)

            ff = function (x, p)
                if x <= (ps.x0 + ps.x1) / 2
                    return p.hL, p.bL
                else
                    return p.hR, p.bR
                end
            end

            return fw, ff, bc, p
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
