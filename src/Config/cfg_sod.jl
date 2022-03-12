"""
$(SIGNATURES)

Initialize Sod shock tube
"""
function ib_sod(
    set::AbstractSetup,
    ps::AbstractPhysicalSpace,
    vs::Union{AbstractVelocitySpace,Nothing},
    gas::AbstractProperty,
)

    if set.nSpecies == 1

        primL = [1.0, 0.0, 0.5]
        primR = [0.125, 0.0, 0.625]

        primL = ifelse(set.space[5:6] == "3v", [primL[1:2]; zeros(2); primL[end]], primL)
        primR = ifelse(set.space[5:6] == "3v", [primR[1:2]; zeros(2); primR[end]], primR)

        wL = prim_conserve(primL, gas.γ)
        wR = prim_conserve(primR, gas.γ)

        p = (
            x0 = ps.x0,
            x1 = ps.x1,
            wL = wL,
            wR = wR,
            primL = primL,
            primR = primR,
            γ = gas.γ,
        )

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
            p = (p..., u = vs.u)
            ff = function (x, p)
                w = ifelse(x <= (p.x0 + p.x1) / 2, p.wL, p.wR)
                prim = conserve_prim(w, p.γ)
                h = maxwellian(p.u, prim)
                return h
            end

            return fw, ff, bc, p
        elseif set.space == "1d2f1v"
            p = (p..., u = vs.u, K = gas.K)
            ff = function (x, p)
                w = ifelse(x <= (p.x0 + p.x1) / 2, p.wL, p.wR)
                prim = conserve_prim(w, p.γ)
                h = maxwellian(p.u, prim)
                b = h * p.K / 2 / prim[end]
                return h, b
            end

            return fw, ff, bc, p
        elseif set.space == "1d1f3v"
            p = (p, u = vs.u, v = vs.v, w = vs.w)
            ff = function (x, p)
                w = ifelse(x <= (p.x0 + p.x1) / 2, p.wL, p.wR)
                prim = conserve_prim(w, p.γ)
                h = maxwellian(p.u, p.v, p.w, prim)
                return h
            end

            return fw, ff, bc, p
        end

    elseif set.nSpecies == 2

        primL = zeros(3, 2)
        primL[:, 1] .= [gas.mi, 0.0, gas.mi / 2]
        primL[:, 2] .= [gas.me, 0.0, gas.me / 2]

        primR = zeros(3, 2)
        primR[:, 1] .= [gas.mi * 0.125, 0.0, gas.mi / 1.6]
        primR[:, 2] .= [gas.me * 0.125, 0.0, gas.me / 1.6]

        wL = mixture_prim_conserve(primL, gas.γ)
        wR = mixture_prim_conserve(primR, gas.γ)

        p = (
            x0 = ps.x0,
            x1 = ps.x1,
            wL = wL,
            wR = wR,
            primL = primL,
            primR = primR,
            γ = gas.γ,
        )

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
            hL = mixture_maxwellian(vs.u, primL)
            hR = mixture_maxwellian(vs.u, primR)
            p = (p..., hL = hL, hR = hR)
            ff = function (x, p)
                if x <= (p.x0 + p.x1) / 2
                    return p.hL
                else
                    return p.hR
                end
            end

            return fw, ff, bc, p
        elseif set.space == "1d2f1v"
            hL = mixture_maxwellian(vs.u, primL)
            hR = mixture_maxwellian(vs.u, primR)
            bL = similar(hL)
            bR = similar(hR)
            for j = 1:2
                bL[:, j] .= hL[:, j] .* K ./ (2.0 .* primL[end, j])
                bR[:, j] .= hR[:, j] .* K ./ (2.0 .* primR[end, j])
            end

            p = (p, hL = hL, hR = hR, bL = bL, bR = bR)

            ff = function (x, p)
                if x <= (p.x0 + p.x1) / 2
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
