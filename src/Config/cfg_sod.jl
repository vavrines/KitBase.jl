"""
    1d0f0v: ib_sod(γ)
    1d1f1v: ib_sod(γ, u::T) where {T<:AbstractArray{<:AbstractFloat,1}}
    1d1f3v: ib_sod(γ, u::T, v::T, w::T) where {T<:AbstractArray{<:AbstractFloat,3}}
    1d2f1v: ib_sod(γ, u::T, K) where {T<:AbstractArray{<:AbstractFloat,1}}

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

        wL = prim_conserve(primL, gas.γ)
        wR = prim_conserve(primR, gas.γ)

        fw = function(x)
            if x <= (ps.x0 + ps.x1) / 2
                return wL
            else
                return wR
            end
        end

        bc = function(x)
            if x <= (ps.x0 + ps.x1) / 2
                return primL
            else
                return primR
            end
        end

        if set.space[1:4] == "1d0f"
            return fw, bc
        elseif set.space == "1d1f1v"
            ff = function(x)
                w = fw(x)
                prim = conserve_prim(w, gas.γ)
                h = maxwellian(vs.u, prim)
                return h
            end

            return fw, ff, bc
        elseif set.space == "1d2f1v"
            ff = function(x)
                w = fw(x)
                prim = conserve_prim(w, gas.γ)
                h = maxwellian(vs.u, prim)
                b = h * gas.K / 2 / prim[end]
                return h, b
            end

            return fw, ff, bc
        elseif set.space == "1d1f3v"
            ff = function(x)
                w = fw(x)
                prim = conserve_prim(w, gas.γ)
                h = maxwellian(vs.u, vs.v, vs.w, prim)
                return h
            end

            return fw, ff, bc
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

        fw = function(x)
            if x <= (ps.x0 + ps.x1) / 2
                return wL
            else
                return wR
            end
        end

        bc = function(x)
            if x <= (ps.x0 + ps.x1) / 2
                return primL
            else
                return primR
            end
        end

        if set.space[1:4] == "1d0f"
            return fw, bc
        elseif set.space == "1d1f1v"
            hL = mixture_maxwellian(vs.u, primL)
            hR = mixture_maxwellian(vs.u, primR)

            ff = function(x)
                if x <= (ps.x0 + ps.x1) / 2
                    return hL
                else
                    return hR
                end
            end

            return fw, ff, bc
        elseif set.space == "1d2f1v"
            hL = mixture_maxwellian(vs.u, primL)
            hR = mixture_maxwellian(vs.u, primR)
            bL = similar(hL)
            bR = similar(hR)
            for j = 1:2
                bL[:, j] .= hL[:, j] .* K ./ (2.0 .* primL[end, j])
                bR[:, j] .= hR[:, j] .* K ./ (2.0 .* primR[end, j])
            end

            ff = function(x)
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
