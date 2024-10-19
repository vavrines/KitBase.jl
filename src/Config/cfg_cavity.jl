"""
$(SIGNATURES)

Initialize lid-driven cavity
"""
function ib_cavity(
    set::AbstractSetup,
    ps::AbstractPhysicalSpace,
    vs::Union{AbstractVelocitySpace,Nothing},
    gas::AbstractProperty,
    Um=0.15,
    Vm=0.0,
    Tm=1.0,
)
    if set.nSpecies == 1
        prim = [1.0, 0.0, 0.0, 1.0]
        w = prim_conserve(prim, gas.γ)
        primU = [1.0, Um, Vm, Tm]

        p = (y1=ps.y1, w=w, prim=prim, primU=primU)

        fw = function (args...)
            p = args[end]
            return p.w
        end

        bc = function (x, y, args...)
            p = args[end]

            if y == p.y1
                return p.primU
            else
                return p.prim
            end
        end

        if set.space[1:4] == "2d0f"
            return fw, bc, p
        elseif set.space == "2d1f2v"
            h = maxwellian(vs.u, vs.v, prim)
            p = (p..., h=h)
            ff = function (args...)
                p = args[end]
                return p.h
            end

            return fw, ff, bc, p
        elseif set.space == "2d2f2v"
            h = maxwellian(vs.u, vs.v, prim)
            b = h .* gas.K / 2.0 / prim[end]
            p = (p..., h=h, b=b)
            ff = function (args...)
                p = args[end]
                return p.h, p.b
            end

            return fw, ff, bc, p
        end
    end

    return nothing
end
