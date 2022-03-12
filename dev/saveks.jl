using KitBase
import KitBase.BSON

cd(@__DIR__)

# save
set = Setup(space = "1d1f1v", boundary = "fix", maxTime = 0.1)
ps = PSpace1D(0, 1, 100, 1)
vs = VSpace1D(-5, 5, 28)
gas = Gas(Kn = 5e-4, K = 0.0, γ = 3.0)
ib = begin
    _fw = function (x)
        ρ, U, T = begin
            if x < 0.5
                1.0, 0.0, 1 / 0.5
            else
                0.125, 0.0, 1 / 0.625
            end
        end
        prim = [ρ, U, 1 / T]
        return prim_conserve(prim, gas.γ)
    end
    _ff = function (x)
        w = _fw(x)
        prim = conserve_prim(w, gas.γ)
        return maxwellian(vs.u, prim)
    end
    _bc = function (x)
        w = _fw(x)
        return conserve_prim(w, gas.γ)
    end
    IB1F(_fw, _ff, _bc)
end
ks = SolverSet(set, ps, vs, gas, ib)

BSON.@save "ks.bson" ks

using OffsetArrays

# load
BSON.@load "ks.bson" ks
gas = ks.gas
ks.ib.fw(1)