using KitBase
import KitBase.BSON
cd(@__DIR__)

# save
set = Setup(space = "1d1f1v", boundary = "fix", maxTime = 0.1)
ps = PSpace1D(0, 1, 100, 1)
vs = VSpace1D(-5, 5, 28)
gas = Gas(Kn = 5e-4, K = 0.0, γ = 3.0)
ib = begin
    primL = [1.0, 0.0, 1 / 0.5]
    primR = [0.125, 0.0, 1 / 0.625]
    p = (γ = gas.γ, u = vs.u, primL = primL, primR = primR)
    _fw = function (x, p)
        prim = ifelse(x < 0.5, p.primL, p.primR)
        return prim_conserve(prim, p.γ)
    end
    _ff = function (x, p)
        prim = ifelse(x < 0.5, p.primL, p.primR)
        return maxwellian(p.u, prim)
    end
    _bc = function (x, p)
        prim = ifelse(x < 0.5, p.primL, p.primR)
        return prim
    end
    IB1F(_fw, _ff, _bc, p)
end
ks = SolverSet(set, ps, vs, gas, ib)
BSON.@save "ks.bson" ks

# load
using OffsetArrays
BSON.@load "ks.bson" ks

w0, f0, prim0 = ks.ib(1)
