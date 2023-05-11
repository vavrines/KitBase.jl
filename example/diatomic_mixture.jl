using KitBase, Plots
using KitBase.ProgressMeter: @showprogress

cf = config_ntuple(
    space = "1d3f1v2s",
    cfl = 0.5,
    maxTime = 5.0,
    umin = -5.0,
    umax = 5.0, 
    nu = 48,
    knudsen = 0.01,
    prandtl = 0.72,
    inK = 2.0,
    inKr = 2.0,
    Tr0 = 91.5 / 273,
    Z0 = 18.1,
)

γ = heat_capacity_ratio(cf.inK, cf.inKr, 1)

set = set_setup(; cf...)
gas = PolyatomicMixture(
    Kn = cf.knudsen,
    Pr = cf.prandtl,
    K = cf.inK,
    Kr = cf.inKr,
)

vs = MVSpace1D(cf.umin, cf.umax, -8, 8, cf.nu)

_h0 = 
    0.5 * (1 / π)^1.5 .*
    (exp.(-(vs.u[:, 1] .- 1) .^ 2) .+ exp.(-(vs.u[:, 1] .+ 1) .^ 2))
_b0 = @. _h0 * gas.K / 2.0 / 1.0
_r0 = @. _h0 * gas.Kr / 2.0 / 5.0

_h1 = 
    0.5 * (1 / π)^1.5 .*
    (exp.(-(vs.u[:, 2] .- 1.5) .^ 2) .+ exp.(-(vs.u[:, 2] .+ 1.5) .^ 2)) .* gas.m[2]
_b1 = @. _h1 * gas.K / 2.0 / 1.0
_r1 = @. _h1 * gas.Kr / 2.0 / 3.0

h0 = hcat(_h0, _h1)
b0 = hcat(_b0, _b1)
r0 = hcat(_r0, _r1)

w0 = mixture_polyatomic_moments_conserve(h0, b0, r0, vs.u, vs.weights)
prim0 = mixture_conserve_prim(w0, gas.K, gas.Kr)

MH = zero(h0)
MB = zero(b0)
MR = zero(r0)
mixture_polyatomic_maxwellian!(MH, MB, MR, vs.u, prim0, gas.K, gas.Kr)

w1 = mixture_polyatomic_moments_conserve(MH, MB, MR, vs.u, vs.weights)
"""w1 == w0"""










plot(vs.u[:, 1], h0[:, 1])
plot!(vs.u[:, 2], h0[:, 2])
