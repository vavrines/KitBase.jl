using KitBase, Plots
using KitBase.ProgressMeter: @showprogress

cf = config_ntuple(
    matter = "gas",
    case = "relaxation",
    space = "1d3f1v",
    collision = "rykov",
    cfl = 0.5,
    maxTime = 5.0,
    umin = -5.0,
    umax = 5.0,
    nu = 80,
    knudsen = 0.1,
    prandtl = 0.72,
    inK = 2.0,
    inKr = 2.0,
    omega = 0.81,
    Tr0 = 91.5 / 273,
    Z0 = 18.1,
    sigma = 1 / 1.55,
    omega1 = 0.2354,
    omega2 = 0.3049,
)

γ = heat_capacity_ratio(cf.inK, cf.inKr, 1)

set = set_setup(; cf...)
ps = set_geometry(; cf...)
vs = set_velocity(; cf...)
gas = DiatomicGas(
    Kn = cf.knudsen,
    Pr = cf.prandtl,
    K = cf.inK,
    Kr = cf.inKr,
    ω = cf.omega,
    T₀ = cf.Tr0,
)

#--- 1. non-equilibrium distributions ---#
h0 = 
    0.5 * (1 / π)^1.5 .*
    (exp.(-(vs.u .- 1) .^ 2) .+ exp.(-(vs.u .+ 1) .^ 2))

b0 = @. h0 * gas.K / 2.0 / 1.0
r0 = @. h0 * gas.Kr / 2.0 / 5.0

plot(vs.u, h0)

w0 = polyatomic_moments_conserve(h0, b0, r0, vs.u, vs.weights)
prim0 = conserve_prim(w0, gas.K, gas.Kr)

MH = zero(h0)
MB = zero(b0)
MR = zero(r0)

polyatomic_maxwellian!(MH, MB, MR, vs.u, prim0, gas.K, gas.Kr)

w1 = polyatomic_moments_conserve(MH, MB, MR, vs.u, vs.weights)
prim1 = conserve_prim(w1, gas.K, gas.Kr)
"""prim0 == prim1"""

#--- relaxation ---#
Zr = rykov_zr(1 / prim0[4], gas.T₀, gas.Z₀)

w = deepcopy(w0)
prim = deepcopy(prim0)
h = deepcopy(h0)
b = deepcopy(b0)
r = deepcopy(r0)

function step(w, prim, h, b, r, dt)
    MH = similar(h)
    MB = similar(b)
    MR = similar(r)
    
    τ = vhs_collision_time(prim[1:end-2], gas.μᵣ, gas.ω)
    Zr = rykov_zr(1.0 / prim[4], gas.T₀, gas.Z₀)
    τr = Zr * τ

    Tr = 1 / prim[end] + dt * (1/prim[3] - 1/prim[end]) / τr
    prim[end] = 1/Tr

    w .= prim_conserve(prim, gas.γ, gas.Kr)
    
    polyatomic_maxwellian!(MH, MB, MR, vs.u, prim, gas.K, gas.Kr)

    for i in eachindex(h)
        h[i] = (h[i] + dt / τ * MH[i]) / (1.0 + dt / τ)
        b[i] = (b[i] + dt / τ * MB[i]) / (1.0 + dt / τ)
        r[i] = (r[i] + dt / τ * MR[i]) / (1.0 + dt / τ)
    end
end

dt = 0.001
nt = 2000

function solve(w, prim, h, b, r, dt, nt)
    whis = zeros(4, nt)

    @showprogress for iter = 1:nt
        step(w, prim, h, b, r, dt)
        whis[:, iter] .= w
    end

    return whis
end

w_his = solve(w, prim, h, b, r, dt, nt)
prim_his = zeros(5, nt)
for i = 1:nt
    prim_his[:, i] .= conserve_prim(w_his[:, i], gas.K, gas.Kr)
end

plot(1 ./ prim_his[4, :])
plot!(1 ./ prim_his[5, :])
