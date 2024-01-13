using KitBase, Plots
using KitBase.ProgressMeter: @showprogress

function step(w, prim, prim1, h, b, r, dt)
    τ = vhs_collision_time(prim[1:end-2], gas.μᵣ, gas.ω)
    Zr = rykov_zr(1.0 / prim[4], gas.T₀, gas.Z₀)
    τr = Zr * τ

    θ = 1 / prim1[end] + dt / τr * (1 / prim1[3] - 1 / prim1[end])
    prim1[end] = 1 / θ
    prim1[4] = 1 / ((2.5 / prim1[3] - θ) * 2 / 3)

    polyatomic_maxwellian!(MH, MB, MR, vs.u, prim1, gas.K, gas.Kr)

    for i in eachindex(h)
        #h[i] = (h[i] + dt / τ * MH[i]) / (1.0 + dt / τ)
        #b[i] = (b[i] + dt / τ * MB[i]) / (1.0 + dt / τ)
        #r[i] = (r[i] + dt / τ * MR[i]) / (1.0 + dt / τ)

        h[i] += dt / τ * (MH[i] - h[i])
        b[i] += dt / τ * (MB[i] - b[i])
        r[i] += dt / τ * (MR[i] - r[i])

        #h[i] = MH[i]
        #b[i] = MB[i]
        #r[i] = MR[i]
    end

    w[end] = polyatomic_moments_conserve(h, b, r, vs.u, vs.weights)[end]
    prim .= conserve_prim(w, gas.K, gas.Kr)

    return nothing
end

function solve(w, prim, prim1, h, b, r, dt, nt)
    whis = zeros(4, nt)
    primhis = zeros(5, nt)
    primhis1 = zeros(5, nt)
    hhis = zeros(vs.nu, nt)
    bhis = zeros(vs.nu, nt)
    rhis = zeros(vs.nu, nt)

    @showprogress for iter = 1:nt
        step(w, prim, prim1, h, b, r, dt)
        whis[:, iter] .= w
        primhis[:, iter] .= prim
        primhis1[:, iter] .= prim1
        hhis[:, iter] .= h
        bhis[:, iter] .= b
        rhis[:, iter] .= r
    end

    return whis, primhis, primhis1, hhis, bhis, rhis
end

cf = config_ntuple(
    maxTime = 1.0,
    umin = -5.0,
    umax = 5.0,
    nu = 80,
    knudsen = 0.01,
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

dt = 0.001
nt = cf.maxTime / dt |> Int
tran = collect(0:dt:cf.maxTime-dt)
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
#h0 = 
#    0.5 * (1 / π)^1.5 .*
#    (exp.(-(vs.u .- 1) .^ 2) .+ exp.(-(vs.u .+ 1) .^ 2))
#b0 = @. h0 * gas.K / 2.0 / 1.0
#r0 = @. h0 * gas.Kr / 2.0 / 5.0

prim0 = [
    0.31830988383259134,
    1.3631773076395243e-17,
    0.9259259887805802,
    0.6000000439882016,
    5.0,
]
w0 = prim_conserve(prim0, gas.γ, gas.Kr)

h0 = maxwellian(vs.u, [prim0[1:2]; prim0[4]])
b0 = energy_maxwellian(h0, prim0[4], gas.K)
r0 = energy_maxwellian(h0, prim0[5], gas.Kr)

begin
    MH = zero(h0)
    MB = zero(b0)
    MR = zero(r0)

    w = deepcopy(w0)
    prim = deepcopy(prim0)
    prim1 = deepcopy(prim0)
    h = deepcopy(h0)
    b = deepcopy(b0)
    r = deepcopy(r0)
end

wh, ph, ph1, hh, bh, rh = solve(w, prim, prim1, h, b, r, dt, nt)

plot(tran, 1 ./ prim_his[4, :]; label = "translation", xlabel = "t", ylabel = "T")
plot!(tran, 1 ./ prim_his[5, :]; label = "internal")
plot!(tran, 1 ./ prim_his[3, :]; label = "total")
savefig("rx_t.pdf")

#plot!(1 ./ prim_his1[4, :])
#plot!(1 ./ prim_his1[5, :])

contourf(tran, vs.u, hh; xlabel = "t", ylabel = "u")
savefig("rx_h.pdf")

contourf(tran, vs.u, rh; xlabel = "t", ylabel = "u")
savefig("rx_r.pdf")
