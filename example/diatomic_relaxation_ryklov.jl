using KitBase, Plots
using KitBase.ProgressMeter: @showprogress

cf = config_ntuple(;
    matter="gas",
    case="relaxation",
    space="1d3f1v",
    collision="rykov",
    cfl=0.5,
    maxTime=5.0,
    umin=-5.0,
    umax=5.0,
    nu=32,
    knudsen=0.1,
    prandtl=0.72,
    inK=2.0,
    inKr=2.0,
    omega=0.81,
    Tr0=91.5 / 273,
    Z0=18.1,
    sigma=1 / 1.55,
    omega1=0.2354,
    omega2=0.3049,
)

γ = heat_capacity_ratio(cf.inK, cf.inKr, 1)

set = set_setup(; cf...)
ps = set_geometry(; cf...)
vs = set_velocity(; cf...)
gas =
    DiatomicGas(; Kn=cf.knudsen, Pr=cf.prandtl, K=cf.inK, Kr=cf.inKr, ω=cf.omega, T₀=cf.Tr0)

prim0 = [1.0, 0.0, 1.6556, 1.0, 100.0]
w0 = prim_conserve(prim0, γ, gas.Kr)

h0 = zeros(vs.nu)
b0 = similar(h0)
r0 = similar(h0)

Ht = similar(h0)
Bt = similar(b0)
Rt = similar(r0)
Hr = similar(h0)
Br = similar(b0)
Rr = similar(r0)

maxwellian!(Ht, Bt, Rt, Hr, Br, Rr, vs.u, prim0, gas.K, gas.Kr)

Zr = rykov_zr(1 / prim0[4], gas.T₀, gas.Z₀)

@. h0 = (1.0 - 1.0 / Zr) * Ht + 1.0 / Zr * Hr
@. b0 = (1.0 - 1.0 / Zr) * Bt + 1.0 / Zr * Br
@. r0 = (1.0 - 1.0 / Zr) * Rt + 1.0 / Zr * Rr

w0 .= diatomic_moments_conserve(h0, b0, r0, vs.u, vs.weights)
prim0 .= conserve_prim(w0, gas.K, gas.Kr)
ks = SolverSet(set, ps, vs, gas, nothing)

w = deepcopy(w0)
prim = deepcopy(prim0)
h = deepcopy(h0)
b = deepcopy(b0)
r = deepcopy(r0)

function step(KS, w, prim, h, b, r, dt)
    w_old = deepcopy(w)
    q = heat_flux(h, b, r, prim, KS.vs.u, KS.vs.weights)

    MHT = similar(h)
    MBT = similar(b)
    MRT = similar(r)
    MHR = similar(h)
    MBR = similar(b)
    MRR = similar(r)

    maxwellian!(MHT, MBT, MRT, MHR, MBR, MRR, KS.vs.u, prim, KS.gas.K, KS.gas.Kr)
    τ_old = vhs_collision_time(prim[1:end-1], KS.gas.μᵣ, KS.gas.ω)
    Zr = rykov_zr(1.0 / prim[4], KS.gas.T₀, KS.gas.Z₀)
    Er0_old = 0.5 * sum(@. KS.vs.weights * ((1.0 / Zr) * MRR + (1.0 - 1.0 / Zr) * MRT))

    w[4] += dt * (Er0_old - w_old[4]) / τ_old
    prim .= conserve_prim(w, KS.gas.K, KS.gas.Kr)

    #--- calculate M^{n+1} and tau^{n+1} ---#
    maxwellian!(MHT, MBT, MRT, MHR, MBR, MRR, KS.vs.u, prim, KS.gas.K, KS.gas.Kr)

    SHT = similar(h)
    SBT = similar(b)
    SRT = similar(r)
    SHR = similar(h)
    SBR = similar(b)
    SRR = similar(r)
    rykov!(
        SHT,
        SBT,
        SRT,
        SHR,
        SBR,
        SRR,
        KS.vs.u,
        MHT,
        MBT,
        MRT,
        MHR,
        MBR,
        MRR,
        q,
        prim,
        KS.gas.Pr,
        KS.gas.K,
        KS.gas.σ,
        KS.gas.ω₁,
        KS.gas.ω₂,
    )

    MH = (1.0 - 1.0 / Zr) * (MHT + SHT) + 1.0 / Zr * (MHR + SHR)
    MB = (1.0 - 1.0 / Zr) * (MBT + SBT) + 1.0 / Zr * (MBR + SBR)
    MR = (1.0 - 1.0 / Zr) * (MRT + SRT) + 1.0 / Zr * (MRR + SRR)

    τ = vhs_collision_time(prim[1:end-1], KS.gas.μᵣ, KS.gas.ω)

    #--- update distribution function ---#
    for i in eachindex(h)
        h[i] = (h[i] + dt / τ * MH[i]) / (1.0 + dt / τ)
        b[i] = (b[i] + dt / τ * MB[i]) / (1.0 + dt / τ)
        r[i] = (r[i] + dt / τ * MR[i]) / (1.0 + dt / τ)
    end
end

dt = 0.001
nt = 2000

function solve(KS, w, prim, h, b, r, dt, nt)
    whis = zeros(4, nt)

    @showprogress for iter in 1:nt
        step(KS, w, prim, h, b, r, dt)
        whis[:, iter] .= w
    end

    return whis
end

w_his = solve(ks, w, prim, h, b, r, dt, nt)
prim_his = zeros(5, nt)
for i in 1:nt
    prim_his[:, i] .= conserve_prim(w_his[:, i], ks.gas.K, ks.gas.Kr)
end

plot(1 ./ prim_his[4, :])
plot!(1 ./ prim_his[5, :])
