using ProgressMeter, Plots
import KitBase

begin
    # case
    case = "relaxation"
    space = "1d3f1v"
    flux = "kfvs"
    collision = "rykov"
    nSpecies = 1
    interpOrder = 1
    limiter = "vanleer"
    boundary = "fix"
    cfl = 0.8
    maxTime = 5.0

    # phase space
    x0 = -1.0
    x1 = 1.0
    nx = 1
    pMeshType = "uniform"
    nxg = 0

    # velocity space
    umin = -5.0
    umax = 5.0
    nu = 32
    vMeshType = "rectangle"
    nug = 0

    # gas
    knudsen = 0.1
    mach = 0.0
    prandtl = 0.72
    inK = 2.0
    inKr = 2.0
    omega = 0.81
    alphaRef = 1.0
    omegaRef = 0.5
    Tr0 = 91.5 / 273
    Z0 = 18.1
    sigma = 1 / 1.55
    omega1 = 0.2354
    omega2 = 0.3049
end

γ = KitBase.heat_capacity_ratio(inK, inKr, 1)
set = KitBase.Setup(case, space, flux, collision, nSpecies, interpOrder, limiter, boundary, cfl, maxTime)
pSpace = KitBase.PSpace1D(x0, x1, nx)
vSpace = KitBase.VSpace1D(umin, umax, nu)
gas = KitBase.DiatomicGas(
    knudsen,
    mach,
    prandtl,
    inK,
    inKr,
    γ,
    omega,
    alphaRef,
    omegaRef,
    Tr0,
    Z0,
    sigma,
    omega1,
    omega2,
)

prim0 = [1.0, 0.0, 1.6556, 1.0, 100.0]
w0 = KitBase.prim_conserve(prim0, γ, inKr)

h0 = zeros(vSpace.nu)
b0 = similar(h0)
r0 = similar(h0)

Ht = similar(h0)
Bt = similar(b0)
Rt = similar(r0)
Hr = similar(h0)
Br = similar(b0)
Rr = similar(r0)

KitBase.maxwellian!(Ht, Bt, Rt, Hr, Br, Rr, vSpace.u, prim0, inK, inKr)

Zr = KitBase.rykov_zr(1/prim0[4], Tr0, Z0)

@. h0 = (1.0 - 1.0 / Zr) * Ht + 1.0 / Zr * Hr
@. b0 = (1.0 - 1.0 / Zr) * Bt + 1.0 / Zr * Br
@. r0 = (1.0 - 1.0 / Zr) * Rt + 1.0 / Zr * Rr

w0 .= KitBase.diatomic_moments_conserve(h0, b0, r0, vSpace.u, vSpace.weights)
prim0 .= KitBase.conserve_prim(w0, inK, inKr)

ib = KitBase.IB3F(
    w0,
    prim0,
    h0,
    b0,
    r0,
    prim0,
    nothing,
    nothing,
    nothing,
    w0,
    prim0,
    h0,
    b0,
    r0,
    prim0,
    nothing,
    nothing,
    nothing,
)

cd(@__DIR__)
folder = pwd()

ks = KitBase.SolverSet(set, pSpace, vSpace, gas, ib, folder)

w = deepcopy(w0)
prim = deepcopy(prim0)
h = deepcopy(h0)
b = deepcopy(b0)
r = deepcopy(r0)

function step(KS, w, prim, h, b, r, dt)

    w_old = deepcopy(w)
    q = KitBase.heat_flux(h, b, r, prim, KS.vSpace.u, KS.vSpace.weights)

    MHT = similar(h)
    MBT = similar(b)
    MRT = similar(r)
    MHR = similar(h)
    MBR = similar(b)
    MRR = similar(r)
    
    KitBase.maxwellian!(MHT, MBT, MRT, MHR, MBR, MRR, KS.vSpace.u, prim, KS.gas.K, KS.gas.Kr)
    τ_old = KitBase.vhs_collision_time(prim[1:end-1], KS.gas.μᵣ, KS.gas.ω)
    Zr = KitBase.rykov_zr(1.0 / prim[4], KS.gas.T₀, KS.gas.Z₀)
    Er0_old = 0.5 * sum(@. KS.vSpace.weights * ((1.0 / Zr) * MRR + (1.0 - 1.0 / Zr) * MRT))
    
    w[4] += dt * (Er0_old - w_old[4]) / τ_old
    prim .= KitBase.conserve_prim(w, KS.gas.K, KS.gas.Kr)

    #--- calculate M^{n+1} and tau^{n+1} ---#
    KitBase.maxwellian!(MHT, MBT, MRT, MHR, MBR, MRR, KS.vSpace.u, prim, KS.gas.K, KS.gas.Kr)

    SHT = similar(h)
    SBT = similar(b)
    SRT = similar(r)
    SHR = similar(h)
    SBR = similar(b)
    SRR = similar(r)
    KitBase.rykov!(
        SHT,
        SBT,
        SRT,
        SHR,
        SBR,
        SRR,
        KS.vSpace.u,
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

    τ = KitBase.vhs_collision_time(prim[1:end-1], KS.gas.μᵣ, KS.gas.ω)

    #--- update distribution function ---#
    for i in eachindex(h)
        h[i] = (h[i] + dt / τ * MH[i]) / (1.0 + dt / τ)
        b[i] = (b[i] + dt / τ * MB[i]) / (1.0 + dt / τ)
        r[i] = (r[i] + dt / τ * MR[i]) / (1.0 + dt / τ)
    end

end

dt = 0.001
nt = 100

function solve(KS, w, prim, h, b, r, dt, nt)
    whis = zeros(4, nt)

    @showprogress for iter = 1:nt
        step(KS, w, prim, h, b, r, dt)
        whis[:, iter] .= w
    end

    return whis
end

w_his = solve(ks, w, prim, h, b, r, dt, nt)
prim_his = zeros(5, nt)
for i in 1:nt
    prim_his[:, i] .= KitBase.conserve_prim(w_his[:, i], ks.gas.K, ks.gas.Kr)
end

plot(1 ./ prim_his[4, :])
plot!(1 ./ prim_his[5, :])