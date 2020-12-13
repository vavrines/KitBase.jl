using OffsetArrays
import KitBase

cd(@__DIR__)
D = KitBase.read_dict("config.txt")
for key in keys(D)
    s = Symbol(key)
    @eval $s = $(D[key])
end

begin
    γ = KitBase.heat_capacity_ratio(inK, 1)
    set = KitBase.Setup(case, space, flux, collision, nSpecies, interpOrder, limiter, cfl, maxTime)
    pSpace = KitBase.PSpace1D(x0, x1, nx, pMeshType, nxg)
    vSpace = KitBase.VSpace1D(umin, umax, nu, vMeshType, nug)
    μᵣ = KitBase.ref_vhs_vis(knudsen, alphaRef, omegaRef)
    gas = KitBase.Particle(
        knudsen,
        mach,
        prandtl,
        inK,
        γ,
        omega,
        alphaRef,
        omegaRef,
        μᵣ,
        mass,
        0,
    )
    wL, primL, hL, bL, bcL, wR, primR, hR, bR, bcR = KitBase.ib_rh(gas.Ma, gas.γ, gas.K, vSpace.u)
    ib = KitBase.IB(wL, primL, bcL, wR, primR, bcR)
    ks = KitBase.SolverSet(set, pSpace, vSpace, gas, ib, pwd())
end

begin
    ctr = OffsetArray{KitBase.ControlVolumeParticle1D}(undef, eachindex(ks.pSpace.x))
    face = Array{KitBase.Interface1D}(undef, ks.pSpace.nx + 1)

    for i in eachindex(ctr)
        if i <= ks.pSpace.nx ÷ 2
            ctr[i] = KitBase.ControlVolumeParticle1D(
                ks.pSpace.x[i],
                ks.pSpace.dx[i],
                ks.ib.wL,
                ks.ib.primL,
            )
        else
            ctr[i] = KitBase.ControlVolumeParticle1D(
                ks.pSpace.x[i],
                ks.pSpace.dx[i],
                ks.ib.wR,
                ks.ib.primR,
            )
        end
    end

    for i = 1:ks.pSpace.nx+1
        face[i] = KitBase.Interface1D(ks.ib.wL)
    end

    ptc = KitBase.init_ptc(ks, ctr)
    ptc_new = deepcopy(ptc)
    ks.gas.np = length(ptc) ÷ 2
end

t = 0.0
dt = KitBase.timestep(ks, ctr, t)

KitBase.reconstruct!(ks, ctr)

@inbounds Threads.@threads for i = 1:ks.pSpace.nx+1
    KitBase.flux_equilibrium!(
        face[i].fw,
        ctr[i-1].w .+ 0.5 .* ctr[i-1].dx .* ctr[i-1].sw,
        ctr[i].w .- 0.5 .* ctr[i].dx .* ctr[i].sw,
        ks.gas.K,
        ks.gas.γ,
        ks.gas.μᵣ,
        ks.gas.ω,
        ks.gas.Pr,
        dt,
        0.5 * ctr[i-1].dx,
        0.5 * ctr[i].dx,
        ctr[i-1].sw,
        ctr[i].sw,
    )
end

KitBase.update_transport!(ks, ctr, ptc, ptc_new, dt)
KitBase.update_collision!(ks, ctr, ptc_new, face, dt, :bgk)
KitBase.update_boundary!(ks, ctr, ptc, ptc_new, face, dt, :bgk, :fix)

KitBase.plot_line(ks, ctr)