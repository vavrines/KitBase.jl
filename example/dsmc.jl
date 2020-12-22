using Revise, OffsetArrays, ProgressMeter
import KitBase

begin
    cd(@__DIR__)
    D = KitBase.read_dict("poiseuille.txt")
    for key in keys(D)
        s = Symbol(key)
        @eval $s = $(D[key])
    end

    γ = KitBase.heat_capacity_ratio(inK, 1)
    set = KitBase.Setup(
        case,
        space,
        flux,
        collision,
        nSpecies,
        interpOrder,
        limiter,
        cfl,
        maxTime,
    )
    pSpace = KitBase.PSpace1D(x0, x1, nx, pMeshType, nxg)
    vSpace = KitBase.VSpace1D(umin, umax, nu, vMeshType, nug)
    μᵣ = KitBase.ref_vhs_vis(knudsen, alphaRef, omegaRef)
    gas =
        KitBase.Gas(knudsen, mach, prandtl, inK, γ, omega, alphaRef, omegaRef, μᵣ, mass, 0)

    primL = [1.0, 0.0, -1.0, 1.0] # left wall
    primR = [1.0, 0.0, 1.0, 1.0] # right wall
    wL = KitBase.prim_conserve(primL, γ)
    wR = KitBase.prim_conserve(primR, γ)
    ib = KitBase.IB(wL, primL, primL, wR, primR, primR)

    ks = KitBase.SolverSet(set, pSpace, vSpace, gas, ib, pwd())
end

begin
    ctr = OffsetArray{KitBase.ControlVolumeParticle1D}(undef, eachindex(ks.pSpace.x))
    face = Array{KitBase.Interface1D}(undef, ks.pSpace.nx + 1)
    for i in eachindex(ctr)
        prim = [1.0, 0.0, primL[3] + 2.0 * (ks.pSpace.x[i] - ks.pSpace.x0), 1.0]

        ctr[i] = KitBase.ControlVolumeParticle1D(
            ks.pSpace.x[i],
            ks.pSpace.dx[i],
            KitBase.prim_conserve(prim, ks.gas.γ),
            prim,
            KitBase.vhs_collision_time(prim, ks.gas.μᵣ, ks.gas.ω),
        )
    end

    for i = 1:ks.pSpace.nx+1
        face[i] = KitBase.Interface1D(ks.ib.wL)
    end

    t = 0.0
    dt = KitBase.timestep(ks, ctr, t)
    res = zeros(4)
end

ptc = KitBase.init_ptc!(ks, ctr)

@showprogress for iter = 1:100
    KitBase.transport!(ks, ptc.x, ptc.v, ptc.flag, dt)
    KitBase.sort!(ks, ctr, ptc.x, ptc.idx, ptc.ref)
    KitBase.dsmc!(ks, ctr, ptc.ref, ptc.v, dt)
    KitBase.stat!(ks, ctr, ptc)
end

begin
    using Plots
    sol = zeros(ks.pSpace.nx, 4)
    for i = 1:ks.pSpace.nx
        sol[i, :] .= ctr[i].prim
    end

    plot(ks.pSpace.x[1:ks.pSpace.nx], sol[:, :])
end
