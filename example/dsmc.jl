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
        matter,
        case,
        space,
        flux,
        collision,
        nSpecies,
        interpOrder,
        limiter,
        boundary,
        cfl,
        maxTime,
    )
    pSpace = KitBase.PSpace1D(x0, x1, nx, nxg)
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
    ctr = OffsetArray{KitBase.ControlVolumeParticle1D}(undef, eachindex(ks.ps.x))
    face = Array{KitBase.Interface1D}(undef, ks.ps.nx + 1)
    for i in eachindex(ctr)
        prim = [1.0, 0.0, primL[3] + 2.0 * (ks.ps.x[i] - ks.ps.x0), 1.0]

        ctr[i] = KitBase.ControlVolumeParticle1D(
            ks.ps.x[i],
            ks.ps.dx[i],
            KitBase.prim_conserve(prim, ks.gas.γ),
            prim,
            KitBase.vhs_collision_time(prim, ks.gas.μᵣ, ks.gas.ω),
        )
    end

    for i in 1:ks.ps.nx+1
        face[i] = KitBase.Interface1D(ks.ib.wL)
    end

    t = 0.0
    dt = KitBase.timestep(ks, ctr, t)
    res = zeros(4)
end

ptc = KitBase.init_ptc!(ks, ctr)

@showprogress for iter in 1:100
    KitBase.free_transport!(ks, ptc.x, ptc.v, ptc.flag, dt)
    KitBase.boundary!(ks, ctr, ptc, face, dt, :maxwell)
    KitBase.sort!(ks, ctr, ptc.x, ptc.idx, ptc.ref)
    KitBase.dsmc!(ks, ctr, ptc.ref, ptc.v, dt)
    KitBase.stat!(ks, ctr, ptc)
end

begin
    using Plots
    sol = zeros(ks.ps.nx, 4)
    for i in 1:ks.ps.nx
        sol[i, :] .= ctr[i].prim
    end

    plot(ks.ps.x[1:ks.ps.nx], sol[:, :])
end
