using Revise, OffsetArrays, ProgressMeter, Plots
import KitBase

begin
    cd(@__DIR__)
    D = KitBase.read_dict("config.txt")
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
    wL, primL, hL, bL, bcL, wR, primR, hR, bR, bcR =
        KitBase.ib_rh(gas.Ma, gas.γ, gas.K, vSpace.u)
    ib = KitBase.IB(wL, primL, bcL, wR, primR, bcR)
    ks = KitBase.SolverSet(set, pSpace, vSpace, gas, ib, pwd())

    ctr = OffsetArray{KitBase.ControlVolumeParticle1D}(undef, eachindex(ks.ps.x))
    face = Array{KitBase.Interface1D}(undef, ks.ps.nx + 1)

    for i in eachindex(ctr)
        if i <= ks.ps.nx ÷ 2
            ctr[i] = KitBase.ControlVolumeParticle1D(
                ks.ps.x[i],
                ks.ps.dx[i],
                ks.ib.wL,
                ks.ib.primL,
                #KitBase.prim_conserve([1.0, 0.0, 0.5], ks.gas.γ),
                #[1.0, 0.0, 0.5],
                KitBase.vhs_collision_time(ks.ib.primL, ks.gas.μᵣ, ks.gas.ω),
            )
        else
            ctr[i] = KitBase.ControlVolumeParticle1D(
                ks.ps.x[i],
                ks.ps.dx[i],
                ks.ib.wR,
                ks.ib.primR,
                #ks.ib.wL,
                #ks.ib.primL,
                #KitBase.prim_conserve([0.8, 0.0, 0.8], ks.gas.γ),
                #[0.8, 0.0, 0.8],
                KitBase.vhs_collision_time(ks.ib.primR, ks.gas.μᵣ, ks.gas.ω),
            )
        end
    end

    for i = 1:ks.ps.nx+1
        face[i] = KitBase.Interface1D(ks.ib.wL)
    end

    t = 0.0
    dt = KitBase.timestep(ks, ctr, t)
    res = zeros(3)
end

ptc = KitBase.init_ptc!(ks, ctr; mode = :soa, factor = 2)
ptc_new = deepcopy(ptc)

KitBase.bgk_transport!(ks, ctr, ptc, ptc_new, dt)
KitBase.bgk_collision!(ks, ctr, ptc_new, face, res)
KitBase.boundary!(ks, ctr, ptc_new, face, dt, :fix)
KitBase.duplicate!(ptc, ptc_new, ks.gas.np)



@showprogress for iter = 1:100#nt

    @inbounds Threads.@threads for i = 1:ks.ps.nx+1
        KitBase.flux_equilibrium!(
            face[i].fw,
            ctr[i-1].w,
            ctr[i].w,
            ks.gas.K,
            ks.gas.γ,
            ks.gas.μᵣ,
            ks.gas.ω,
            ks.gas.Pr,
            dt,
            0.5 * ctr[i-1].dx,
            0.5 * ctr[i].dx,
        )
    end

    #KitBase.update!(ks, ctr, ptc, ptc_new, face, dt, res; coll = :bgk, bc = :fix)

    KitBase.bgk_transport!(ks, ctr, ptc, ptc_new, dt)
    KitBase.bgk_collision!(ks, ctr, ptc_new, face, res)
    KitBase.boundary!(ks, ctr, ptc_new, face, dt, :fix)
    KitBase.duplicate!(ptc, ptc_new, ks.gas.np)

    t += dt
end

KitBase.plot_line(ks, ctr)
