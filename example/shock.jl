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
    gas = KitBase.Gas(
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
    wL, primL, hL, bL, bcL, wR, primR, hR, bR, bcR =
        KitBase.ib_rh(gas.Ma, gas.γ, gas.K, vSpace.u)
    ib = KitBase.IB(wL, primL, bcL, wR, primR, bcR)
    ks = KitBase.SolverSet(set, pSpace, vSpace, gas, ib, pwd())

    ctr = OffsetArray{KitBase.ControlVolumeParticle1D}(undef, eachindex(ks.pSpace.x))
    face = Array{KitBase.Interface1D}(undef, ks.pSpace.nx + 1)

    for i in eachindex(ctr)
        if i <= ks.pSpace.nx ÷ 2
            ctr[i] = KitBase.ControlVolumeParticle1D(
                ks.pSpace.x[i],
                ks.pSpace.dx[i],
                ks.ib.wL,
                ks.ib.primL,
                #KitBase.prim_conserve([1.0, 0.0, 0.5], ks.gas.γ),
                #[1.0, 0.0, 0.5],
                KitBase.vhs_collision_time(ks.ib.primL, ks.gas.μᵣ, ks.gas.ω),
            )
        else
            ctr[i] = KitBase.ControlVolumeParticle1D(
                ks.pSpace.x[i],
                ks.pSpace.dx[i],
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

    for i = 1:ks.pSpace.nx+1
        face[i] = KitBase.Interface1D(ks.ib.wL)
    end

    t = 0.0
    dt = KitBase.timestep(ks, ctr, t)
    res = zeros(3)
end

ptc = KitBase.init_ptc!(ks, ctr; mode=:soa, factor=2)
ptc_new = deepcopy(ptc)

