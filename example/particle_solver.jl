using OffsetArrays, ProgressMeter
import KitBase

begin
    cd(@__DIR__)
    D = KitBase.read_dict("config.txt")
    for key in keys(D)
        s = Symbol(key)
        @eval $s = $(D[key])
    end

    γ = KitBase.heat_capacity_ratio(inK, 1)
    set = KitBase.Setup(case, space, flux, collision, nSpecies, interpOrder, limiter, cfl, maxTime)
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
    wL, primL, hL, bL, bcL, wR, primR, hR, bR, bcR = KitBase.ib_rh(gas.Ma, gas.γ, gas.K, vSpace.u)
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
                KitBase.vhs_collision_time(ks.ib.primR, ks.gas.μᵣ, ks.gas.ω),
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
res = zeros(3)

function solve!(ks, ctr, ptc, ptc_new, t, dt, res, nt=10)
    @showprogress for iter in 1:nt
    #for iter in 1:10
        #KitBase.reconstruct!(ks, ctr)
        
        @inbounds Threads.@threads for i = 1:ks.pSpace.nx+1
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
        
        KitBase.particle_transport!(ks, ctr, ptc, ptc_new, dt)
        KitBase.particle_collision!(ks, ctr, ptc_new, face, res, :bgk)
        KitBase.particle_boundary!(ks, ctr, ptc, ptc_new, face, dt, :bgk, :fix)
        KitBase.duplicate!(ptc, ptc_new, ks.gas.np)
        
        t += dt
    end
    @show res
    return t
end

t = solve!(ks, ctr, ptc, ptc_new, t, dt, res, 2)

KitBase.plot_line(ks, ctr)

#for i = 1:1
    KitBase.particle_transport!(ks, ctr, ptc, ptc_new, dt)
    KitBase.particle_collision!(ks, ctr, ptc_new, face, res, :bgk)
    KitBase.particle_boundary!(ks, ctr, ptc, ptc_new, face, dt, :bgk, :fix)
    KitBase.duplicate!(ptc, ptc_new, ks.gas.np)
#end

begin
    using Plots
    pltx = ks.pSpace.x[1:ks.pSpace.nx]
    plty = zeros(ks.pSpace.nx, 6)
    for i in eachindex(pltx)
        for j = 1:2
            plty[i, j] = ctr[i].wf[j] .+ ctr[i].wp[j]
        end
        plty[i, 3] = ctr[i].wf[end] + ctr[i].wp[end]
    end
    plot(pltx, plty[:, 1:3], lw = 2, legend=:none)
end
