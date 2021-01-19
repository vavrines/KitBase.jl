KitBase.sample_maxwell(0.8, 1.0)
KitBase.sample_maxwell([1.0, 0.0, 1.0])
KitBase.sample_maxwell([1.0, 0.0, 0.0, 1.0])
KitBase.sample_maxwell([1.0, 0.0, 0.0, 0.0, 1.0])
KitBase.sample_maxwell(1.0, -2.0, 2.0)
KitBase.sample_maxwell(1.0, 0.3, -2.0, 2.0)
KitBase.sample_maxwell([1.0, 0.0, 1.0], -1.0, 1.0)
KitBase.sample_maxwell([1.0, 0.0, 0.0, 1.0], -1.0, 1.0)
KitBase.sample_maxwell([1.0, 0.0, 0.0, 0.0, 1.0], -1.0, 1.0)

KitBase.next_collision_time(1.0)

KitBase.boundary_time(1.0, 1.0, 0.5, 1.5)
KitBase.boundary_time(1.0, -1.0, 0.5, 1.5)
KitBase.boundary_time(1.0, 0.0, 0.5, 1.5)
KitBase.boundary_time(1.0, [1.0, 0.0, 0.0], 0.5, 1.5)
KitBase.boundary_time(1.0, [-1.0, 0.0, 0.0], 0.5, 1.5)
KitBase.boundary_time(1.0, [0.0, 0.0, 0.0], 0.5, 1.5)

ptc1 = KitBase.Particle1D(1e-3, 0.0, randn(3), 0.1, 1) |> show
ptc2 = KitBase.Particle2D(1e-3, 0.0, 0.0, randn(3), 0.1, 1, 1) |> show

KitBase.sample_particle!(ptc1, 1e-4, rand(), randn(3), rand(), 2, 0, 0.1)
KitBase.sample_particle!(ptc1, 1e-4, rand(), rand(), [1.0, 0.0, 1.0], 2, 1e-3, 0.72, 0)
KitBase.sample_particle!(
    ptc1,
    1e-4,
    rand(),
    rand(),
    [1.0, 0.0, 1.0],
    -1.0,
    1.0,
    2,
    1e-3,
    0.72,
    0,
)

cd(@__DIR__)
ks, ctr, face, simTime = KitBase.initialize("config.txt")

KitBase.sample_particle!(ptc1, ks, ctr[1], 1)

KitBase.Particle(
    ones(50) .* 1e-3,
    randn(50),
    randn(50, 3),
    rand(50),
    collect(1:50),
    zeros(Int64, 50),
    zeros(50),
)
KitBase.ControlVolumeParticle1D(
    rand(),
    rand(),
    KitBase.prim_conserve([1.0, 0.0, 1.0], 1.67),
    [1.0, 0.0, 1.0],
)
KitBase.ControlVolumeParticle2D(
    rand(),
    rand(),
    rand(),
    rand(),
    KitBase.prim_conserve([1.0, 0.0, 0.0, 1.0], 1.67),
    [1.0, 0.0, 0.0, 1.0],
)

using OffsetArrays

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

for iter = 1:1
    KitBase.free_transport!(ks, ptc.x, ptc.v, ptc.flag, dt)
    KitBase.boundary!(ks, ctr, ptc, face, dt, :maxwell)
    KitBase.sort!(ks, ctr, ptc.x, ptc.idx, ptc.ref)
    KitBase.dsmc!(ks, ctr, ptc.ref, ptc.v, dt)
    KitBase.stat!(ks, ctr, ptc)
end
