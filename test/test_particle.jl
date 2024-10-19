KB.sample_maxwell(0.8, 1.0)
KB.sample_maxwell([1.0, 0.0, 1.0])
KB.sample_maxwell([1.0, 0.0, 0.0, 1.0])
KB.sample_maxwell([1.0, 0.0, 0.0, 0.0, 1.0])
KB.sample_maxwell(1.0, -2.0, 2.0)
KB.sample_maxwell(1.0, 0.3, -2.0, 2.0)
KB.sample_maxwell([1.0, 0.0, 1.0], -1.0, 1.0)
KB.sample_maxwell([1.0, 0.0, 0.0, 1.0], -1.0, 1.0)
KB.sample_maxwell([1.0, 0.0, 0.0, 0.0, 1.0], -1.0, 1.0)

KB.next_collision_time(1.0)

KB.boundary_time(1.0, 1.0, 0.5, 1.5)
KB.boundary_time(1.0, -1.0, 0.5, 1.5)
KB.boundary_time(1.0, 0.0, 0.5, 1.5)
KB.boundary_time(1.0, [1.0, 0.0, 0.0], 0.5, 1.5)
KB.boundary_time(1.0, [-1.0, 0.0, 0.0], 0.5, 1.5)
KB.boundary_time(1.0, [0.0, 0.0, 0.0], 0.5, 1.5)

ptc1 = KB.Particle1D(1e-3, 0.0, randn(3), 0.1, 1)
ptc2 = KB.Particle2D(1e-3, 0.0, 0.0, randn(3), 0.1, 1, 1)

KB.sample_particle!(ptc1, 1e-4, rand(), randn(3), rand(), 2, 0, 0.1)
KB.sample_particle!(ptc1, 1e-4, rand(), rand(), [1.0, 0.0, 1.0], 2, 1e-3, 0.72, 0)
KB.sample_particle!(
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
ks, ctr, face, simTime = KB.initialize("config.txt")

KB.sample_particle!(ptc1, ks, ctr[1], 1, 2.0, 0.05)

KB.Particle(
    ones(50) .* 1e-3,
    randn(50),
    randn(50, 3),
    rand(50),
    collect(1:50),
    zeros(Int64, 50),
    zeros(50),
)
KB.ControlVolumeParticle1D(
    rand(),
    rand(),
    KB.prim_conserve([1.0, 0.0, 1.0], 1.67),
    [1.0, 0.0, 1.0],
)
KB.ControlVolumeParticle2D(
    rand(),
    rand(),
    rand(),
    rand(),
    KB.prim_conserve([1.0, 0.0, 0.0, 1.0], 1.67),
    [1.0, 0.0, 0.0, 1.0],
)

using OffsetArrays

begin
    cd(@__DIR__)
    D = KB.read_dict("poiseuille.txt")
    for key in keys(D)
        s = Symbol(key)
        @eval $s = $(D[key])
    end

    γ = KB.heat_capacity_ratio(inK, 1)
    set = KB.Setup(
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
        hasForce,
    )
    pSpace = KB.PSpace1D(x0, x1, nx, nxg)
    vSpace = KB.VSpace1D(umin, umax, nu; type=vMeshType, ng=nug)
    μᵣ = KB.ref_vhs_vis(knudsen, alphaRef, omegaRef)
    gas = KB.Gas(knudsen, mach, prandtl, inK, γ, omega, alphaRef, omegaRef, μᵣ, mass, 0, ())

    primL = [1.0, 0.0, -1.0, 1.0] # left wall
    primR = [1.0, 0.0, 1.0, 1.0] # right wall
    wL = KB.prim_conserve(primL, γ)
    wR = KB.prim_conserve(primR, γ)

    fw = function (x, p)
        if x <= (pSpace.x0 + pSpace.x1) / 2
            wL
        else
            wR
        end
    end
    bc = function (x, p)
        if x <= (pSpace.x0 + pSpace.x1) / 2
            primL
        else
            primR
        end
    end
    ib = KB.IB(fw, bc, nothing)

    ks = KB.SolverSet(set, pSpace, vSpace, gas, ib, pwd())
end

begin
    ctr = OffsetArray{KB.ControlVolumeParticle1D}(undef, eachindex(ks.ps.x))
    face = Array{KB.Interface1D}(undef, ks.ps.nx + 1)
    for i in eachindex(ctr)
        prim = [1.0, 0.0, primL[3] + 2.0 * (ks.ps.x[i] - ks.ps.x0), 1.0]

        ctr[i] = KB.ControlVolumeParticle1D(
            ks.ps.x[i],
            ks.ps.dx[i],
            KB.prim_conserve(prim, ks.gas.γ),
            prim,
            KB.vhs_collision_time(prim, ks.gas.μᵣ, ks.gas.ω),
        )
    end

    for i in 1:ks.ps.nx+1
        face[i] = KB.Interface1D(ctr[1].w)
    end

    t = 0.0
    dt = KB.timestep(ks, ctr, t)
    res = zeros(4)
end

ptc = KB.init_ptc!(ks, ctr)

for iter in 1:1
    KB.free_transport!(ks, ptc.x, ptc.v, ptc.flag, dt)
    KB.boundary!(ks, ctr, ptc, face, dt, :maxwell)
    KB.sort!(ks, ctr, ptc.x, ptc.idx, ptc.ref)
    KB.dsmc!(ks, ctr, ptc.ref, ptc.v, dt)
    KB.stat!(ks, ctr, ptc)
end
