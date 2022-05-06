using KitBase, OffsetArrays, Plots
using KitBase.ProgressMeter: @showprogress

set = KitBase.Setup(
    "scalar", # matter
    "advection", # case
    "1d0f0v", # space
    "gks", # flux
    "", # collision: for scalar conservation laws there are none
    1, # species
    2, # interpolation order
    "vanleer", # limiter
    "period", # boundary
    0.5, # cfl
    1.0, # simulation time
)

pSpace = KitBase.PSpace1D(0.0, 1.0, 100, 1)
vSpace = nothing
property = KitBase.Scalar(0, 1e-6)

w0 = 1.0
prim0 = KitBase.conserve_prim(w0)
ib = nothing
ks = KitBase.SolverSet(set, pSpace, vSpace, property, ib, @__DIR__)

ctr = OffsetArray{KitBase.ControlVolume1D}(undef, eachindex(ks.pSpace.x))
for i in eachindex(ctr)
    u = sin(2π * ks.pSpace.x[i])
    ctr[i] = KitBase.ControlVolume1D(
        u,
        KitBase.conserve_prim(u),
    )
end

face = Array{KitBase.Interface1D}(undef, ks.pSpace.nx + 1)
for i = 1:ks.pSpace.nx+1
    face[i] = KitBase.Interface1D(0.0)
end

t = 0.0
dt = KitBase.timestep(ks, ctr, t)
nt = ks.set.maxTime ÷ dt |> Int

anim = @animate for iter = 1:nt
    KitBase.reconstruct!(ks, ctr)

    for i in eachindex(face)
        face[i].fw = KitBase.flux_gks(
            ctr[i-1].w,
            ctr[i].w,
            ks.gas.μᵣ,
            dt,
            0.5 * ks.ps.dx[i-1],
            0.5 * ks.ps.dx[i],
        )
    end

    for i = 1:ks.pSpace.nx
        ctr[i].w += (face[i].fw - face[i+1].fw) / ks.ps.dx[i]
        ctr[i].prim .= KitBase.conserve_prim(ctr[i].w)
    end
    ctr[0].w = ctr[ks.pSpace.nx].w
    ctr[ks.pSpace.nx+1].w = ctr[1].w

    sol = zeros(ks.pSpace.nx)
    for i = 1:ks.pSpace.nx
        sol[i] = ctr[i].w
    end
    plot(ks.pSpace.x[1:ks.pSpace.nx], sol, xlabel = "x", label = "u", ylims = [-1, 1])
end

gif(anim, "burgers.gif", fps = 45)
