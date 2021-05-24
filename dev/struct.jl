using StructArrays, ProgressMeter, KitBase

cd(@__DIR__)
D = Dict{Symbol,Any}()
begin
    D[:matter] = "gas"
    D[:case] = "sod"
    D[:space] = "1d1f1v"
    D[:flux] = "kfvs"
    D[:collision] = "bgk"
    D[:nSpecies] = 1
    D[:interpOrder] = 2
    D[:limiter] = "vanleer"
    D[:boundary] = "fix"
    D[:cfl] = 0.5
    D[:maxTime] = 5.0

    D[:x0] = 0.0
    D[:x1] = 1.0
    D[:nx] = 100
    D[:pMeshType] = "uniform"
    D[:nxg] = 0

    D[:umin] = -5.0
    D[:umax] = 5.0
    D[:nu] = 64
    D[:vMeshType] = "rectangle"
    D[:nug] = 0

    D[:knudsen] = 0.01
    D[:mach] = 0.0
    D[:prandtl] = 1.0
    D[:inK] = 0.0
    D[:omega] = 0.81
    D[:alphaRef] = 1.0
    D[:omegaRef] = 0.5
end

ks = SolverSet(D)

ctr, face = init_fvm(ks, ks.ps)



ctr = Array{ControlVolume1D1F}(undef, ks.pSpace.nx)
for i in eachindex(ctr)
    prim = [10.0 * rand(), randn(), 1 / rand()]
    w = prim_conserve(prim, ks.gas.γ)
    g = maxwellian(ks.vs.u, prim)

    ctr[i] = ControlVolume1D1F(
        ks.pSpace.x[i],
        ks.pSpace.dx[i],
        w,
        prim,
        g,
    )
end
ctr = ctr |> StructArray

face = Array{Interface1D1F}(undef, ks.pSpace.nx + 1)
for i = 1:ks.pSpace.nx+1
    face[i] = Interface1D1F(deepcopy(ks.ib.wL), deepcopy(ks.ib.fL))
    #face[i] = Interface1D1F(ks.ib.wL, ks.ib.fL)
end
face = face |> StructArray

plot_line(ks, ctr)

t = 0.0
dt = timestep(ks, ctr, t)
nt = Int(ks.set.maxTime ÷ dt) + 1
res = zero(ks.ib.wL)

@showprogress for iter = 1:100#nt
    reconstruct!(ks, ctr)
    evolve!(ks, ctr, face, dt; mode = Symbol(ks.set.flux), bc = Symbol(ks.set.boundary))
    KitBase.update!(ks, ctr, face, dt, res; coll = Symbol(ks.set.collision), bc = Symbol(ks.set.boundary))

    t += dt
    #if t > ks.set.maxTime || maximum(res) < 5.e-7
    #    break
    #end
end

plot_line(ks, ctr)


###

cd(@__DIR__)
D = Dict{Symbol,Any}()
begin
    D[:matter] = "gas"
    D[:case] = "sod"
    D[:space] = "1d2f1v"
    D[:flux] = "kfvs"
    D[:collision] = "bgk"
    D[:nSpecies] = 1
    D[:interpOrder] = 2
    D[:limiter] = "vanleer"
    D[:boundary] = "fix"
    D[:cfl] = 0.5
    D[:maxTime] = 5.0

    D[:x0] = 0.0
    D[:x1] = 1.0
    D[:nx] = 100
    D[:pMeshType] = "uniform"
    D[:nxg] = 0

    D[:umin] = -5.0
    D[:umax] = 5.0
    D[:nu] = 64
    D[:vMeshType] = "rectangle"
    D[:nug] = 0

    D[:knudsen] = 0.01
    D[:mach] = 0.0
    D[:prandtl] = 1.0
    D[:inK] = 2.0
    D[:omega] = 0.81
    D[:alphaRef] = 1.0
    D[:omegaRef] = 0.5
end

ks = SolverSet(D)

ctr, face = init_fvm(ks, ks.ps)



ctr = Array{ControlVolume1D2F}(undef, ks.pSpace.nx)
for i in eachindex(ctr)
    prim = [10.0 * rand(), randn(), 1 / rand()]
    w = prim_conserve(prim, ks.gas.γ)
    g = maxwellian(ks.vs.u, prim)
    b = g ./ prim[end]

    ctr[i] = ControlVolume1D2F(
        ks.pSpace.x[i],
        ks.pSpace.dx[i],
        w,
        prim,
        g,
        b,
    )
end
ctr = ctr |> StructArray

face = Array{Interface1D2F}(undef, ks.pSpace.nx + 1)
for i = 1:ks.pSpace.nx+1
    #face[i] = Interface1D1F(deepcopy(ks.ib.wL), deepcopy(ks.ib.fL))
    face[i] = Interface1D2F(ks.ib.wL, ks.ib.hL)
end
face = face |> StructArray

plot_line(ks, ctr)

t = 0.0
dt = timestep(ks, ctr, t)
nt = Int(ks.set.maxTime ÷ dt) + 1
res = zero(ks.ib.wL)

@showprogress for iter = 1:100#nt
    reconstruct!(ks, ctr)
    evolve!(ks, ctr, face, dt; mode = Symbol(ks.set.flux), bc = Symbol(ks.set.boundary))
    KitBase.update!(ks, ctr, face, dt, res; coll = Symbol(ks.set.collision), bc = Symbol(ks.set.boundary))

    t += dt
    #if t > ks.set.maxTime || maximum(res) < 5.e-7
    #    break
    #end
end

plot_line(ks, ctr)