using KitBase

D = Dict{Symbol,Any}()
begin
    D[:matter] = "gas"
    D[:case] = "cavity"
    D[:space] = "2d2f2v"
    D[:flux] = "kfvs"
    D[:collision] = "bgk"
    D[:nSpecies] = 1
    D[:interpOrder] = 2
    D[:limiter] = "vanleer"
    D[:boundary] = "maxwell"
    D[:cfl] = 0.8
    D[:maxTime] = 0.01

    D[:x0] = 0.0
    D[:x1] = 1.0
    D[:nx] = 45
    D[:y0] = 0.0
    D[:y1] = 1.0
    D[:ny] = 45
    D[:pMeshType] = "uniform"
    D[:nxg] = 0
    D[:nyg] = 0

    D[:umin] = -5.0
    D[:umax] = 5.0
    D[:nu] = 28
    D[:vmin] = -5.0
    D[:vmax] = 5.0
    D[:nv] = 28
    D[:vMeshType] = "rectangle"
    D[:nug] = 0
    D[:nvg] = 0

    D[:knudsen] = 0.075
    D[:mach] = 2.0
    D[:prandtl] = 1.0
    D[:inK] = 1.0
    D[:omega] = 0.72
    D[:alphaRef] = 1.0
    D[:omegaRef] = 0.5

    D[:uLid] = 0.15
    D[:vLid] = 0.0
    D[:tLid] = 1.0
end

ks = SolverSet(D)

ctr, a1face, a2face = KitBase.init_fvm(ks, ks.pSpace; structarray = true)
ctr0, a1face0, a2face0 = KitBase.init_fvm(ks, ks.pSpace; structarray = false)

using BenchmarkTools

@btime KitBase.solve!(ks, ctr, a1face, a2face, 0.0)
```~0.010509362```

@btime KitBase.solve!(ks, ctr0, a1face0, a2face0, 0.0)
```~0.010509326```
