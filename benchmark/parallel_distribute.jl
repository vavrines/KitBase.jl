using Distributed, SharedArrays, Plots
addprocs(1)
@everywhere using KitBase

begin
    vars = Dict{String,Any}()
    vars["case"] = "sod"
    vars["space"] = "1d0f0v"
    vars["flux"] = "kfvs"
    vars["collision"] = "bgk"
    vars["nSpecies"] = 1
    vars["interpOrder"] = 1
    vars["limiter"] = "vanleer"
    vars["boundary"] = "fix"
    vars["cfl"] = 0.5
    vars["maxTime"] = 0.2
    vars["x0"] = 0.0
    vars["x1"] = 1.0
    vars["nx"] = 1000
    vars["pMeshType"] = "uniform"
    vars["nxg"] = 0
    vars["knudsen"] = 0.001
    vars["mach"] = 0.0
    vars["prandtl"] = 1.0
    vars["inK"] = 0.0
    vars["omega"] = 0.81
    vars["alphaRef"] = 1.0
    vars["omegaRef"] = 0.5
end

for key in keys(vars)
    s = Symbol(key)
    @eval $s = $(vars[key])
end

set = Setup(
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


pSpace = KitBase.set_geometry(vars)
vSpace = KitBase.set_velocity(vars)
gas = KitBase.set_property(vars)
ib = KitBase.set_ib(vars, set, vSpace, gas)
folder = @__DIR__

ks = KitBase.SolverSet(set, pSpace, vSpace, gas, ib, folder)


wp = SharedArray{Float64}((ks.pSpace.nx, 3), init=A->(A=zeros(ks.pSpace.nx, 3)))

