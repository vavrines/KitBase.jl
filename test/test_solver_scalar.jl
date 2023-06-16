# scalar
set = Setup(
    "scalar", # matter
    "advection", # case
    "1d0f0v", # space
    "gks", # flux
    "", # collision: for scalar conservation laws there are none
    1, # species
    1, # interpolation order
    "vanleer", # limiter
    "period", # boundary
    0.5, # cfl
    1.0, # simulation time
    false,
)
pSpace = PSpace1D(0.0, 1.0, 100, 1)
vSpace = nothing
property = Scalar(1.0, 1e-6)
w0 = 1.0
prim0 = conserve_prim(w0, property.a)
ib = IB((x, y) -> sin(2π * x), nothing, nothing)
ks = SolverSet(set, pSpace, vSpace, property, ib)
SolverSet(set, pSpace, property)

ctr, face = KitBase.init_fvm(ks)
for i in eachindex(ctr)
    ctr[i].w = sin(2π * ks.ps.x[i])
end

dt = 1e-3
reconstruct!(ks, ctr)
evolve!(ks, ctr, face, dt)
update!(ks, ctr, face, dt, 0.0)
