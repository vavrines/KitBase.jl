using KitBase, Plots

set = Setup(
    "gas", # matter
    "sod", # case
    "1d2f1v", # space
    "kfvs", # flux
    "bgk", # collision
    1, # species
    2, # interpolation order
    "vanleer", # limiter
    "fix", # boundary
    0.5, # cfl
    0.15, # simulation time
)
ps = PSpace1D(0.0, 1.0, 100, 1)
vs = VSpace1D(-5.0, 5.0, 72)
gas = Gas(1e-4, 0.0, 1.0, 2.0)
fw, ff, bc = ib_sod(set, ps, vs, gas)
ib = IB1F(fw, ff, bc)

ks = SolverSet(set, ps, vs, gas, ib)
ctr, face = init_fvm(ks, ks.ps, structarray = true)

t = 0.0
dt = KitBase.timestep(ks, ctr, t)
nt = ks.set.maxTime ÷ dt |> Int
for iter = 1:nt
    reconstruct!(ks, ctr)
    evolve!(ks, ctr, face, dt)
    update!(ks, ctr, face, dt, zeros(3))
end

plot(ks, ctr)
