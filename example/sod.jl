using KitBase, Plots
using KitBase.ProgressMeter: @showprogress

set = Setup(case = "sod", space = "1d1f1v", maxTime = 0.12)
ps = PSpace1D(0.0, 1.0, 100, 1)
vs = VSpace1D(-5.0, 5.0, 72)
gas = Gas(Kn = 1e-4, K = 0.0, γ = 3.0)
fw, ff, bc, p = config_ib(set, ps, vs, gas)
ib = IB1F(fw, ff, bc, p)

ks = SolverSet(set, ps, vs, gas, ib)
ctr, face = init_fvm(ks, :static_array; structarray = true)

t = 0.0
dt = KitBase.timestep(ks, ctr, t)
nt = ks.set.maxTime ÷ dt |> Int
@showprogress for iter = 1:nt
    reconstruct!(ks, ctr)
    evolve!(ks, ctr, face, dt)
    update!(ks, ctr, face, dt, zeros(3))
end

plot(ks, ctr)
