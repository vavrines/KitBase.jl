cd(@__DIR__)
ks, ctr, face, simTime = KitBase.initialize("config_1d0f.txt")
dt = KitBase.timestep(ks, ctr, 0.0)
KitBase.update!(ks, ctr, face, dt, zeros(3); bc = :extra)
KitBase.update!(ks, ctr, face, dt, zeros(3); bc = :period)
KitBase.update!(ks, ctr, face, dt, zeros(3); bc = :balance)
t = KitBase.solve!(ks, ctr, face, 0.0)
KitBase.evolve!(ks, ctr, face, dt; mode = :roe)

ks, ctr, face, simTime = KitBase.initialize("config_1d1f.txt")
dt = KitBase.timestep(ks, ctr, 0.0)
KitBase.update!(ks, ctr, face, dt, zeros(3); bc = :extra)
KitBase.update!(ks, ctr, face, dt, zeros(3); bc = :period)
KitBase.update!(ks, ctr, face, dt, zeros(3); bc = :balance)
t = KitBase.solve!(ks, ctr, face, 0.0)
KitBase.evolve!(ks, ctr, face, dt; mode = :kcu)

ks, ctr, face, simTime = KitBase.initialize("config_1d2f.txt")
t = KitBase.solve!(ks, ctr, face, 0.0)

gas = KitBase.Gas(
    ks.gas.Kn,
    ks.gas.Ma,
    ks.gas.Pr,
    ks.gas.K,
    ks.gas.γ,
    ks.gas.ω,
    ks.gas.αᵣ,
    ks.gas.ωᵣ,
    ks.gas.μᵣ,
    1e-3,
    1000,
)
ks1 = KitBase.SolverSet(ks.set, ks.pSpace, ks.vSpace, gas, ks.ib, ks.outputFolder)
KitBase.init_ptc!(ks1, ctr, mode = :soa)
KitBase.init_ptc!(ks1, ctr, mode = :aos)