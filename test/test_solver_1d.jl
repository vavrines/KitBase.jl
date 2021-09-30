ks, ctr, face, simTime = KitBase.initialize("config_1d0f.txt")

using KitBase.JLD2
set = ks
t = 0.0
cd(@__DIR__)
@save "restart.jld2" set ctr t
ks, ctr, face, simTime = KitBase.initialize("restart.jld2")

dt = KitBase.timestep(ks, ctr, 0.0)
KitBase.update!(ks, ctr, face, dt, zeros(3); bc = :period)
KitBase.update!(ks, ctr, face, dt, zeros(3); bc = :extra)
KitBase.update!(ks, ctr, face, dt, zeros(3); bc = :balance)
KitBase.update!(ks, ctr, face, dt, zeros(3); bc = :mirror)
t = KitBase.solve!(ks, ctr, face, 0.0)
KitBase.evolve!(ks, ctr, face, dt; mode = :roe)

ks, ctr, face, simTime = KitBase.initialize("config_1d0f2s.txt")
dt = KitBase.timestep(ks, ctr, 0.0)
KitBase.update!(ks, ctr, face, dt, zeros(3, 2); bc = :extra)
KitBase.update!(ks, ctr, face, dt, zeros(3, 2); bc = :period)
KitBase.update!(ks, ctr, face, dt, zeros(3, 2); bc = :balance)
KitBase.update!(ks, ctr, face, dt, zeros(3, 2); bc = :mirror)
KitBase.evolve!(ks, ctr, face, dt; mode = :gks)

ks, ctr, face, simTime = KitBase.initialize("config_1d1f1v.txt")
KitBase.reconstruct!(ks, ctr)
KitBase.update!(ks, ctr, face, dt, zeros(3); bc = :extra)
KitBase.update!(ks, ctr, face, dt, zeros(3); bc = :period)
KitBase.update!(ks, ctr, face, dt, zeros(3); bc = :balance)
KitBase.update!(ks, ctr, face, dt, zeros(3); bc = :mirror)
KitBase.evolve!(ks, ctr, face, dt; mode = :kfvs)
KitBase.evolve!(ks, ctr, face, dt; mode = :kcu)

ks, ctr, face, simTime = KitBase.initialize("config_1d1f3v.txt")
KitBase.reconstruct!(ks, ctr)
KitBase.update!(ks, ctr, face, dt, zeros(5); bc = :extra)
KitBase.update!(ks, ctr, face, dt, zeros(5); bc = :period)
KitBase.update!(ks, ctr, face, dt, zeros(5); bc = :balance)
KitBase.update!(ks, ctr, face, dt, zeros(5); bc = :mirror)
KitBase.evolve!(ks, ctr, face, dt; mode = :kfvs)
KitBase.evolve!(ks, ctr, face, dt; mode = :kcu)

ks, ctr, face, simTime = KitBase.initialize("config_1d2f.txt")
KitBase.reconstruct!(ks, ctr)
KitBase.update!(ks, ctr, face, dt, zeros(3); bc = :extra)
KitBase.update!(ks, ctr, face, dt, zeros(3); bc = :period)
KitBase.update!(ks, ctr, face, dt, zeros(3); bc = :balance)
KitBase.update!(ks, ctr, face, dt, zeros(3); bc = :mirror)
KitBase.evolve!(ks, ctr, face, dt; mode = :kfvs)
KitBase.evolve!(ks, ctr, face, dt; mode = :kcu)

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

ks, ctr, face, simTime = KitBase.initialize("config_1d2f2s.txt")
KitBase.reconstruct!(ks, ctr)
KitBase.update!(ks, ctr, face, dt, zeros(3); bc = :extra)
KitBase.update!(ks, ctr, face, dt, zeros(3); bc = :period)
KitBase.update!(ks, ctr, face, dt, zeros(3); bc = :balance)
KitBase.update!(ks, ctr, face, dt, zeros(3); bc = :mirror)
KitBase.evolve!(ks, ctr, face, dt; mode = :kfvs)
KitBase.evolve!(ks, ctr, face, dt; mode = :kcu)

ks, ctr, face, simTime = KitBase.initialize("config_1d4f2s.txt")
KitBase.reconstruct!(ks, ctr)
KitBase.update!(ks, ctr, face, dt, zeros(3); bc = :extra)
KitBase.update!(ks, ctr, face, dt, zeros(3); bc = :period)
KitBase.update!(ks, ctr, face, dt, zeros(3); bc = :balance)
KitBase.update!(ks, ctr, face, dt, zeros(3); bc = :mirror)
KitBase.evolve!(ks, ctr, face, dt; mode = :kfvs)
KitBase.evolve!(ks, ctr, face, dt; mode = :kcu)

ks, ctr, face, simTime = KitBase.initialize("config_1d3f2s.txt")
KitBase.reconstruct!(ks, ctr)
KitBase.update!(ks, ctr, face, dt, zeros(3); bc = :extra)
KitBase.update!(ks, ctr, face, dt, zeros(3); bc = :period)
KitBase.update!(ks, ctr, face, dt, zeros(3); bc = :balance)
KitBase.update!(ks, ctr, face, dt, zeros(3); bc = :mirror)
KitBase.evolve!(ks, ctr, face, dt; mode = :kfvs)
KitBase.evolve!(ks, ctr, face, dt; mode = :kcu)
