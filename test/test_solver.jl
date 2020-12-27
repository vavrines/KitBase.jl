cd(@__DIR__)
ks, ctr, face, simTime = KitBase.initialize("config.txt")

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
ptc = KitBase.init_ptc!(ks1, ctr)
