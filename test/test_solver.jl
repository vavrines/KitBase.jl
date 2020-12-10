cd(@__DIR__)
ks, ctr, face, simTime = KitBase.initialize("config.txt")

particle = Particle(
    ks.gas.Kn, 
    ks.gas.Ma, 
    ks.gas.Pr,
    ks.gas.K,
    ks.gas.γ,
    1e-3,
    ks.gas.ω,
    ks.gas.αᵣ,
    ks.gas.ωᵣ,
    ks.gas.μᵣ,
)

ks1 = SolverSet(ks.set, ks.pSpace, ks.vSpace, particle, ks.ib, ks.outputFolder)
KitBase.init_ptc(ks1, ctr)