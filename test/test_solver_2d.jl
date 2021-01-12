import KitBase

cd(@__DIR__)
ks, ctr, a1face, a2face, simTime = KitBase.initialize("config_2d2f.txt")

dt = KitBase.timestep(ks, ctr, simTime)
KitBase.reconstruct!(ks, ctr)
KitBase.evolve!(ks, ctr, a1face, a2face, dt; mode = :kfvs, bc = :maxwell)