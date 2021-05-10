cd(@__DIR__)
ks, ctr, a1face, a2face, simTime = KitBase.initialize("config_2d2f.txt")
simTime = KitBase.solve!(ks, ctr, a1face, a2face, simTime)

res = zeros(4)
dt = KitBase.timestep(ks, ctr, simTime)
KitBase.evolve!(ks, ctr, a1face, a2face, dt; mode = :kcu, bc = :maxwell)

ks, ctr, a1face, a2face, simTime = KitBase.initialize("config_2d1f.txt")
KitBase.reconstruct!(ks, ctr)
KitBase.update!(ks, ctr, a1face, a2face, dt, zeros(4); coll = :bgk, bc = :maxwell)
KitBase.evolve!(ks, ctr, a1face, a2face, dt; mode = :kfvs, bc = :maxwell)
KitBase.evolve!(ks, ctr, a1face, a2face, dt; mode = :kcu, bc = :maxwell)

KitBase.plot_contour(ks, ctr)
