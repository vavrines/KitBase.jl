cd(@__DIR__)
ks, ctr, a1face, a2face, simTime = KitBase.initialize("config_2d2f.txt")
simTime = KitBase.solve!(ks, ctr, a1face, a2face, simTime)
KitBase.extract_sol(ks, ctr)

res = zeros(4)
dt = KitBase.timestep(ks, ctr, 0.0)
KitBase.evolve!(ks, ctr, a1face, a2face, dt; mode = :kcu, bc = :maxwell)
KitBase.update!(ks, ctr, a1face, a2face, dt, zeros(4); coll = :bgk, bc = :extra)
KitBase.update!(ks, ctr, a1face, a2face, dt, zeros(4); coll = :bgk, bc = :period)
KitBase.update!(ks, ctr, a1face, a2face, dt, zeros(4); coll = :bgk, bc = :mirror)

# boundary
c1 = deepcopy(ctr[1])
KB.bc_riemann!(c1, ctr[1], ks.vs.u, ks.vs.v, [1.0, 0.5, 0.0, 1.0], 5 / 3, 2, [1.0 / √2, 1.0 / √2])

KB.extract_wall(c1, [1.0, 0.5, 0.0, 1.0], ks.vs, 1, 5/3, [1.0 / √2, 1.0 / √2], 1)

ks, ctr, a1face, a2face, simTime = KitBase.initialize("config_2d1f.txt")
KitBase.reconstruct!(ks, ctr)
KitBase.update!(ks, ctr, a1face, a2face, dt, zeros(4); coll = :bgk, bc = :maxwell)
KitBase.evolve!(ks, ctr, a1face, a2face, dt; mode = :kfvs, bc = :maxwell)
KitBase.evolve!(ks, ctr, a1face, a2face, dt; mode = :kcu, bc = :maxwell)
KitBase.update!(ks, ctr, a1face, a2face, dt, zeros(4); coll = :bgk, bc = :extra)
KitBase.update!(ks, ctr, a1face, a2face, dt, zeros(4); coll = :bgk, bc = :period)
KitBase.update!(ks, ctr, a1face, a2face, dt, zeros(4); coll = :bgk, bc = :mirror)

ks, ctr, a1face, a2face, simTime = KitBase.initialize("config_2d.txt")
KitBase.reconstruct!(ks, ctr)
KitBase.update!(ks, ctr, a1face, a2face, dt, zeros(4); coll = :bgk, bc = :maxwell)
KitBase.evolve!(ks, ctr, a1face, a2face, dt; mode = :hll, bc = :maxwell)
KitBase.update!(ks, ctr, a1face, a2face, dt, zeros(4); coll = :bgk, bc = :extra)
KitBase.update!(ks, ctr, a1face, a2face, dt, zeros(4); coll = :bgk, bc = :period)
KitBase.update!(ks, ctr, a1face, a2face, dt, zeros(4); coll = :bgk, bc = :mirror)

KitBase.plot_contour(ks, ctr)
plot(ks, ctr)
