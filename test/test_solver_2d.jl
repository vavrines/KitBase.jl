cd(@__DIR__)
ks, ctr, a1face, a2face, simTime = KB.initialize("config_2d2f.txt")
simTime = KB.solve!(ks, ctr, a1face, a2face, simTime)
KB.extract_sol(ks, ctr)

res = zeros(4)
dt = KB.timestep(ks, ctr, 0.0)
KB.evolve!(ks, ctr, a1face, a2face, dt; mode = :kcu, bc = :maxwell)
KB.update!(ks, ctr, a1face, a2face, dt, zeros(4); coll = :bgk, bc = :extra)
KB.update!(ks, ctr, a1face, a2face, dt, zeros(4); coll = :bgk, bc = :period)
KB.update!(ks, ctr, a1face, a2face, dt, zeros(4); coll = :bgk, bc = :mirror)

# boundary
c1 = deepcopy(ctr[1])
KB.bc_riemann!(
    c1,
    ctr[1],
    ks.vs.u,
    ks.vs.v,
    [1.0, 0.5, 0.0, 1.0],
    5 / 3,
    2,
    [1.0 / √2, 1.0 / √2],
)

KB.extract_wall(c1, [1.0, 0.5, 0.0, 1.0], ks.vs, 1, 5 / 3, [1.0 / √2, 1.0 / √2], 1)

ks, ctr, a1face, a2face, simTime = KB.initialize("config_2d1f.txt")
KB.reconstruct!(ks, ctr)
KB.update!(ks, ctr, a1face, a2face, dt, zeros(4); coll = :bgk, bc = :maxwell)
KB.evolve!(ks, ctr, a1face, a2face, dt; mode = :kfvs, bc = :maxwell)
KB.evolve!(ks, ctr, a1face, a2face, dt; mode = :kcu, bc = :maxwell)
KB.update!(ks, ctr, a1face, a2face, dt, zeros(4); coll = :bgk, bc = :extra)
KB.update!(ks, ctr, a1face, a2face, dt, zeros(4); coll = :bgk, bc = :period)
KB.update!(ks, ctr, a1face, a2face, dt, zeros(4); coll = :bgk, bc = :mirror)

ks, ctr, a1face, a2face, simTime = KB.initialize("config_2d.txt")
KB.reconstruct!(ks, ctr)
KB.update!(ks, ctr, a1face, a2face, dt, zeros(4); coll = :bgk, bc = :maxwell)
KB.evolve!(ks, ctr, a1face, a2face, dt; mode = :hll, bc = :maxwell)
KB.update!(ks, ctr, a1face, a2face, dt, zeros(4); coll = :bgk, bc = :extra)
KB.update!(ks, ctr, a1face, a2face, dt, zeros(4); coll = :bgk, bc = :period)
KB.update!(ks, ctr, a1face, a2face, dt, zeros(4); coll = :bgk, bc = :mirror)

KB.plot_contour(ks, ctr)
plot(ks, ctr)
