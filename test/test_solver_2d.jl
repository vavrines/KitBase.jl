import KitBase

cd(@__DIR__)
ks, ctr, a1face, a2face, simTime = KitBase.initialize("config_2d2f.txt")

res = zeros(4)
dt = KitBase.timestep(ks, ctr, simTime)

for iter = 1:2
    KitBase.reconstruct!(ks, ctr)
    KitBase.evolve!(ks, ctr, a1face, a2face, dt; mode = :kfvs, bc = :maxwell)
    KitBase.update!(ks, ctr, a1face, a2face, dt, res; coll = :bgk, bc = :maxwell)
end