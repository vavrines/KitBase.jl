cd(@__DIR__)

config = read_cfg("config_3d.txt")
set = KitBase.set_setup(; config...)
ps = KitBase.set_geometry(; config...)
vs = KitBase.set_velocity(; config...)
gas = KitBase.set_property(config)

bc = (x, y, z, args...) -> [1.0, 1.0, 0.0, 0.0, 1.0]
fw = (x, y, z, args...) -> prim_conserve([1.0, 1.0, 0.0, 0.0, 1.0], gas.γ)

ib0 = IB(fw, bc, NamedTuple())

ks = SolverSet(set, ps, vs, gas, ib0)

ctr, a1face, a2face, a3face = init_fvm(ks; structarray = true)
simTime = 0.0
simTime = KB.solve!(ks, ctr, a1face, a2face, a3face, simTime)
KB.extract_sol(ks, ctr)

res = zeros(5)
dt = KB.timestep(ks, ctr, 0.0)
KB.update!(ks, ctr, a1face, a2face, a3face, dt, zeros(4); coll=:bgk, bc=:extra)
KB.update!(ks, ctr, a1face, a2face, a3face, dt, zeros(4); coll=:bgk, bc=:period)
KB.update!(ks, ctr, a1face, a2face, a3face, dt, zeros(4); coll=:bgk, bc=:mirror)
