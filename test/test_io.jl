cd(@__DIR__)
allowed = ["case", "space", "flux", "collision"]
D = KitBase.read_dict("config.txt", allowed)

ks, ctr, face, simTime = KitBase.initialize("config.txt")
KitBase.write_jld(ks, ctr)

KitBase.plot_line(ks, ctr; backend = :plots)
""":gr mode is not suitable for remote test due to Qt and GKS issues"""