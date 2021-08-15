cd(@__DIR__)
allowed = ["case", "space", "flux", "collision"]
D = KitBase.read_dict("config.txt", allowed)

D = KitBase.read_dict("config.txt")
ks, ctr, face, simTime = KitBase.initialize(D)

ks, ctr, face, simTime = KitBase.initialize("config.txt")
KitBase.write_jld(ks, ctr)

using KitBase.Plots
a = plot(ks, ctr, legend=:none)
b = plot_line(ks, ctr)
