cd(@__DIR__)
allowed = ["case", "space", "flux", "collision"]
D = KitBase.read_dict("config.txt", allowed)

D = KitBase.read_dict("config.txt")
ks, ctr, face, simTime = KitBase.initialize(D)

ks, ctr, face, simTime = KitBase.initialize("config.txt")
KitBase.write_jld(ks, ctr)

KitBase.plot_line(ks, ctr)
plot(ks, ctr, legend = :none)

# tecplot writer
x = zeros(5, 3)
y = zero(x)
for i = 1:5, j = 1:3
    x[i, j] = 0.25 * (i-1)
    y[i, j] = 0.5 * (j-1)
end
sol = randn(4, 2, 2)

write_tec(x, y, sol)
