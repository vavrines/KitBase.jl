cd(@__DIR__)
allowed = ["case", "space", "flux", "collision"]
D = KB.read_dict("config.txt", allowed)

D = KB.read_dict("config.txt")
ks, ctr, face, simTime = KB.initialize(D)

ks, ctr, face, simTime = KB.initialize("config.txt")
KB.write_sol(ks, ctr)
KB.write_sol(ks, ctr; mode=:jld)

KB.plot_line(ks, ctr)
plot(ks, ctr; legend=:none)

# tecplot writer
x = rand(5)
sol = randn(5, 2)
KB.write_tec(x, sol)

x = zeros(4, 3)
y = zero(x)
for i in axes(x, 1), j in axes(y, 2)
    x[i, j] = 0.25 * (i - 1)
    y[i, j] = 0.5 * (j - 1)
end
sol = randn((size(x)..., 2))
KB.write_tec(x, y, sol)
