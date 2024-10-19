using KitBase

np = 10000 # number of particles
ptc = Array{Particle1D}(undef, np)

for i in 1:np
    rd1 = rand(3)
    rd2 = rand(3)
    rd = rand()

    m = 1e-5
    v = @. sqrt(-log(rd1) / 1.0) * sin(2.0 * π * rd2)
    x = 0.0 + (rd - 0.5) * 1.0
    idx = i

    if v[1] < 0
        tb = (0.0 - 0.5 * 1.0 - x) / v[1]
    elseif v[1] > 0
        tb = (0.0 + 0.5 * 1.0 - x) / v[1]
    else
        tb = 1.0
    end

    e = 0.0
    flag = 1

    ptc[i] = Particle1D(m, x, v, e, idx, flag, tb)
end

v = zeros(np, 3)
for i in 1:np
    v[i, :] .= ptc[i].v
end

using Plots
histogram(v[:, 1]; label="u")
histogram!(v[:, 2]; label="v")
histogram!(v[:, 3]; label="w")
