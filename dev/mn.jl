using KitBase, Plots

set = Setup(space = "1d1f1v", boundary = "period", maxTime = 0.5)
ps = PSpace1D(0, 1, 100, 1)
vs = VSpace1D(-5, 5, 100)
gas = Gas(Kn = 1e-3, K = 0.0, γ = 3.0)
ib = begin
    fw = function (x)
        ρ = 1 + 0.1 * sin(2π * x)
        T = 1 / ρ
        prim = [ρ, 1.0, 1 / T]
        return prim_conserve(prim, gas.γ)
    end
    ff = function (x)
        w = fw(x)
        prim = conserve_prim(w, gas.γ)
        return maxwellian(vs.u, prim)
    end
    bc = function (x)
        w = fw(x)
        return conserve_prim(w, gas.γ)
    end
    IB1F(fw, ff, bc)
end
ks = SolverSet(set, ps, vs, gas, ib)

#=ctr, face = init_fvm(ks)
plot(ks, ctr)

t = 0.0
dt = timestep(ks, ctr, t)
nt = ks.set.maxTime ÷ dt |> Int
@showprogress for iter = 1:nt
    evolve!(ks, ctr, face, dt)
    update!(ks, ctr, face, dt, zeros(3))
end

plot!(ks, ctr)=#

# basis
m = zeros(5, length(vs.u))
m = zeros(7, length(vs.u))
for j in axes(m, 2)
    for i in axes(m, 1)
        m[i, j] = vs.u[j]^(i-1)
    end
end

# moments
sol = zeros(axes(ps.x, 1), size(m, 1))
for i in axes(sol, 1)
    f = ff(ps.x[i])
    for j in axes(sol, 2)
        sol[i, j] = sum(vs.weights .* m[j, :] .* f)
    end
end

# minimizer
α = zeros(size(m, 1))

# test
#f0 = maxwellian(vs.u, [1.0, 0.0, 1.0])
f0 = 0.4 .* maxwellian(vs.u, [1.0, 1.0, 1.0]) .+
    0.6 .* maxwellian(vs.u, [1.0, -1.0, 1.0])

u0 = [sum(vs.weights .* m[i, :] .* f0) for i in axes(m, 1)]
α = zeros(size(m, 1))
res = KitBase.optimize_closure(α, m, vs.weights, u0, exp)
f = exp.(res.minimizer' * m)'

plot(vs.u, f0)
plot!(vs.u, f)
