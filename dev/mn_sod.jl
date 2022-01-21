using KitBase, Plots, OffsetArrays
using ProgressMeter: @showprogress
using Base.Threads: @threads

set = Setup(space = "1d1f1v", boundary = "fix", maxTime = 0.05)
ps = PSpace1D(0, 1, 100, 1)
vs = VSpace1D(-5, 5, 36)
gas = Gas(Kn = 1e-3, K = 0.0, γ = 3.0)
ib = begin
    _fw = function (x)
        # wave
        #ρ = 1 + 0.1 * sin(2π * x)
        #U = 1.0
        #T = 1 / ρ

        # sod
        ρ, U, T = begin
            if x < 0.5
                1.0, 0.0, 1 / 0.5
            else
                0.125, 0.0, 1 / 0.625
            end
        end
        prim = [ρ, U, 1 / T]
        return prim_conserve(prim, gas.γ)
    end
    _ff = function (x)
        w = _fw(x)
        prim = conserve_prim(w, gas.γ)
        return maxwellian(vs.u, prim)
    end
    _bc = function (x)
        w = _fw(x)
        return conserve_prim(w, gas.γ)
    end
    IB1F(_fw, _ff, _bc)
end
ks = SolverSet(set, ps, vs, gas, ib)
#=
ctr, face = init_fvm(ks)
plot(ks, ctr)

t = 0.0
dt = timestep(ks, ctr, t)
nt = ks.set.maxTime ÷ dt |> Int
@showprogress for iter = 1:nt
    evolve!(ks, ctr, face, dt)
    update!(ks, ctr, face, dt, zeros(3))
end

plot!(ks, ctr)
=#

# basis
#m = zeros(5, length(vs.u))
m = zeros(7, length(vs.u))
for j in axes(m, 2)
    for i in axes(m, 1)
        m[i, j] = vs.u[j]^(i-1)
    end
end
α = zeros(size(m, 1)) # minimizer

###
# one-cell test
###
#f0 = maxwellian(vs.u, [1.0, 0.0, 1.0])
#=f0 = 0.4 .* maxwellian(vs.u, [1.0, 1.0, 1.0]) .+
    0.6 .* maxwellian(vs.u, [1.0, -1.0, 1.0])

u0 = [sum(vs.weights .* m[i, :] .* f0) for i in axes(m, 1)]
α = zeros(size(m, 1))
res = KitBase.optimize_closure(α, m, vs.weights, u0, exp)
f = exp.(res.minimizer' * m)'

plot(vs.u, f0)
plot!(vs.u, f)=#

ctr = OffsetArray{ControlVolume1F}(undef, eachindex(ks.ps.x))
for i in eachindex(ctr)
    cons = ks.ib.fw(ks.ps.x[i])
    prim = conserve_prim(cons, ks.gas.γ)
    f = ks.ib.ff(ks.ps.x[i])

    w = zeros(size(m, 1))
    for j in axes(w, 1)
        w[j] = sum(ks.vs.weights .* m[j, :] .* f)
    end

    ctr[i] = ControlVolume(w, prim, f, 1)
end

face = Array{Interface1F}(undef, ks.ps.nx + 1)
for i = 1:ks.ps.nx+1
    fw = zeros(size(m, 1))
    ff = zero(ks.vs.u)
    face[i] = Interface(fw, ff, 1)
end

function rc!(ks, ctr, p)
    m, α = p

    @inbounds @threads for i = 0:ks.ps.nx+1
        res = KitBase.optimize_closure(α, m, ks.vs.weights, ctr[i].w, exp)
        ctr[i].f .= exp.(res.minimizer' * m)'
    end
end

function ev!(ks, ctr, face, p)
    m, dt = p

    @inbounds @threads for i = 1:ks.ps.nx+1
        flux_kfvs!(face[i].ff, ctr[i-1].f, ctr[i].f, ks.vs.u, dt)
        for j in axes(face[i].fw, 1)
            face[i].fw[j] = sum(ks.vs.weights .* m[j, :] .* face[i].ff)
        end
    end
end

function up!(ks, ctr, face, p)
    m, dt = p

    @inbounds @threads for i = 1:ks.ps.nx
        m = maxwellian(ks.vs.u, ctr[i].prim)
        τ = vhs_collision_time(ctr[i].prim, ks.gas.μᵣ, ks.gas.ω)
        q = @. (m - ctr[i].f) / τ
        Q = [sum(ks.vs.weights .* m[j, :] .* q) for j = 1:7]

        @. ctr[i].w += (face[i].fw - face[i+1].fw) / ks.ps.dx[i] + Q * dt
        ctr[i].prim .= conserve_prim([ctr[i].w[1:2]; ctr[i].w[3] / 2], ks.gas.γ)
    end
end

t = 0.0
dt = timestep(ks, ctr, t)
nt = ks.set.maxTime ÷ dt |> Int

@showprogress for iter = 1:nt
    rc!(ks, ctr, (m, α))
    ev!(ks, ctr, face, (m, dt))
    up!(ks, ctr, face, (m, dt))
end

plot(ks, ctr)
