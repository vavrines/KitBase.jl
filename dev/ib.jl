using KitBase, Plots
using Base.Threads: @threads
using KitBase.ProgressMeter: @showprogress
using PyCall
itp = pyimport("scipy.interpolate")

set = Setup(;
    case="cylinder",
    space="2d0f0v",
    boundary=["fix", "extra", "mirror", "extra"],
    limiter="minmod",
    cfl=0.1,
    maxTime=2.0, # time
    flux="hll",
    hasForce=true,
)
ps = PSpace2D(0, 3, 60, 0, 2, 40, 1, 1)
vs = VSpace2D(-3, 3, 28, -3, 3, 28)
gas = Gas(; Kn=1e-3, Ma=0.1, K=1.0)

prim0 = [1.0, 0.0, 0.0, 1.0]
prim1 = [1.0, gas.Ma * sound_speed(1.0, gas.γ), 0.0, 1.0]
fw = (args...) -> prim_conserve(prim1, gas.γ)
ff = function (args...)
    prim = conserve_prim(fw(args...), gas.γ)
    h = maxwellian(vs.u, vs.v, prim)
    b = h .* gas.K / 2 / prim[end]
    return h, b
end
bc = function (x, y, args...)
    if abs((x - 3)^2 + y^2 - 1) < 1e-3
        return prim0
    else
        return prim1
    end
end
ib = IB2F(fw, ff, bc, NamedTuple())

#ks = SolverSet(set, ps, vs, gas, ib)
ks = SolverSet(set, ps, nothing, gas, ib)
ctr, a1face, a2face = init_fvm(ks; structarray=true)

radius = 1
θs = linspace(0, π / 2, 30)
xs = @. ps.x1 - radius * cos(θs)
ys = @. sin(θs)
lps = hcat(xs, ys)
h = ps.dx[1]
ΔV = 2π * radius / 4 / length(θs) * h

scatter(xs, ys; xticks=0:0.05:3, yticks=0:0.05:1, ratio=1)

struct IB{T1,T2,T3,TF,ND}
    xlp::T1 # position of Lagrangian points
    lpid::T2
    nlpid::T3
    ΔV::TF
end

ids = ones(Int, 30, 2)
for i in axes(ids, 1)
    ids[i, 1] = ceil(xs[i] / ps.dx[1]) |> Int
    ids[i, 2] = ceil(ys[i] / ps.dy[1]) |> Int
    if ids[i, 1] == 0
        ids[i, 1] += 1
    end
    if ids[i, 1] == ps.nx + 1
        ids[i, 1] -= 1
    end
    if ids[i, 2] == 0
        ids[i, 2] += 1
    end
    if ids[i, 2] == ps.ny + 1
        ids[i, 2] -= 1
    end
end
ids = [ids[i, :] for i in axes(ids, 1)]

neighbor_ids = [[ids[i]] for i in axes(ids, 1)]
for i in axes(ids, 1)
    push!(neighbor_ids[i], [ids[i][1] - 1, ids[i][2]])
    push!(neighbor_ids[i], [ids[i][1] + 1, ids[i][2]])
    push!(neighbor_ids[i], [ids[i][1], ids[i][2] - 1])
    push!(neighbor_ids[i], [ids[i][1], ids[i][2] + 1])
end

ib = IB{typeof(lps),typeof(ids),typeof(neighbor_ids),Float64,2}(lps, ids, neighbor_ids, ΔV)

t = 0.0
dt = timestep(ks, ctr, 0.0)
nt = ks.set.maxTime ÷ dt |> Int
res = zeros(4)

function ib_force!(ps, ctr, lps)
    for i in axes(ctr, 1), j in axes(ctr, 2)
        ctr[i, j].a .= 0.0
    end

    ufield = zero(ps.x)
    vfield = zero(ps.x)
    for i in axes(ufield, 1), j in axes(ufield, 2)
        ufield[i, j] = ctr[i, j].prim[2]
        vfield[i, j] = ctr[i, j].prim[3]
    end

    ucurve = itp.RegularGridInterpolator((ps.x[:, 1], ps.y[1, :]), ufield)
    vcurve = itp.RegularGridInterpolator((ps.x[:, 1], ps.y[1, :]), vfield)
    vls = hcat(ucurve(lps), vcurve(lps))

    fls = zero(vls)
    for i in axes(fls, 1)
        fls[i, 1] = -vls[i, 1] / dt
        fls[i, 2] = -vls[i, 2] / dt
    end

    for i in axes(ctr, 1), j in axes(ctr, 2)
        for k in axes(fls, 1)
            δx = ps.x[i, j] - xs[k, 1]
            δy = ps.y[i, j] - ys[k, 1]
            δ = KB.dirac_delta(δx, δy, h, h)
            ctr[i, j].a[1] += fls[k, 1] * δ * ΔV
            ctr[i, j].a[2] += fls[k, 2] * δ * ΔV
        end
    end

    return nothing
end

ib_force!(ks.ps, ctr, ib.xlp)

function stepf!(w, prim, fwL, fwR, fwD, fwU, γ, Δs, RES, AVG, f, dt)
    w_old = deepcopy(w)

    @. w += (fwL - fwR + fwD - fwU) / Δs
    prim .= conserve_prim(w, γ)

    # force
    prim[2] += f[1] * dt
    prim[3] += f[2] * dt
    w .= prim_conserve(prim, γ)

    @. RES += (w - w_old)^2
    @. AVG += abs(w)
end

function updatef!(KS, ctr, a1face, a2face, dt, residual;)
    nx, ny, dx, dy = KS.ps.nx, KS.ps.ny, KS.ps.dx, KS.ps.dy

    sumRes = zero(ctr[1].w)
    sumAvg = zero(ctr[1].w)

    @inbounds @threads for j in 1:ny
        for i in 1:nx
            stepf!(
                ctr[i, j].w,
                ctr[i, j].prim,
                a1face[i, j].fw,
                a1face[i+1, j].fw,
                a2face[i, j].fw,
                a2face[i, j+1].fw,
                KS.gas.γ,
                KS.ps.dx[i, j] * KS.ps.dy[i, j],
                sumRes,
                sumAvg,
                ctr[i, j].a,
                dt,
            )
        end
    end

    for i in eachindex(residual)
        residual[i] = sqrt(sumRes[i] * nx * ny) / (sumAvg[i] + 1.e-7)
    end

    KitBase.bc_extra!(ctr; dirc=:xr)
    KitBase.bc_extra!(ctr; dirc=:yr)
    KitBase.bc_mirror!(ctr; dirc=:yl)

    return nothing
end

@showprogress for iter in 1:50#nt
    evolve!(ks, ctr, a1face, a2face, dt)
    ib_force!(ks.ps, ctr, ib.xlp)
    updatef!(ks, ctr, a1face, a2face, dt, res)

    global t += dt
end

plot(ks, ctr)
