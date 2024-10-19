using KitBase, Plots
using KitBase.ProgressMeter: @showprogress
pyplot()

begin
    set = Setup(
        "gas",
        "cylinder",
        "2d2f2v",
        "kfvs",
        "bgk",
        1, # species
        2, # order of accuracy
        "vanleer", # limiter
        ["maxwell", "extra", "mirrror", "mirror"],
        0.3, # cfl
        10.0, # time
    )
    ps = CSpace2D(1.0, 6.0, 20, 0.0, π, 40, 1, 1)
    vs = VSpace2D(-10.0, 10.0, 36, -10.0, 10.0, 36)
    gas = Gas(; Kn=5e-2, Ma=5.0, Pr=1.0, K=1.0)

    prim0 = [1.0, 0.0, 0.0, 1.0]
    prim1 = [1.0, gas.Ma * sound_speed(1.0, gas.γ), 0.0, 1.0]
    fw = (args...) -> prim_conserve(prim1, gas.γ)
    ff = function (args...)
        prim = conserve_prim(fw(args...), gas.γ)
        h = maxwellian(vs.u, vs.v, prim)
        b = h .* gas.K / 2 / prim[end]
        return h, b
    end
    bc = function (x, y)
        if abs(x^2 + y^2 - 1) < 1e-3
            return prim0
        else
            return prim1
        end
    end
    ib = IB2F(fw, ff, bc)

    ks = SolverSet(set, ps, vs, gas, ib)
end

ctr, a1face, a2face = KitBase.init_fvm(ks, ks.ps)

t = 0.0
dt = timestep(ks, ctr, 0.0)
nt = ks.set.maxTime ÷ dt |> Int
@showprogress for iter in 1:500#nt
    evolve!(ks, ctr, a1face, a2face, dt; bc=[:maxwell, :extra, :mirror, :mirror])
    update!(ks, ctr, a1face, a2face, dt, zeros(4); bc=[:maxwell, :extra, :mirror, :mirror])

    for j in ks.ps.nθ÷2+1:ks.ps.nθ
        ctr[ks.ps.nr+1, j].w .= ks.ib.fw(6, 0)
        ctr[ks.ps.nr+1, j].prim .= conserve_prim(ctr[ks.ps.nr+1, j].w, ks.gas.γ)
        ctr[ks.ps.nr+1, j].sw .= 0.0
        ctr[ks.ps.nr+1, j].h .= maxwellian(ks.vs.u, ks.vs.v, ctr[ks.ps.nr+1, j].prim)
        ctr[ks.ps.nr+1, j].b .=
            ctr[ks.ps.nr+1, j].h .* ks.gas.K ./ 2 ./ ctr[ks.ps.nr+1, j].prim[end]
    end

    global t += dt
end

begin
    sol = zeros(ks.ps.nr, ks.ps.nθ, 4)
    for i in axes(sol, 1), j in axes(sol, 2)
        sol[i, j, :] .= ctr[i, j].prim
        sol[i, j, end] = 1 / sol[i, j, end]
    end
    contourf(
        ps.x[1:ks.ps.nr, 1:ks.ps.nθ],
        ps.y[1:ks.ps.nr, 1:ks.ps.nθ],
        sol[:, :, 1];
        ratio=1,
    )
end
