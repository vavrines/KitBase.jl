using KitBase, Plots
pyplot()

begin
    set = Setup(
        "gas",
        "cylinder",
        "2d2f",
        "kfvs",
        "bgk",
        1, # species
        2, # order of accuracy
        "vanleer", # limiter
        "maxwell",
        0.5, # cfl
        10.0, # time
    )

    ps = CSpace2D(1.0, 6.0, 30, 0.0, π, 50, 0, 1)
    vs = VSpace2D(-8.0, 8.0, 48, -8.0, 8.0, 48)

    gas = Gas(
        5e-2,
        2.0, # Mach
        1.0,
        1.0, # K
        5 / 3,
        0.81,
        1.0,
        0.5,
    )
    
    prim0 = [1.0, 0.0, 0.0, 1.0]
    w0 = prim_conserve(prim0, gas.γ)
    h0 = maxwellian(vs.u, vs.v, prim0)
    b0 = h0 * gas.K / 2 / prim0[end]
    ib = IB2F(w0, prim0, h0, b0, prim0, w0, prim0, h0, b0, prim0)

    ks = SolverSet(set, ps, vs, gas, ib)
end

ctr, a1face, a2face = init_fvm(ks, ks.pSpace; structarray = true)

dt = timestep(ks, ctr, 0.0)

evolve!(ks, ctr, a1face, a2face, dt; bc = :fix)

update!(ks, ctr, a1face, a2face, dt, zeros(4))


begin
    sol = zeros(axes(ps.x)..., 4)
    for i in axes(sol, 1), j in axes(sol, 2)
        sol[i, j, :] .= ctr[i, j].prim
    end
    contourf(ps.x, ps.y, sol[:, :, 4], ratio=1)
end
