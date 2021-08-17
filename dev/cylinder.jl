using KitBase, Plots
using KitBase.ProgressMeter: @showprogress
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
    
    prim0 = [1.0, gas.Ma * sound_speed(1.0, gas.γ), 0.0, 1.0]
    w0 = prim_conserve(prim0, gas.γ)
    h0 = maxwellian(vs.u, vs.v, prim0)
    b0 = h0 * gas.K / 2 / prim0[end]
    ib = IB2F(w0, prim0, h0, b0, prim0, w0, prim0, h0, b0, prim0)

    ks = SolverSet(set, ps, vs, gas, ib)
end

function boundary!(
    KS::X,
    ctr::Y,
    a1face::Z,
    a2face::Z,
    dt,
    coll::Symbol,
) where {
    X<:AbstractSolverSet,
    Y<:AbstractArray{ControlVolume2D2F,2},
    Z<:AbstractArray{Interface2D2F,2},
}

    nx, ny, dx, dy = begin
        if KS.ps isa CSpace2D
            KS.ps.nr, KS.ps.nθ, KS.ps.dr, KS.ps.darc
        else
            KS.ps.nx, KS.ps.ny, KS.ps.dx, KS.ps.dy
        end
    end

    resL = zero(KS.ib.wL)
    avgL = zero(KS.ib.wL)
    resR = zero(KS.ib.wL)
    avgR = zero(KS.ib.wL)
    resU = zero(KS.ib.wL)
    avgU = zero(KS.ib.wL)
    resD = zero(KS.ib.wL)
    avgD = zero(KS.ib.wL)

    @inbounds for j = 1:ny
        KitBase.step!(
            ctr[1, j].w,
            ctr[1, j].prim,
            ctr[1, j].h,
            ctr[1, j].b,
            a1face[1, j].fw,
            a1face[1, j].fh,
            a1face[1, j].fb,
            a1face[2, j].fw,
            a1face[2, j].fh,
            a1face[2, j].fb,
            a2face[1, j].fw,
            a2face[1, j].fh,
            a2face[1, j].fb,
            a2face[1, j+1].fw,
            a2face[1, j+1].fh,
            a2face[1, j+1].fb,
            KS.vSpace.u,
            KS.vSpace.v,
            KS.vSpace.weights,
            KS.gas.K,
            KS.gas.γ,
            KS.gas.μᵣ,
            KS.gas.ω,
            KS.gas.Pr,
            dx[1, j] * dy[1, j],
            dt,
            resL,
            avgL,
            coll,
        )
    end

    @inbounds for i = 2:nx-1 # skip overlap
        KitBase.step!(
            ctr[i, 1].w,
            ctr[i, 1].prim,
            ctr[i, 1].h,
            ctr[i, 1].b,
            a1face[i, 1].fw,
            a1face[i, 1].fh,
            a1face[i, 1].fb,
            a1face[i+1, 1].fw,
            a1face[i+1, 1].fh,
            a1face[i+1, 1].fb,
            a2face[i, 1].fw,
            a2face[i, 1].fh,
            a2face[i, 1].fb,
            a2face[i, 2].fw,
            a2face[i, 2].fh,
            a2face[i, 2].fb,
            KS.vSpace.u,
            KS.vSpace.v,
            KS.vSpace.weights,
            KS.gas.K,
            KS.gas.γ,
            KS.gas.μᵣ,
            KS.gas.ω,
            KS.gas.Pr,
            dx[i, 1] * dy[i, 1],
            dt,
            resD,
            avgD,
            coll,
        )

        KitBase.step!(
            ctr[i, ny].w,
            ctr[i, ny].prim,
            ctr[i, ny].h,
            ctr[i, ny].b,
            a1face[i, ny].fw,
            a1face[i, ny].fh,
            a1face[i, ny].fb,
            a1face[i+1, ny].fw,
            a1face[i+1, ny].fh,
            a1face[i+1, ny].fb,
            a2face[i, ny].fw,
            a2face[i, ny].fh,
            a2face[i, ny].fb,
            a2face[i, ny+1].fw,
            a2face[i, ny+1].fh,
            a2face[i, ny+1].fb,
            KS.vSpace.u,
            KS.vSpace.v,
            KS.vSpace.weights,
            KS.gas.K,
            KS.gas.γ,
            KS.gas.μᵣ,
            KS.gas.ω,
            KS.gas.Pr,
            dx[i, ny] * dy[i, ny],
            dt,
            resU,
            avgU,
            coll,
        )
    end

    for i = 1:nx
        ctr[i, 0].w .= ctr[i, 1].w
        ctr[i, 0].w[3] = -ctr[i, 0].w[3]
        ctr[i, 0].prim .= ctr[i, 1].prim
        ctr[i, 0].prim[3] = -ctr[i, 0].prim[3]
        ctr[i, ny+1].w .= ctr[i, ny].w
        ctr[i, ny+1].w[3] = -ctr[i, ny+1].w[3]
        ctr[i, ny+1].prim .= ctr[i, ny].prim
        ctr[i, ny+1].prim[3] = -ctr[i, ny+1].prim[3]

        for p = 1:KS.vs.nu, q = 1:KS.vs.nv
            ctr[i, 0].h[p, q] = ctr[i, 1].h[p, KS.vs.nv+1-q]
            ctr[i, 0].b[p, q] = ctr[i, 1].b[p, KS.vs.nv+1-q]
            ctr[i, ny+1].h[p, q] = ctr[i, ny].h[p, KS.vs.nv+1-q]
            ctr[i, ny+1].b[p, q] = ctr[i, ny].b[p, KS.vs.nv+1-q]
        end
    end
end

ctr, a1face, a2face = init_fvm(ks, ks.pSpace; structarray = true)

t = 0.0

dt = timestep(ks, ctr, 0.0)


@showprogress for iter = 1:50#nt
    evolve!(ks, ctr, a1face, a2face, dt; bc = :maxwell)
    update!(ks, ctr, a1face, a2face, dt, zeros(4))
    boundary!(ks, ctr, a1face, a2face, dt, Symbol(ks.set.collision))
end

begin
    sol = zeros(axes(ps.x)..., 4)
    for i in axes(sol, 1), j in axes(sol, 2)
        sol[i, j, :] .= ctr[i, j].prim
    end
    contourf(ps.x, ps.y, sol[:, :, 2], ratio=1)
end
