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
        ["maxwell", "fix", "slip", "slip"],
        0.5, # cfl
        10.0, # time
    )
    ps = CSpace2D(1.0, 6.0, 30, 0.0, π, 50, 0, 1)
    vs = VSpace2D(-8.0, 8.0, 48, -8.0, 8.0, 48)
    gas = Gas(5e-2, 2.0, 1.0, 1.0)
    
    prim0 = [1.0, 0.0, 0.0, 1.0]
    prim1 = [1.0, gas.Ma * sound_speed(1.0, gas.γ), 0.0, 1.0]
    fw = (args...) -> prim_conserve(prim1, gas.γ)
    ff = function(args...)
        prim = conserve_prim(fw(args...), gas.γ)
        h = maxwellian(vs.u, vs.v, prim)
        b = h .* gas.K / 2 / prim[end]
        return h, b
    end
    bc = function(x, y)
        if abs(x^2 + y^2 - 1) < 1e-3
            return prim0
        else
            return prim1
        end
    end
    ib = IB2F(fw, ff, bc)

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

    resL = zero(ctr[1].w)
    avgL = zero(ctr[1].w)
    resR = zero(ctr[1].w)
    avgR = zero(ctr[1].w)
    resU = zero(ctr[1].w)
    avgU = zero(ctr[1].w)
    resD = zero(ctr[1].w)
    avgD = zero(ctr[1].w)

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

ctr, a1face, a2face = init_fvm(ks, ks.pSpace)

t = 0.0
dt = timestep(ks, ctr, 0.0)
@showprogress for iter = 1:50#nt
    evolve!(ks, ctr, a1face, a2face, dt; bc = [:maxwell, :fix, :slip, :slip])
    update!(ks, ctr, a1face, a2face, dt, zeros(4))
    boundary!(ks, ctr, a1face, a2face, dt, Symbol(ks.set.collision))

    global t += dt
end

begin
    sol = zeros(axes(ps.x)..., 4)
    for i in axes(sol, 1), j in axes(sol, 2)
        sol[i, j, :] .= ctr[i, j].prim
        sol[i, j, end] = 1 / sol[i, j, end]
    end
    contourf(ps.x, ps.y, sol[:, :, 4], ratio=1)
end
