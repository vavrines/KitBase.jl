using KitBase

cd(@__DIR__)
ks, ctr, a1face, a2face, simTime = KitBase.initialize("cavity.txt")
simTime = KitBase.solve!(ks, ctr, a1face, a2face, simTime)

# equivalent low-level procedures
using KitBase.ProgressMeter: @showprogress
simTime = 0.0
res = zeros(4)
dt = KitBase.timestep(ks, ctr, simTime)
nt = floor(ks.set.maxTime / dt) |> Int
@showprogress for iter = 1:nt
    KitBase.reconstruct!(ks, ctr)
    KitBase.evolve!(ks, ctr, a1face, a2face, dt; mode = :kfvs, bc = :maxwell)
    KitBase.update!(ks, ctr, a1face, a2face, dt, res; coll = :bgk, bc = :maxwell)
end

# lower-level backend 
@showprogress for iter = 1:nt
    @inbounds Threads.@threads for j = 1:ks.ps.ny
        for i = 2:ks.ps.nx
            KitBase.flux_kfvs!(
                a1face[i, j].fw,
                a1face[i, j].fh,
                a1face[i, j].fb,
                ctr[i-1, j].h,
                ctr[i-1, j].b,
                ctr[i, j].h,
                ctr[i, j].b,
                ks.vs.u,
                ks.vs.v,
                ks.vs.weights,
                dt,
                ks.ps.dy[i, j],
            )
        end
    end

    # vertical flux
    vn = ks.vs.v
    vt = -ks.vs.u
    @inbounds Threads.@threads for j = 2:ks.ps.ny
        for i = 1:ks.ps.nx
            KitBase.flux_kfvs!(
                a2face[i, j].fw,
                a2face[i, j].fh,
                a2face[i, j].fb,
                ctr[i, j-1].h,
                ctr[i, j-1].b,
                ctr[i, j].h,
                ctr[i, j].b,
                vn,
                vt,
                ks.vs.weights,
                dt,
                ks.ps.dy[i, j],
            )
            a2face[i, j].fw .= KitBase.global_frame(a2face[i, j].fw, 0.0, 1.0)
        end
    end

    # boundary flux
    @inbounds Threads.@threads for j = 1:ks.ps.ny
        KitBase.flux_boundary_maxwell!(
            a1face[1, j].fw,
            a1face[1, j].fh,
            a1face[1, j].fb,
            [1.0, 0.0, 0.0, 1.0],
            ctr[1, j].h,
            ctr[1, j].b,
            ks.vs.u,
            ks.vs.v,
            ks.vs.weights,
            ks.gas.K,
            dt,
            ks.ps.dy[1, j],
            1.0,
        )

        KitBase.flux_boundary_maxwell!(
            a1face[ks.ps.nx+1, j].fw,
            a1face[ks.ps.nx+1, j].fh,
            a1face[ks.ps.nx+1, j].fb,
            [1.0, 0.0, 0.0, 1.0],
            ctr[ks.ps.nx, j].h,
            ctr[ks.ps.nx, j].b,
            ks.vs.u,
            ks.vs.v,
            ks.vs.weights,
            ks.gas.K,
            dt,
            ks.ps.dy[ks.ps.nx, j],
            -1.0,
        )
    end

    @inbounds Threads.@threads for i = 1:ks.ps.nx
        KitBase.flux_boundary_maxwell!(
            a2face[i, 1].fw,
            a2face[i, 1].fh,
            a2face[i, 1].fb,
            [1.0, 0.0, 0.0, 1.0],
            ctr[i, 1].h,
            ctr[i, 1].b,
            vn,
            vt,
            ks.vs.weights,
            ks.gas.K,
            dt,
            ks.ps.dx[i, 1],
            1,
        )
        a2face[i, 1].fw .= KitBase.global_frame(a2face[i, 1].fw, 0.0, 1.0)

        KitBase.flux_boundary_maxwell!(
            a2face[i, ks.ps.ny+1].fw,
            a2face[i, ks.ps.ny+1].fh,
            a2face[i, ks.ps.ny+1].fb,
            [1.0, 0.0, -0.15, 1.0],
            ctr[i, ks.ps.ny].h,
            ctr[i, ks.ps.ny].b,
            vn,
            vt,
            ks.vs.weights,
            ks.gas.K,
            dt,
            ks.ps.dy[i, ks.ps.ny],
            -1,
        )
        a2face[i, ks.ps.ny+1].fw .= KitBase.global_frame(a2face[i, ks.ps.ny+1].fw, 0.0, 1.0)
    end

    # update
    @inbounds Threads.@threads for j = 1:ks.ps.ny
        for i = 1:ks.ps.nx
            KitBase.step!(
                ctr[i, j].w,
                ctr[i, j].prim,
                ctr[i, j].h,
                ctr[i, j].b,
                a1face[i, j].fw,
                a1face[i, j].fh,
                a1face[i, j].fb,
                a1face[i+1, j].fw,
                a1face[i+1, j].fh,
                a1face[i+1, j].fb,
                a2face[i, j].fw,
                a2face[i, j].fh,
                a2face[i, j].fb,
                a2face[i, j+1].fw,
                a2face[i, j+1].fh,
                a2face[i, j+1].fb,
                ks.vs.u,
                ks.vs.v,
                ks.vs.weights,
                ks.gas.K,
                ks.gas.γ,
                ks.gas.μᵣ,
                ks.gas.ω,
                ks.gas.Pr,
                ks.ps.dx[i, j] * ks.ps.dy[i, j],
                dt,
                zeros(4),
                zeros(4),
                :bgk,
            )
        end
    end
end

# visulization
using Plots
plot(ks, ctr)
#savefig("cavity.png")

# low-level backend
begin
    using Plots
    sol = zeros(4, ks.ps.nx, ks.ps.ny)
    for i in axes(sol, 2)
        for j in axes(sol, 3)
            sol[1:3, i, j] .= ctr[i, j].prim[1:3]
            sol[4, i, j] = 1.0 / ctr[i, j].prim[4]
        end
    end
    contourf(ks.ps.x[1:ks.ps.nx, 1], ks.ps.y[1, 1:ks.ps.ny], sol[4, :, :]')
end
