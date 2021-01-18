using ProgressMeter, Plots
import KitBase

cd(@__DIR__)
ks = KitBase.SolverSet("cavity.txt")

w = zeros(4, ks.pSpace.nx, ks.pSpace.ny)
prim = similar(w)
h = zeros(ks.vSpace.nu, ks.vSpace.nv, ks.pSpace.nx, ks.pSpace.ny)
b = similar(h)
for i = 1:ks.pSpace.nx, j = 1:ks.pSpace.ny
    prim[:, i, j] .= [1., 0., 0., 1.]
    w[:, i, j] .= KitBase.prim_conserve(prim[:, i, j], ks.gas.γ)
    h[:, :, i, j] .= KitBase.maxwellian(ks.vSpace.u, ks.vSpace.v, prim[:, i, j])
    b[:, :, i, j] .= h[:, :, i, j] .* ks.gas.K ./ (2.0 * prim[end, i, j])
end

fw1 = zeros(4, ks.pSpace.nx+1, ks.pSpace.ny)
fw2 = zeros(4, ks.pSpace.nx, ks.pSpace.ny+1)
fh1 = zeros(ks.vSpace.nu, ks.vSpace.nv, ks.pSpace.nx+1, ks.pSpace.ny)
fb1 = similar(fh1)
fh2 = zeros(ks.vSpace.nu, ks.vSpace.nv, ks.pSpace.nx, ks.pSpace.ny+1)
fb2 = similar(fh2)

res = zeros(4)
dt = 0.0015
nt = 100

@showprogress for iter = 1:10
    
    @inbounds for j = 1:ks.pSpace.ny
        for i = 2:ks.pSpace.nx
            fw = @view fw1[:, i, j]
            fh = @view fh1[:, :, i, j]
            fb = @view fb1[:, :, i, j]

            KitBase.flux_kfvs!(
                fw,
                fh,
                fb,
                h[:, :, i-1, j],
                b[:, :, i-1, j],
                h[:, :, i, j],
                b[:, :, i, j],
                ks.vSpace.u,
                ks.vSpace.v,
                ks.vSpace.weights,
                dt,
                ks.pSpace.dy[i, j],
            )
        end
    end
    
    vn = ks.vSpace.v
    vt = -ks.vSpace.u
    @inbounds for j = 2:ks.pSpace.ny
        for i = 1:ks.pSpace.nx
            fw = @view fw2[:, i, j]
            fh = @view fh2[:, :, i, j]
            fb = @view fb2[:, :, i, j]

            KitBase.flux_kfvs!(
                fw,
                fh,
                fb,
                h[:, :, i, j-1],
                b[:, :, i, j-1],
                h[:, :, i, j],
                b[:, :, i, j],
                vn,
                vt,
                ks.vSpace.weights,
                dt,
                ks.pSpace.dx[i, j],
            )
            fw .= KitBase.global_frame(fw, 0., 1.)
        end
    end
    
    #=
    for j = 1:ks.pSpace.ny
        KitBase.flux_boundary_maxwell!(
            a1face[1, j].fw,
            a1face[1, j].fh,
            a1face[1, j].fb,
            ks.ib.bcL,
            ctr[1, j].h,
            ctr[1, j].b,
            ks.vSpace.u,
            ks.vSpace.v,
            ks.vSpace.weights,
            ks.gas.K,
            dt,
            ctr[1, j].dy,
            1.,
        )

        KitBase.flux_boundary_maxwell!(
            a1face[ks.pSpace.nx+1, j].fw,
            a1face[ks.pSpace.nx+1, j].fh,
            a1face[ks.pSpace.nx+1, j].fb,
            ks.ib.bcR,
            ctr[ks.pSpace.nx, j].h,
            ctr[ks.pSpace.nx, j].b,
            ks.vSpace.u,
            ks.vSpace.v,
            ks.vSpace.weights,
            ks.gas.K,
            dt,
            ctr[ks.pSpace.nx, j].dy,
            -1.,
        )
    end
    
    
    for i = 1:ks.pSpace.nx
        KitBase.flux_boundary_maxwell!(
            a2face[i, 1].fw,
            a2face[i, 1].fh,
            a2face[i, 1].fb,
            ks.ib.bcD,
            ctr[i, 1].h,
            ctr[i, 1].b,
            vn,
            vt,
            ks.vSpace.weights,
            ks.gas.K,
            dt,
            ctr[i, 1].dx,
            1,
        )
        a2face[i, 1].fw .= KitBase.global_frame(a2face[i, 1].fw, 0., 1.)
        
        KitBase.flux_boundary_maxwell!(
            a2face[i, ks.pSpace.ny+1].fw,
            a2face[i, ks.pSpace.ny+1].fh,
            a2face[i, ks.pSpace.ny+1].fb,
            [1., 0.0, -1., 1.0],
            ctr[i, ks.pSpace.ny].h,
            ctr[i, ks.pSpace.ny].b,
            vn,
            vt,
            ks.vSpace.weights,
            ks.gas.K,
            dt,
            ctr[i, ks.pSpace.ny].dy,
            -1,
        )
        a2face[i, ks.pSpace.ny+1].fw .= KitBase.global_frame(
            a2face[i, ks.pSpace.ny+1].fw,
            0.,
            1.,
        )
    end
    =#

    #KitBase.update!(ks, ctr, a1face, a2face, dt, res; coll = :bgk, bc = :maxwell)

    for j in 1:ks.pSpace.ny, i in 1:ks.pSpace.nx
        @. w[:, i, j] += (fw1[:, i, j] - fw1[:, i+1, j] + fw2[:, i, j] - fw2[:, i, j+1]) / ks.pSpace.dx[1, 1] / ks.pSpace.dy[1, 1]
        prim[:, i, j] .= KitBase.conserve_prim(w[:, i, j], ks.gas.γ)

        MH = KitBase.maxwellian(ks.vSpace.u, ks.vSpace.v, prim[:, i, j])
        MB = MH .* ks.gas.K ./ (2.0 * prim[end, i, j])
        τ = KitBase.vhs_collision_time(prim[:, i, j], ks.gas.μᵣ, ks.gas.ω)

        for q in axes(ks.vSpace.v, 2), p in axes(ks.vSpace.u, 1)
            h[p, q, i, j] = (h[p, q, i, j] + (fh1[p, q, i, j] - fh1[p, q, i+1, j] + 
            fh2[p, q, i, j] - fh2[p, q, i, j+1]) / ks.pSpace.dx[1, 1] / ks.pSpace.dy[1, 1] + dt / τ * MH[p, q]) / (1.0 + dt / τ)
            
        end
    end

end

contourf(ks.pSpace.x[1:ks.pSpace.nx, 1], ks.pSpace.y[1, 1:ks.pSpace.ny], prim[4, :, :]')