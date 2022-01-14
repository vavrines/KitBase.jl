using LinearAlgebra, WriteVTK, ProgressMeter, JLD2, PyCall
import KitBase

cd(@__DIR__)
meshio = pyimport("meshio")
m = meshio.read("../assets/mesh/cylinder_circle.msh")

D = KitBase.read_dict("cylinder.txt")
set = KitBase.set_setup(D)

ps = KitBase.set_geometry(D)
for i in eachindex(ps.faceType)
    i1 = ps.facePoints[i, 1]
    i2 = ps.facePoints[i, 2]

    if i1 in [ps.cells.index[1]; ps.cells.index[2]] &&
       i2 in [ps.cells.index[1]; ps.cells.index[2]]
        ps.faceType[i] = 2
        c1 = ps.faceCells[i, 1]
        c2 = ps.faceCells[i, 2]

        if c1 != -1
            ps.cellType[c1] = 2
        elseif c2 != -1
            ps.cellType[c2] = 2
        else
            throw("index error")
        end
    end
end
for i in eachindex(ps.cellType)
    if ps.cellType[i] == 1 && ps.cellCenter[i, 1] > 0.0
        ps.cellType[i] = 3
    end
end

vs = KitBase.set_velocity(D)
gas = KitBase.set_property(D)

begin
    primL = [1.0, KitBase.sound_speed(1.0, gas.γ) * gas.Ma, 0.0, 1.0]
    wL = KitBase.prim_conserve(primL, gas.γ)
    hL = KitBase.maxwellian(vs.u, vs.v, primL)
    bL = @. hL * gas.K / 2 / primL[end]
    primR = [1.0, 0.0, 0.0, 1.0]
    wR = KitBase.prim_conserve(primR, gas.γ)
    hR = KitBase.maxwellian(vs.u, vs.v, primR)
    bR = @. hR * gas.K / 2 / primR[end]
    ib = KitBase.IB2F(wL, primL, hL, bL, primL, wL, primL, hL, bL, primR)
end

ks = KitBase.SolverSet(set, ps, vs, gas, ib, @__DIR__)

ctr, face = KitBase.init_fvm(ks, ks.pSpace)

dt = KitBase.timestep(ks, ctr, 0.0)
nt = ks.set.maxTime ÷ dt |> Int
#=
@showprogress for iter = 1:2000#nt
    @inbounds Threads.@threads for i in eachindex(face)
        vn = ks.vSpace.u .* face[i].n[1] .+ ks.vSpace.v .* face[i].n[2]
        vt = ks.vSpace.v .* face[i].n[1] .- ks.vSpace.u .* face[i].n[2]

        if !(-1 in ps.faceCells[i, :])
            KitBase.flux_kfvs!(
                face[i].fw,
                face[i].fh,
                face[i].fb,
                ctr[ps.faceCells[i, 1]].h,
                ctr[ps.faceCells[i, 1]].b,
                ctr[ps.faceCells[i, 2]].h,
                ctr[ps.faceCells[i, 2]].b,
                vn,
                vt,
                ks.vSpace.weights,
                dt,
                face[i].len,
            )
            face[i].fw .= KitBase.global_frame(face[i].fw, face[i].n[1], face[i].n[2])
        else
            idx = ifelse(ps.faceCells[i, 1] != -1, 1, 2)

            if ps.cellType[ps.faceCells[i, idx]] == 2
                bc = KitBase.local_frame(ks.ib.primR, face[i].n[1], face[i].n[2])

                KitBase.flux_boundary_maxwell!(
                    face[i].fw,
                    face[i].fh,
                    face[i].fb,
                    bc,
                    ctr[ps.faceCells[i, idx]].h,
                    ctr[ps.faceCells[i, idx]].b,
                    vn,
                    vt,
                    ks.vSpace.weights,
                    ks.gas.K,
                    dt,
                    face[i].len,
                )

                face[i].fw .= KitBase.global_frame(face[i].fw, face[i].n[1], face[i].n[2])
            end
        end
    end

    @inbounds Threads.@threads for i in eachindex(ctr)
        if ps.cellType[i] in (0, 2)
            MH = KitBase.maxwellian(vs.u, vs.v, ctr[i].prim)
            MB = MH .* ks.gas.K ./ (2.0 * ctr[i].prim[end])
            τ = KitBase.vhs_collision_time(ctr[i].prim, ks.gas.μᵣ, ks.gas.ω)

            for j in 1:3
                dirc = sign(dot(ctr[i].n[j], face[ps.cellFaces[i, j]].n))
                @. ctr[i].w -= dirc * face[ps.cellFaces[i, j]].fw / ps.cellArea[i]
                @. ctr[i].h -= dirc * face[ps.cellFaces[i, j]].fh / ps.cellArea[i]
                @. ctr[i].b -= dirc * face[ps.cellFaces[i, j]].fb / ps.cellArea[i]
            end
            ctr[i].prim .= KitBase.conserve_prim(ctr[i].w, ks.gas.γ)

            @. ctr[i].h += (MH - ctr[i].h) * dt / τ
            @. ctr[i].b += (MB - ctr[i].b) * dt / τ
            #=for q in axes(vs.v, 2), p in axes(vs.u, 1)
                ctr[i].h[p, q] = (ctr[i].h[p, q] - (face[ps.cellFaces[i, 1]].fh[p, q] + face[ps.cellFaces[i, 2]].fh[p, q] + face[ps.cellFaces[i, 3]].fh[p, q]) / ps.cellArea[i] + 
                    dt / τ * MH[p, q]) / (1.0 + dt / τ)
                ctr[i].b[p, q] = (ctr[i].b[p, q] - (face[ps.cellFaces[i, 1]].fb[p, q] + face[ps.cellFaces[i, 2]].fb[p, q] + face[ps.cellFaces[i, 3]].fb[p, q]) / ps.cellArea[i] + 
                    dt / τ * MB[p, q]) / (1.0 + dt / τ)
            end=#
        end
    end

    for i in eachindex(ps.cellType)
        if ps.cellType[i] == 3
            ids = ps.cellNeighbors[i, :]
            deleteat!(ids, findall(x->x==-1, ids))
            id1, id2 = ids
            ctr[i].w .= 0.5 .* (ctr[id1].w .+ ctr[id2].w)
            ctr[i].h .= 0.5 .* (ctr[id1].h .+ ctr[id2].h)
            ctr[i].b .= 0.5 .* (ctr[id1].b .+ ctr[id2].b)
            ctr[i].prim .= KitBase.conserve_prim(ctr[i].w, ks.gas.γ)
        end
    end
end
=#

@showprogress for iter = 1:nt
    @inbounds Threads.@threads for i in eachindex(face)
        vn = ks.vSpace.u .* face[i].n[1] .+ ks.vSpace.v .* face[i].n[2]
        vt = ks.vSpace.v .* face[i].n[1] .- ks.vSpace.u .* face[i].n[2]

        if !(-1 in ps.faceCells[i, :])
            KitBase.flux_kfvs!(
                face[i].fw,
                face[i].fh,
                face[i].fb,
                ctr[ps.faceCells[i, 1]].h,
                ctr[ps.faceCells[i, 1]].b,
                ctr[ps.faceCells[i, 2]].h,
                ctr[ps.faceCells[i, 2]].b,
                vn,
                vt,
                ks.vSpace.weights,
                dt,
                face[i].len,
            )
            face[i].fw .= KitBase.global_frame(face[i].fw, face[i].n[1], face[i].n[2])
        else
            idx = ifelse(ps.faceCells[i, 1] != -1, 1, 2)

            if ps.cellType[ps.faceCells[i, idx]] == 2
                bc = KitBase.local_frame(ks.ib.primR, face[i].n[1], face[i].n[2])

                KitBase.flux_boundary_maxwell!(
                    face[i].fw,
                    face[i].fh,
                    face[i].fb,
                    bc,
                    ctr[ps.faceCells[i, idx]].h,
                    ctr[ps.faceCells[i, idx]].b,
                    vn,
                    vt,
                    ks.vSpace.weights,
                    ks.gas.K,
                    dt,
                    face[i].len,
                )

                face[i].fw .= KitBase.global_frame(face[i].fw, face[i].n[1], face[i].n[2])
            end
        end
    end

    sumres = zeros(4)
    sumavg = zeros(4)
    @inbounds Threads.@threads for i in eachindex(ctr)
        if ps.cellType[i] in (0, 2)
            dirc = [sign(dot(ctr[i].n[j], face[ps.cellFaces[i, j]].n)) for j = 1:3]

            KitBase.step!(
                ctr[i].w,
                ctr[i].prim,
                ctr[i].h,
                ctr[i].b,
                face[ps.cellFaces[i, 1]].fw,
                face[ps.cellFaces[i, 1]].fh,
                face[ps.cellFaces[i, 1]].fb,
                face[ps.cellFaces[i, 2]].fw,
                face[ps.cellFaces[i, 2]].fh,
                face[ps.cellFaces[i, 2]].fb,
                face[ps.cellFaces[i, 3]].fw,
                face[ps.cellFaces[i, 3]].fh,
                face[ps.cellFaces[i, 3]].fb,
                ks.vSpace.u,
                ks.vSpace.v,
                ks.vSpace.weights,
                ks.gas.K,
                ks.gas.γ,
                ks.gas.μᵣ,
                ks.gas.ω,
                ks.gas.Pr,
                ks.pSpace.cellArea[i],
                dirc,
                dt,
                sumres,
                sumavg,
                :bgk,
            )
        end
    end

    for i in eachindex(ps.cellType)
        if ps.cellType[i] == 3
            ids = ps.cellNeighbors[i, :]
            deleteat!(ids, findall(x -> x == -1, ids))
            id1, id2 = ids
            ctr[i].w .= 0.5 .* (ctr[id1].w .+ ctr[id2].w)
            ctr[i].h .= 0.5 .* (ctr[id1].h .+ ctr[id2].h)
            ctr[i].b .= 0.5 .* (ctr[id1].b .+ ctr[id2].b)
            ctr[i].prim .= KitBase.conserve_prim(ctr[i].w, ks.gas.γ)
        end
    end
end

KitBase.write_vtk(ks, ctr)
