using LinearAlgebra
cd(@__DIR__)
D = KitBase.read_dict("t1_2f.txt")
set = KitBase.set_setup(D)
ps = KitBase.set_geometry(D)
vs = KitBase.set_velocity(D)
gas = KitBase.set_property(D)
begin
    primL = [1.0, KitBase.sound_speed(1.0, gas.γ) * gas.Ma, 0.0, 1.0]
    wL = KitBase.prim_conserve(primL, gas.γ)
    hL = KitBase.maxwellian(vs.u, vs.v, primL)
    bL = @. hL * gas.K / 2 / primL[end]
    ib = KitBase.IB2F(wL, primL, hL, bL, primL, wL, primL, hL, bL, primL)
end
ks = KitBase.SolverSet(set, ps, vs, gas, ib)
ctr, face = KitBase.init_fvm(ks, ks.pSpace)
dt = KitBase.timestep(ks, ctr, 0.0)
nt = ks.set.maxTime ÷ dt |> Int
for iter = 1:nt
    KitBase.reconstruct!(ks, ctr)
    KitBase.evolve!(ks, ctr, face, dt)
    #=
    @inbounds for i in eachindex(face)
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
    =#
    res = zeros(4)
    KitBase.update!(ks, ctr, face, dt, res)
    #=
    sumres = zeros(4)
    sumavg = zeros(4)
    @inbounds for i in eachindex(ctr)
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
    =#
end
KitBase.write_vtk(ks, ctr)

D = KitBase.read_dict("t1_1f.txt")
set = KitBase.set_setup(D)
begin
    primL = [1.0, KitBase.sound_speed(1.0, gas.γ) * gas.Ma, 0.0, 1.0]
    wL = KitBase.prim_conserve(primL, gas.γ)
    hL = KitBase.maxwellian(vs.u, vs.v, primL)
    ib = KitBase.IB1F(wL, primL, hL, primL, wL, primL, hL, primL)
end
ks = KitBase.SolverSet(set, ps, vs, gas, ib, @__DIR__)
ctr, face = KitBase.init_fvm(ks, ks.pSpace)
KitBase.reconstruct!(ks, ctr)
KitBase.evolve!(ks, ctr, face, dt)
res = zeros(4)
KitBase.update!(ks, ctr, face, dt, res)

D = KitBase.read_dict("t1_0f.txt")
set = KitBase.set_setup(D)
ks = KitBase.SolverSet(set, ps, vs, gas, ib, @__DIR__)
ctr, face = KitBase.init_fvm(ks, ks.pSpace)
KitBase.reconstruct!(ks, ctr)
KitBase.evolve!(ks, ctr, face, dt)
res = zeros(4)
KitBase.update!(ks, ctr, face, dt, res)
