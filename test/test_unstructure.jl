using LinearAlgebra
cd(@__DIR__)
D = KB.read_dict("t1_2f.txt")
set = KB.set_setup(D)
ps = KB.set_geometry(D)
vs = KB.set_velocity(D)
gas = KB.set_property(D)
begin
    primL = [1.0, KB.sound_speed(1.0, gas.γ) * gas.Ma, 0.0, 1.0]
    wL = KB.prim_conserve(primL, gas.γ)
    hL = KB.maxwellian(vs.u, vs.v, primL)
    bL = @. hL * gas.K / 2 / primL[end]

    fw = function (x, y, p)
        return wL
    end
    ff = function (x, y, p)
        h = maxwellian(vs.u, vs.v, primL)
        b = h * gas.K / 2 / primL[end]
        return h, b
    end
    bc = function (x, y, p)
        return primL
    end
    ib = KB.IB2F(fw, ff, bc, NamedTuple())
end
ks = KB.SolverSet(set, ps, vs, gas, ib)
ctr, face = KB.init_fvm(ks, ks.ps)
dt = KB.timestep(ks, ctr, 0.0)
nt = ks.set.maxTime ÷ dt |> Int
for iter in 1:nt
    KB.reconstruct!(ks, ctr)
    KB.evolve!(ks, ctr, face, dt)
    #=
    @inbounds for i in eachindex(face)
        vn = ks.vs.u .* face[i].n[1] .+ ks.vs.v .* face[i].n[2]
        vt = ks.vs.v .* face[i].n[1] .- ks.vs.u .* face[i].n[2]

        if !(-1 in ps.faceCells[i, :])
            KB.flux_kfvs!(
                face[i].fw,
                face[i].fh,
                face[i].fb,
                ctr[ps.faceCells[i, 1]].h,
                ctr[ps.faceCells[i, 1]].b,
                ctr[ps.faceCells[i, 2]].h,
                ctr[ps.faceCells[i, 2]].b,
                vn,
                vt,
                ks.vs.weights,
                dt,
                face[i].len,
            )
            face[i].fw .= KB.global_frame(face[i].fw, face[i].n[1], face[i].n[2])
        else
            idx = ifelse(ps.faceCells[i, 1] != -1, 1, 2)

            if ps.cellType[ps.faceCells[i, idx]] == 2
                bc = KB.local_frame(ks.ib.primR, face[i].n[1], face[i].n[2])

                KB.flux_boundary_maxwell!(
                    face[i].fw,
                    face[i].fh,
                    face[i].fb,
                    bc,
                    ctr[ps.faceCells[i, idx]].h,
                    ctr[ps.faceCells[i, idx]].b,
                    vn,
                    vt,
                    ks.vs.weights,
                    ks.gas.K,
                    dt,
                    face[i].len,
                )

                face[i].fw .= KB.global_frame(face[i].fw, face[i].n[1], face[i].n[2])
            end
        end
    end
    =#
    res = zeros(4)
    KB.update!(ks, ctr, face, dt, res)
    #=
    sumres = zeros(4)
    sumavg = zeros(4)
    @inbounds for i in eachindex(ctr)
        if ps.cellType[i] in (0, 2)
            dirc = [sign(dot(ctr[i].n[j], face[ps.cellFaces[i, j]].n)) for j = 1:3]

            KB.step!(
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
                ks.vs.u,
                ks.vs.v,
                ks.vs.weights,
                ks.gas.K,
                ks.gas.γ,
                ks.gas.μᵣ,
                ks.gas.ω,
                ks.gas.Pr,
                ks.ps.cellArea[i],
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
            ctr[i].prim .= KB.conserve_prim(ctr[i].w, ks.gas.γ)
        end
    end
    =#
end
KB.write_vtk(ks, ctr)

D = KB.read_dict("t1_1f.txt")
set = KB.set_setup(D)
begin
    primL = [1.0, KB.sound_speed(1.0, gas.γ) * gas.Ma, 0.0, 1.0]
    wL = KB.prim_conserve(primL, gas.γ)
    hL = KB.maxwellian(vs.u, vs.v, primL)
    fw = function (x, y, p)
        return wL
    end
    ff = function (x, y, p)
        h = maxwellian(vs.u, vs.v, primL)
        return h
    end
    bc = function (x, y, p)
        return primL
    end
    ib = KB.IB1F(fw, ff, bc, NamedTuple())
end
ks = KB.SolverSet(set, ps, vs, gas, ib, @__DIR__)
ctr, face = KB.init_fvm(ks, ks.ps)
KB.reconstruct!(ks, ctr)
KB.evolve!(ks, ctr, face, dt)
res = zeros(4)
KB.update!(ks, ctr, face, dt, res)

D = KB.read_dict("t1_0f.txt")
set = KB.set_setup(D)
ks = KB.SolverSet(set, ps, vs, gas, ib, @__DIR__)
ctr, face = KB.init_fvm(ks, ks.ps)
KB.reconstruct!(ks, ctr)
KB.evolve!(ks, ctr, face, dt)
res = zeros(4)
KB.update!(ks, ctr, face, dt, res)
