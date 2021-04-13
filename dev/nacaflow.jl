using LinearAlgebra, ProgressMeter, JLD2, WriteVTK, PyCall
import KitBase
#=function KitBase.write_vtk(ps::KitBase.AbstractPhysicalSpace, ctr)
    mcells = [MeshCell(VTKCellTypes.VTK_TRIANGLE, ps.cells[i, :]) for i in axes(ps.cells, 1)]
    vtkfile = vtk_grid("sol", permutedims(ps.points), mcells)
    
    cdata = zeros(length(ctr), 4)
    for i in eachindex(ctr)
        cdata[i, :] .= ctr[i].prim
    end

    vtkfile["Density", VTKCellData()] = cdata[:, 1]
    vtkfile["U", VTKCellData()] = cdata[:, 2]
    vtkfile["V", VTKCellData()] = cdata[:, 3]
    vtkfile["Temperature", VTKCellData()] = 1.0 ./ cdata[:, 4]

    outfiles = vtk_save(vtkfile)

    return nothing
end=#

cd(@__DIR__)
#@load "ctr.jld2" ctr
D = KitBase.read_dict("nacaflow.txt")
set = KitBase.set_setup(D)

ps = KitBase.set_geometry(D)
for i in eachindex(ps.cellType)
    if ps.cellType[i] == 1 && -1.0 < ps.cellCenter[i, 1] < 1.0
        ps.cellType[i] = 2

        for j = 1:3
            if ps.edgeType[ps.cellEdges[i, j]] == 1
                ps.edgeType[ps.cellEdges[i, j]] = 2
            end
        end
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
    ib = KitBase.IB2F(
        wL,
        primL,
        hL,
        bL,
        primL,
        wL,
        primL,
        hL,
        bL,
        primR,
    )
end

ks = KitBase.SolverSet(set, ps, vs, gas, ib, @__DIR__)

ctr, face = KitBase.init_fvm(ks, ks.pSpace)

dt = KitBase.timestep(ks, ctr, 0.0)
nt = ks.set.maxTime ÷ dt |> Int

@showprogress for iter = 1:nt
    @inbounds Threads.@threads for i in eachindex(face)
        vn = ks.vSpace.u .* face[i].n[1] .+ ks.vSpace.v .* face[i].n[2]
        vt = ks.vSpace.v .* face[i].n[1] .- ks.vSpace.u .* face[i].n[2]
        
        if !(-1 in ps.edgeCells[i, :])
            KitBase.flux_kfvs!(
                face[i].fw,
                face[i].fh,
                face[i].fb,
                ctr[ps.edgeCells[i, 1]].h,
                ctr[ps.edgeCells[i, 1]].b,
                ctr[ps.edgeCells[i, 2]].h,
                ctr[ps.edgeCells[i, 2]].b,
                vn,
                vt,
                ks.vSpace.weights,
                dt,
                face[i].len,
            )
            face[i].fw .= KitBase.global_frame(face[i].fw, face[i].n[1], face[i].n[2])
        else
            idx = ifelse(ps.edgeCells[i, 1] != -1, 1, 2)

            if ps.cellType[ps.edgeCells[i, idx]] == 2
                bc = KitBase.local_frame(ks.ib.primR, face[i].n[1], face[i].n[2])

                KitBase.flux_boundary_maxwell!(
                    face[i].fw,
                    face[i].fh,
                    face[i].fb,
                    bc,
                    ctr[ps.edgeCells[i, idx]].h,
                    ctr[ps.edgeCells[i, idx]].b,
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
            dirc = [sign(dot(ctr[i].n[j], face[ps.cellEdges[i, j]].n)) for j = 1:3]

            KitBase.step!(
                ctr[i].w,
                ctr[i].prim,
                ctr[i].h,
                ctr[i].b,
                face[ps.cellEdges[i, 1]].fw,
                face[ps.cellEdges[i, 1]].fh,
                face[ps.cellEdges[i, 1]].fb,
                face[ps.cellEdges[i, 2]].fw,
                face[ps.cellEdges[i, 2]].fh,
                face[ps.cellEdges[i, 2]].fb,
                face[ps.cellEdges[i, 3]].fw,
                face[ps.cellEdges[i, 3]].fh,
                face[ps.cellEdges[i, 3]].fb,
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

    if iter % 1000 == 0
        @save "ctr.jld2" ctr
        KitBase.write_vtk(ks, ctr)
    end
end

#KitBase.write_vtk(ks, ctr)
#@save "ctr.jld2" ctr