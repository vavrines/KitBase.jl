using LinearAlgebra, WriteVTK, ProgressMeter, JLD2
import KitBase
#=
function KitBase.write_vtk(ps::KitBase.AbstractPhysicalSpace, ctr)
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

begin
    cd(@__DIR__)
    cells, points = KitBase.read_mesh("../assets/mesh/naca0012.su2")
    cellid = KitBase.extract_cell(cells)
    edgePoints, edgeCells, cellNeighbors = KitBase.mesh_connectivity_2D(cellid)
    cellType = KitBase.mesh_cell_type(cellNeighbors)
    cellArea = KitBase.mesh_area_2D(points, cellid)
    cellCenter = KitBase.mesh_center_2D(points, cellid)
    edgeCenter = KitBase.mesh_edge_center(points, edgePoints)
    cellEdges = KitBase.mesh_cell_edge(cellid, edgeCells)

    for i in eachindex(cellType)
        if cellType[i] == 1 && -0.5 < cellCenter[i, 2] < 0.5
            cellType[i] = 2
        end
    end

    ps = KitBase.UnstructPSpace(cells, points, cellid, cellType, cellNeighbors, cellEdges, cellCenter, cellArea, edgePoints, edgeCells, edgeCenter)
end

set = KitBase.Setup(
    "gas",
    "naca",
    "2d2f2v",
    "kfvs",
    "bgk",
    1,
    1,
    "vanleer",
    "maxwell",
    0.8,
    5.0,
)

vs = KitBase.VSpace2D(
    -5.0,
    5.0,
    24,
    -5.0,
    5.0,
    24,
)

Kn = 0.075
Ma = 0.8
gas = KitBase.Gas(
    Kn,
    Ma,
    1.0,
    1.0,
    5 / 3,
    0.5,
    1.0,
    0.81,
    KitBase.ref_vhs_vis(Kn, 1.0, 0.5),
)

primL = [1.0, 0.7, 0.0, 1.0]
wL = KitBase.prim_conserve(primL, 5/3)
hL = KitBase.maxwellian(vs.u, vs.v, primL)
bL = @. hL * gas.K / 2 / primL[end]
primR = [1.0, 0.0, 0.0, 1.0]
wR = KitBase.prim_conserve(primR, 5/3)
hR = KitBase.maxwellian(vs.u, vs.v, primR)
bR = @. hR * gas.K / 2 / primR[end]

ib = KitBase.IB2F(
    wL,
    primL,
    hL,
    bL,
    primL,
    wR,
    primR,
    hR,
    bR,
    primR,
)

ks = KitBase.SolverSet(set, ps, vs, gas, ib, @__DIR__)

ctr = Array{KitBase.ControlVolumeUS2F}(undef, size(ps.cells, 1))
for i in eachindex(ctr)
    n = Vector{Float64}[]
    for j = 1:3
        push!(n, KitBase.unit_normal(ps.points[edgePoints[cellEdges[i, j], 1], :], ps.points[edgePoints[cellEdges[i, j], 2], :]))

        if dot(ps.edgeCenter[ps.cellEdges[i, j], :] .- ps.cellCenter[i, :], n[j]) < 0
            n[j] .= -n[j]
        end
    end

    dx = [
        KitBase.point_distance(cellCenter[i, :], ps.points[ps.cells[i, 1], :], ps.points[ps.cells[i, 2], :]),
        KitBase.point_distance(cellCenter[i, :], ps.points[ps.cells[i, 2], :], ps.points[ps.cells[i, 3], :]),
        KitBase.point_distance(cellCenter[i, :], ps.points[ps.cells[i, 3], :], ps.points[ps.cells[i, 1], :]),
    ]

    ctr[i] = KitBase.ControlVolumeUS2F(
        n,
        cellCenter[i, :],
        dx,
        wL,
        primL,
        hL,
        bL,
    )
end

face = Array{KitBase.Interface2D2F}(undef, size(ps.edgePoints, 1))
for i in eachindex(face)
    len = norm(ps.points[edgePoints[i, 1], :] .- ps.points[edgePoints[i, 2], :])
    n = KitBase.unit_normal(ps.points[edgePoints[i, 1], :], ps.points[edgePoints[i, 2], :])
    
    if !(-1 in ps.edgeCells[i, :])
        n0 = ps.cellCenter[ps.edgeCells[i, 2], :] .- ps.cellCenter[ps.edgeCells[i, 1], :]
    else
        idx = ifelse(ps.edgeCells[i, 1] != -1, ps.edgeCells[i, 1], ps.edgeCells[i, 2])
        n0 = ps.cellCenter[idx, :] .- ps.edgeCenter[i, :]
    end
    if dot(n, n0) < 0
        n .= -n
    end
    
    fw = zero(ks.ib.wL)
    fh = zero(ks.ib.hL)

    face[i] = KitBase.Interface2D2F(
        len,
        n[1],
        n[2],
        fw,
        fh,
    )
end

@load "ctr.jld2" ctr

dt = KitBase.timestep(ks, ctr, 0.0)
nt = ks.set.maxTime ÷ dt |> Int

@showprogress for iter = 1:10000#nt
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

    @inbounds Threads.@threads for i in eachindex(ctr)
        if ps.cellType[i] in (0, 2)
            MH = KitBase.maxwellian(vs.u, vs.v, ctr[i].prim)
            MB = MH .* ks.gas.K ./ (2.0 * ctr[i].prim[end])
            τ = KitBase.vhs_collision_time(ctr[i].prim, ks.gas.μᵣ, ks.gas.ω)

            for j in 1:3
                dirc = sign(dot(ctr[i].n[j], face[ps.cellEdges[i, j]].n))
                @. ctr[i].w -= dirc * face[ps.cellEdges[i, j]].fw / ps.cellArea[i]
                @. ctr[i].h -= dirc * face[ps.cellEdges[i, j]].fh / ps.cellArea[i]
                @. ctr[i].b -= dirc * face[ps.cellEdges[i, j]].fb / ps.cellArea[i]
            end
            ctr[i].prim .= KitBase.conserve_prim(ctr[i].w, ks.gas.γ)

            @. ctr[i].h += (MH - ctr[i].h) * dt / τ
            @. ctr[i].b += (MB - ctr[i].b) * dt / τ
            #=for q in axes(vs.v, 2), p in axes(vs.u, 1)
                ctr[i].h[p, q] = (ctr[i].h[p, q] - (face[ps.cellEdges[i, 1]].fh[p, q] + face[ps.cellEdges[i, 2]].fh[p, q] + face[ps.cellEdges[i, 3]].fh[p, q]) / ps.cellArea[i] + 
                    dt / τ * MH[p, q]) / (1.0 + dt / τ)
                ctr[i].b[p, q] = (ctr[i].b[p, q] - (face[ps.cellEdges[i, 1]].fb[p, q] + face[ps.cellEdges[i, 2]].fb[p, q] + face[ps.cellEdges[i, 3]].fb[p, q]) / ps.cellArea[i] + 
                    dt / τ * MB[p, q]) / (1.0 + dt / τ)
            end=#
        end
    end
end

KitBase.write_vtk(ks, ctr)
#KitBase.write_vtk(ps, ctr)

@save "ctr.jld2" ctr