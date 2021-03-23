using LinearAlgebra, WriteVTK
import KitBase

begin
    cd(@__DIR__)

    # time
    tspan = (0.0, 1.0)
    cfl = 0.8

    # quadrature
    quadratureorder = 6
    points, triangulation = KitBase.octa_quadrature(quadratureorder)
    weights = KitBase.quadrature_weights(points, triangulation)
    nq = size(points, 1)

    vs = KitBase.UnstructVSpace(-1.0, 1.0, nq, points, ones(nq) .* 2 / nq)

    # IC
    s2 = 0.03^2
    init_field(x, y) = max(1e-4, 1.0 / (4.0 * π * s2) * exp(-(x^2 + y^2) / 4.0 / s2))

    # geometry
    cells, points = KitBase.read_mesh("../assets/mesh/linesource.su2")
    cellid = KitBase.extract_cell(cells)
    edgePoints, edgeCells, cellNeighbors = KitBase.mesh_connectivity_2D(cellid)
    cellType = KitBase.mesh_cell_type(cellNeighbors)
    cellArea = KitBase.mesh_area_2D(points, cellid)
    cellCenter = KitBase.mesh_center_2D(points, cellid)
    edgeCenter = KitBase.mesh_edge_center(points, edgePoints)

    ps = KitBase.UnstructPSpace(cells, points, cellid, cellType, cellNeighbors, cellArea, cellCenter, edgePoints, edgeCells, edgeCenter)

    # particle
    SigmaS = ones(size(cellid, 1))
    SigmaA = zeros(size(cellid, 1))
    SigmaT = SigmaS + SigmaA
end

KitBase.write_vtk(ps.points, ps.cellid, randn(8442, 3))

ctr = Array{KitBase.ControlVolumeUS1F}(undef, size(ps.cells, 1))
for i in eachindex(ctr)
    phi = zeros(nq)
    phi .= init_field(ps.cellCenter[i, 1], ps.cellCenter[i, 2]) / 4.0 / π
    w = sum(weights .* phi)
    dx = [
        KitBase.pl_distance(cellCenter[i, :], nodes[cells[i, 1], :], nodes[cells[i, 2], :]),
        KitBase.pl_distance(cellCenter[i, :], nodes[cells[i, 2], :], nodes[cells[i, 3], :]),
        KitBase.pl_distance(cellCenter[i, :], nodes[cells[i, 3], :], nodes[cells[i, 1], :]),
    ]

    ctr[i] = KitBase.ControlVolumeUS1F(
        cellCenter[i, :],
        dx,
        w,
        w,
        phi
    )
end

face = Array{KitBase.Interface2D1F}(undef, size(edgeNodes, 1))
for i in eachindex(face)
    len = norm(nodes[edgeNodes[i, 1], :] .- nodes[edgeNodes[i, 1], :])
    n = KitBase.unit_normal(nodes[edgeNodes[i, 1], :], nodes[edgeNodes[i, 2], :])

    fw = 0.0
    ff = zeros(nq)

    face[i] = KitBase.Interface2D1F(
        len,
        n[1],
        n[2],
        fw,
        ff,
    )
end
