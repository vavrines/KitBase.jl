using LinearAlgebra
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

    # IC
    s2 = 0.03^2
    init_field(x, y) = max(1e-4, 1.0 / (4.0 * π * s2) * exp(-(x^2 + y^2) / 4.0 / s2))

    # geometry
    nodes, cells = KitBase.read_mesh("../assets/mesh/linesource.su2")
    edgeNodes, edgeCells, cellNeighbors = KitBase.mesh_connectivity_2D(cells)
    cellid = KitBase.mesh_cell_type(cellNeighbors)
    cellArea = KitBase.mesh_area_2D(nodes, cells)
    cellCenter = KitBase.mesh_center_2D(nodes, cells)
    edgeCenter = KitBase.mesh_edge_center(nodes, edgeNodes)

    ps = KitBase.UnstructMesh(nodes, cells, cellid, cellNeighbors, cellArea, cellCenter, edgeNodes, edgeCells, edgeCenter)

    # particle
    SigmaS = ones(size(cells, 1))
    SigmaA = zeros(size(cells, 1))
    SigmaT = SigmaS + SigmaA
end

using WriteVTK

vtkfile = vtk_grid("my_vtk_file", points, cells)


MeshCell(VTKCellTypes.VTK_TRIANGLE, cells)


cd(@__DIR__)
mesh = meshio.read("../assets/mesh/linesource.su2")
mesh.write("foo.vtk")




using PyCall
meshio = pyimport("meshio")
[("triangle", cells)]


meshio.Mesh(
    points,
    [("triangle", cells)]
    # Optionally provide extra data on points, cells, etc.
    # point_data=point_data,
    # cell_data=cell_data,
    # field_data=field_data
).write(
    "foo.vtk",  # str, os.PathLike, or buffer/open file
    # file_format="vtk",  # optional if first argument is a path; inferred from extension
)



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

