w = [1.0, 0.1, 0.3, 1.0]
cosa = cos(0.5)
sina = sin(0.5)
KitBase.global_frame(randn(2), cosa, sina)
KitBase.local_frame(randn(2), cosa, sina)
KitBase.global_frame(w, cosa, sina)
KitBase.local_frame(w, cosa, sina)

KitBase.global_frame(rand(3), rand(3, 3))
KitBase.local_frame(rand(3), rand(3, 3))
KitBase.global_frame(rand(5), rand(3, 3))
KitBase.local_frame(rand(5), rand(3, 3))

KitBase.unit_normal(rand(2), rand(2))
KitBase.unit_normal(rand(3), rand(3), rand(3))
KitBase.point_distance(rand(2), rand(2), rand(2))

#--- structure mesh ---#
KitBase.PSpace1D()
KitBase.PSpace1D(0.0, 1.0) |> show
KitBase.PSpace2D()
KitBase.PSpace2D(0.0, 1.0, 0.0, 1.0) |> show

x0 = 0.0
nx = 15
dx = 0.1
KitBase.uniform_mesh(x0, nx, dx)

x = collect(0:0.1:1)
y = collect(0:0.1:1)
z = collect(0:0.1:1)
KitBase.ndgrid(x)
KitBase.ndgrid(x, y)
KitBase.ndgrid(x, y, z)
KitBase.meshgrid(x)
KitBase.meshgrid(x, y)
KitBase.meshgrid(x, y, z)

KitBase.find_idx(randn(20), 0.13, mode = :uniform)
KitBase.find_idx(randn(20), 0.13, mode = :nonuniform)

#--- unstructure mesh ---#
cd(@__DIR__)
cellid, nodes = KitBase.read_mesh("t1.msh")
edgeNodes, edgeCells, cellNeighbors = KitBase.mesh_connectivity_2D(cellid)
cellType = KitBase.mesh_cell_type(cellNeighbors)
cellArea = KitBase.mesh_area_2D(nodes, cellid)
cellCenter = KitBase.mesh_center_2D(nodes, cellid)
edgeCenter = KitBase.mesh_edge_center(nodes, edgeNodes)
cellEdges = KitBase.mesh_cell_edge(cellid, edgeCells)
KitBase.UnstructPSpace(cells, nodes, cellType, cellNeighbors, cellEdges, cellCenter, cellArea, edgeNodes, edgeCells, edgeCenter)
KitBase.write_vtk(nodes, cellid, randn(size(cellid, 1)))

cells, nodes = KitBase.read_mesh("square.msh")
edgeNodes, edgeCells, cellNeighbors = KitBase.mesh_connectivity_2D(cells)
cellArea = KitBase.mesh_area_2D(nodes, cells)
