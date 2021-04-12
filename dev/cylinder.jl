using LinearAlgebra, WriteVTK, ProgressMeter, JLD2, PyCall
import KitBase

cd(@__DIR__)
meshio = pyimport("meshio")
m = meshio.read("../assets/mesh/cylinder.msh")

cells, points = KitBase.read_mesh("../assets/mesh/cylinder.msh")
cellid = KitBase.extract_cell(cells)
edgePoints, edgeCells, cellNeighbors = KitBase.mesh_connectivity_2D(cellid)
cellType = KitBase.mesh_cell_type(cellNeighbors)
cellArea = KitBase.mesh_area_2D(points, cellid)
cellCenter = KitBase.mesh_center_2D(points, cellid)
edgeCenter = KitBase.mesh_edge_center(points, edgePoints)
cellEdges = KitBase.mesh_cell_edge(cellid, edgeCells)
edgeType = KitBase.mesh_edge_type(edgeCells, cellType)

ps = KitBase.UnstructPSpace(cells, points, cellid, cellType, cellNeighbors, cellEdges, cellCenter, cellArea, edgePoints, edgeCells, edgeCenter, edgeType)
