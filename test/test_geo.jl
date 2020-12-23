w = [1.0, 0.1, 0.3, 1.0]
cosa = cos(0.5)
sina = sin(0.5)
KitBase.global_frame(w, cosa, sina)
KitBase.local_frame(w, cosa, sina)

#--- structure mesh ---#
KitBase.PSpace1D()
KitBase.PSpace2D()

x0 = 0.0
nx = 15
dx = 0.1
KitBase.uniform_mesh(x0, nx, dx)

x = collect(0:0.1:1)
y = collect(0:0.1:1)
z = collect(0:0.1:1)
KitBase.meshgrid(x, y)
KitBase.meshgrid(x, y, z)

KitBase.find_idx(randn(20), 0.13, mode = :uniform)
KitBase.find_idx(randn(20), 0.13, mode = :nonuniform)

#--- unstructure mesh ---#
cd(@__DIR__)
#=
using PyCall
try
    meshio = pyimport("meshio")
catch
    using Conda
    Conda.add_channel("conda-forge")
    Conda.add("meshio")
    meshio = pyimport("meshio")
end
nodes, cells = read_mesh("t1.msh")
=#
using JLD2
@load "t1.jld2" nodes cells

KitBase.UnstructMesh(nodes, cells)
KitBase.mesh_connectivity_2D(cells)
KitBase.mesh_area_2D(nodes, cells)
KitBase.mesh_center_2D(nodes, cells)
