w = [1.0, 0.1, 0.3, 1.0]
cosa = cos(0.5)
sina = sin(0.5)
global_frame(w, cosa, sina)
local_frame(w, cosa, sina)

PSpace1D()
PSpace2D()

x0 = 0.0
nx = 15
dx = 0.1
uniform_mesh(x0, nx, dx)

x = collect(0:0.1:1)
y = collect(0:0.1:1)
z = collect(0:0.1:1)
meshgrid(x, y)
meshgrid(x, y, z)

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

UnstructMesh(nodes, cells)
mesh_connectivity_2D(cells)
mesh_area_2D(nodes, cells)
mesh_center_2D(nodes, cells)