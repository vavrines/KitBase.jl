import KitBase.AA

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
KitBase.ps1D()
KitBase.ps1D(0.0, 1.0) |> show
KitBase.ps2D()
KitBase.ps2D(0.0, 1.0, 0.0, 1.0) |> show
KitBase.CSpace2D(0.0, 1.0, 10, 0.0, π, 10) |> show
KitBase.ps3D{Int,Int,AA{Float64,3},AA{Float64,5},AA{Float64,4},AA{Float64,4}}(
    0,
    1,
    10,
    0,
    1,
    10,
    0,
    1,
    10,
    zeros(10, 10, 10),
    zeros(10, 10, 10),
    zeros(10, 10, 10),
    zeros(10, 10, 10),
    zeros(10, 10, 10),
    zeros(10, 10, 10),
    zeros(10, 10, 10, 8, 3),
    zeros(10, 10, 10, 6),
    zeros(10, 10, 10, 6),
) |> show

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
ps = KitBase.UnstructPSpace("t1.msh")
KitBase.write_vtk(ps.points, ps.cellid, randn(size(ps.cellid, 1)))
