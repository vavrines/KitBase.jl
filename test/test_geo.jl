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
KitBase.PSpace1D()
KitBase.PSpace1D(0.0, 1.0) |> show
KitBase.PSpace2D()
KitBase.PSpace2D(0.0, 1.0, 0.0, 1.0) |> show
KitBase.CSpace2D(0.0, 1.0, 10, 0.0, π, 10) |> show
KitBase.PSpace3D{Int,Int,AA{Float64,3},AA{Float64,5},AA{Float64,4},AA{Float64,4}}(
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

#--- IB ---#
ps = KB.PSpace2D(0, 3, 15, 0, 2, 10, 1, 1)
radius = 1

flags = ones(Int, axes(ps.x))
for i in axes(flags, 1), j in axes(flags, 2)
    if (ps.x[i, j] - 3)^2 + ps.y[i, j]^2 < radius # (x-3)^2+y^2=1
        flags[i, j] = 0
    end
end
flags[0, :] .= -1
flags[ps.nx+1, :] .= -1
flags[:, 0] .= -1
flags[:, ps.ny+1] .= -1
KB.ghost_flag!(ps, flags)

ghost_ids = findall(flags .== -2)
xbis = [Vector{Float64}(undef, 2) for iter = 1:length(ghost_ids)]
nbis = zero.(xbis)
for iter in axes(xbis, 1)
    idx = ghost_ids[iter]

    θ = atan(ps.y[idx] / (3 - ps.x[idx]))
    xbis[iter][1] = 3 - radius * cos(θ)
    xbis[iter][2] = radius * sin(θ)

    nbis[iter] .= [- radius * cos(θ), radius * sin(θ)]
end
xips = KB.ip_location(ps, ghost_ids, xbis)
ip_cids, ip_nids, ip_bids = KB.ip_connectivity(ps, xips, flags)
ib = SharpIB(flags, ghost_ids, xbis, nbis, xips, ip_cids, ip_nids, ip_bids)
KB.bilinear_coeffs(ps, ib, 1, rand(2))
