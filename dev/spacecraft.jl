"""
3D simulation of spacecraft

sphere domain with R = 5 and center = (2, 0, 0)

structured mesh: 65×37×41
- i: index along positive x direction
- j: index around the annulus perpendicular to x
- k: index from spacecraft to the outer surface

"""

using KitBase, CSV, DataFrames

cd(@__DIR__)

set = (
    ni = 65,
    nj = 37,
    nk = 41,
    xc = [2., 0., 0.],
)

function point_coords(df)
    ni, nj, nk = set.ni, set.nj, set.nk
    x = reshape(collect(df.x), ni, nj, nk)
    y = reshape(collect(df.y), ni, nj, nk)
    z = reshape(collect(df.z), ni, nj, nk)
    points = cat(x, y, z, dims = 4)

    return points
end

df = CSV.read("../assets/mesh/spacecraft.csv", DataFrame)
points = point_coords(df)

#--- geometric illustration ---#
[@show point_distance(points[i, j, end, :], set.xc) for i in axes(points, 1), j in axes(points, 2)]
[@show points[i, 1, end, :] for i in axes(points, 1)]
[@show points[30, j, end, :] for j in axes(points, 2)]

#--- physical space ---#
function physical_space(points)
    ni, nj, nk = size(points)[1:3] .- 1
    x0, y0, z0 = minimum(points[:, :, :, 1]), minimum(points[:, :, :, 2]), minimum(points[:, :, :, 3])
    x1, y1, z1 = maximum(points[:, :, :, 1]), maximum(points[:, :, :, 2]), maximum(points[:, :, :, 3])
end

physical_space(points)

ni, nj, nk = size(points)[1:3] .- 1
x0, y0, z0 = minimum(points[:, :, :, 1]), minimum(points[:, :, :, 2]), minimum(points[:, :, :, 3])
x1, y1, z1 = maximum(points[:, :, :, 1]), maximum(points[:, :, :, 2]), maximum(points[:, :, :, 3])
x = zeros(ni, nj, nk)
y = zero(x)
z = zero(x)
dx = zero(x)
dy = zero(x)
dz = zero(x)
for i = 1:ni, j = 1:nj, k = 1:nk
    x[i, j, k] = points[i, j, k, 1]
end