"""
3D simulation of spacecraft

sphere domain with R = 5 and center = (2, 0, 0)

structured mesh: 65×37×41
- i: index along positive x direction
- j: index around the annulus perpendicular to x
- k: index from spacecraft to the outer surface

"""

using CSV, DataFrames, KitBase

cd(@__DIR__)

ni, nj, nk = 65, 37, 41
df = CSV.read("../assets/mesh/spacecraft.csv", DataFrame)

x = reshape(collect(df.x), ni, nj, nk)
y = reshape(collect(df.y), ni, nj, nk)
z = reshape(collect(df.z), ni, nj, nk)

points = [[x[i, j, k], y[i, j, k], z[i, j, k]] for i = 1:65, j = 1:37, k = 1:41]
pc = [2., 0., 0.]

#--- geometric illustration ---#
[@show point_distance(points[i, j, end], pc) for i in axes(points, 1), j in axes(points, 2)]
[@show points[i, 1, end] for i in axes(points, 1)]
[@show points[30, j, end] for j in axes(points, 2)]
