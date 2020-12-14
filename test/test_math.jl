x0 = 0
x1 = 1
nx = 5

linspace(x0, x1, nx)
heaviside(x0)
fortsign(x0, x1)

m = randn(2, 2)
mat_split(m)

x = randn(16)
y = randn(16)
central_diff(y, x)
central_diff(y, x0)

res = similar(x)
central_diff!(res, y, x)
central_diff!(res, y, x0)

upwind_diff(y, x)
upwind_diff(y, x0)
upwind_diff!(res, y, x)
upwind_diff!(res, y, x0)

unstruct_diff(y, x, 4)
unstruct_diff(sin, randn(12), 4, 1)

KitBase.lgwt(12, -1, 1)
KitBase.extract_last(randn(2, 3), 2; mode = :view)

find_idx(randn(20), 0.13)