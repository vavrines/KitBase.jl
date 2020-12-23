x0 = 0
x1 = 1
nx = 5

KitBase.linspace(x0, x1, nx)
KitBase.heaviside(x0)
KitBase.fortsign(x0, x1)

m = randn(2, 2)
KitBase.mat_split(m)

x = randn(16)
y = randn(16)
KitBase.central_diff(y, x)
KitBase.central_diff(y, x0)

res = similar(x)
KitBase.central_diff!(res, y, x)
KitBase.central_diff!(res, y, x0)

KitBase.upwind_diff(y, x)
KitBase.upwind_diff(y, x0)
KitBase.upwind_diff!(res, y, x)
KitBase.upwind_diff!(res, y, x0)

KitBase.unstruct_diff(y, x, 4)
KitBase.unstruct_diff(sin, randn(12), 4, 1)

KitBase.lgwt(12, -1, 1)
KitBase.extract_last(randn(2, 3), 2; mode = :view)

KitBase.find_idx(randn(20), 0.13)
