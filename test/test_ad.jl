using KitBase

u = -5:0.2:5 |> collect
prim = [1.0, 0.0, 1.0]

KitBase.∂maxwellian(u[1], prim)
KitBase.∂maxwellian(u, prim)
