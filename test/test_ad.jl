using KitBase

u = -5:0.2:5 |> collect
prim = [1.0, 0.0, 1.0]

KB.∂maxwellian(u[1], prim)
KB.∂maxwellian(u, prim)
