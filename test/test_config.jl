u = collect(-5:0.1:5)

KitBase.ib_rh(2.0, 5 / 3)
KitBase.ib_rh(2.0, 5 / 3, u)
KitBase.ib_rh(2.0, 5 / 3, 2.0, u)
KitBase.ib_rh(2.0, 5 / 3, rand(21, 21, 21), rand(21, 21, 21), rand(21, 21, 21))
KitBase.ib_rh(2.0, 5 / 3, 2.0, 1.0, 0.5, 1.0, 1.0, hcat(u, u))

KitBase.ib_sod(5 / 3, u)
KitBase.ib_sod(5 / 3, rand(21, 21, 21), rand(21, 21, 21), rand(21, 21, 21))
KitBase.ib_sod(5 / 3, 2.0, u)

KitBase.ib_cavity(5 / 3, 1.0, 0.0, 1.0, rand(21, 21), rand(21, 21))
KitBase.ib_cavity(5 / 3, 2.0, 1.0, 0.0, 1.0, rand(21, 21), rand(21, 21))

KitBase.ib_briowu(5 / 3, 1.0, 0.5, hcat(u, u))
