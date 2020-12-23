u = collect(-5:0.1:5)

KitBase.ib_rh(2., 5 / 3, u)
KitBase.ib_rh(2., 5 / 3, 2., u)
KitBase.ib_rh(2., 5 / 3, rand(21, 21, 21), rand(21, 21, 21), rand(21, 21, 21))
KitBase.ib_rh(2., 5 / 3, 2., 1., 0.5, 1., 1., hcat(u, u))

KitBase.ib_sod(5 / 3, u)
KitBase.ib_sod(5 / 3, rand(21, 21, 21), rand(21, 21, 21), rand(21, 21, 21))
KitBase.ib_sod(5 / 3, 2., u)

KitBase.ib_cavity(5 / 3, 1., 0., 1., rand(21, 21), rand(21, 21))
KitBase.ib_cavity(5 / 3, 2., 1., 0., 1., rand(21, 21), rand(21, 21))

KitBase.ib_briowu(5 / 3, 1., 0.5, hcat(u, u))
