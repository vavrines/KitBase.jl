sL = 0.13
sR = -0.27

KitBase.linear(sL, sR)
KitBase.vanleer(sL, sR)
KitBase.vanleer(sL, sR, sR)
KitBase.minmod(sL, sR)
KitBase.minmod(sL, sR, sR)
KitBase.vanalbaba(sL, sR)
KitBase.superbee(0.51, 1.99)
KitBase.superbee(0.49, 2.01)
KitBase.superbee(0.51, 2.01)

wL = rand(3)
wR = randn(3)
dx = 1e-2

KitBase.reconstruct2(sL, sR, dx)
KitBase.reconstruct2(wL, wR, dx)
KitBase.reconstruct2(randn(3, 2), randn(3, 2), dx)
KitBase.reconstruct2(randn(3, 2, 2), randn(3, 2, 2), dx)

KitBase.reconstruct2!(similar(wL), wL, wR, dx)
KitBase.reconstruct2!(randn(3, 2), randn(3, 2), randn(3, 2), dx)
KitBase.reconstruct2!(randn(3, 2, 2), randn(3, 2, 2), randn(3, 2, 2), dx)

wN = rand(3)

KitBase.reconstruct3(sL, rand(), sR, dx, dx, :linear)
KitBase.reconstruct3(sL, rand(), sR, dx, dx, :vanleer)
KitBase.reconstruct3(sL, rand(), sR, dx, dx, :minmod)
KitBase.reconstruct3(wL, wN, wR, dx, dx, :linear)
KitBase.reconstruct3(wL, wN, wR, dx, dx, :vanleer)
KitBase.reconstruct3(wL, wN, wR, dx, dx, :minmod)
KitBase.reconstruct3(randn(3, 2), randn(3, 2), randn(3, 2), dx, dx)
KitBase.reconstruct3(randn(3, 2, 2), randn(3, 2, 2), randn(3, 2, 2), dx, dx)

KitBase.reconstruct3!(similar(wL), wL, wN, wR, dx, dx, :linear)
KitBase.reconstruct3!(similar(wL), wL, wN, wR, dx, dx, :vanleer)
KitBase.reconstruct3!(similar(wL), wL, wN, wR, dx, dx, :minmod)
KitBase.reconstruct3!(randn(3, 2), randn(3, 2), randn(3, 2), randn(3, 2), dx, dx)
KitBase.reconstruct3!(
    randn(3, 2, 2),
    randn(3, 2, 2),
    randn(3, 2, 2),
    randn(3, 2, 2),
    dx,
    dx,
)
KitBase.reconstruct3!(
    randn(3, 2, 2, 2),
    randn(3, 2, 2, 2),
    randn(3, 2, 2, 2),
    randn(3, 2, 2, 2),
    dx,
    dx,
)

KitBase.weno5(-2.0, -1.0, 0.0, 1.0, 2.0)
