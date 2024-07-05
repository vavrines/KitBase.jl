sL = 0.13
sR = -0.27

KitBase.linear(sL, sR)
KitBase.vanleer(sL, sR)
KitBase.vanleer(sL, sR, sR)
KitBase.vanleer(sL, sR, sR, 3)
KitBase.minmod(sL, sR)
KitBase.minmod(sL, sR, sR)
KitBase.minmod(sL, sR, sR, 3)
KitBase.vanalbaba(sL, sR)
KitBase.superbee(0.5, 0.9)
KitBase.superbee(-1.0, -0.8)
KitBase.superbee(0.51, 2.01)
KitBase.weno5(-2.0, -1.0, 0.0, 1.0, 2.0)

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

#--- filter ---#
let deg = 2, u = rand(deg + 1)
    ℓ = rand(deg + 1)

    modal_filter!(u, 1e-6; filter = :l2)
    modal_filter!(u, 1e-6; filter = :l2opt)
    modal_filter!(u, 1e-6, ℓ; filter = :l1)
    modal_filter!(u, ℓ; filter = :lasso)
    modal_filter!(u, 10; filter = :exp)
    modal_filter!(u, 10; filter = :houli)
end

let deg = 2, u = rand(deg + 1, deg + 1)
    ℓ = rand(deg + 1, deg + 1)

    modal_filter!(u, 1e-6, 1e-6; filter = :l2)
    modal_filter!(u, 1e-6, 1e-6; filter = :l2opt)
    modal_filter!(u, 1e-6, 1e-6, ℓ; filter = :l1)
    modal_filter!(u, ℓ; filter = :lasso)
    modal_filter!(u, 2; filter = :exp)
end
