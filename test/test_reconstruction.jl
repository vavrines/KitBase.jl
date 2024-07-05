sL = 0.13
sR = -0.27

KB.linear(sL, sR)
KB.vanleer(sL, sR)
KB.vanleer(sL, sR, sR)
KB.vanleer(sL, sR, sR, 3)
KB.minmod(sL, sR)
KB.minmod(sL, sR, sR)
KB.minmod(sL, sR, sR, 3)
KB.vanalbaba(sL, sR)
KB.superbee(0.5, 0.9)
KB.superbee(-1.0, -0.8)
KB.superbee(0.51, 2.01)
KB.weno5(-2.0, -1.0, 0.0, 1.0, 2.0)

wL = rand(3)
wR = randn(3)
dx = 1e-2

KB.reconstruct2(sL, sR, dx)
KB.reconstruct2(wL, wR, dx)
KB.reconstruct2(randn(3, 2), randn(3, 2), dx)
KB.reconstruct2(randn(3, 2, 2), randn(3, 2, 2), dx)

KB.reconstruct2!(similar(wL), wL, wR, dx)
KB.reconstruct2!(randn(3, 2), randn(3, 2), randn(3, 2), dx)
KB.reconstruct2!(randn(3, 2, 2), randn(3, 2, 2), randn(3, 2, 2), dx)

wN = rand(3)

KB.reconstruct3(sL, rand(), sR, dx, dx, :linear)
KB.reconstruct3(sL, rand(), sR, dx, dx, :vanleer)
KB.reconstruct3(sL, rand(), sR, dx, dx, :minmod)
KB.reconstruct3(wL, wN, wR, dx, dx, :linear)
KB.reconstruct3(wL, wN, wR, dx, dx, :vanleer)
KB.reconstruct3(wL, wN, wR, dx, dx, :minmod)
KB.reconstruct3(randn(3, 2), randn(3, 2), randn(3, 2), dx, dx)
KB.reconstruct3(randn(3, 2, 2), randn(3, 2, 2), randn(3, 2, 2), dx, dx)

KB.reconstruct3!(similar(wL), wL, wN, wR, dx, dx, :linear)
KB.reconstruct3!(similar(wL), wL, wN, wR, dx, dx, :vanleer)
KB.reconstruct3!(similar(wL), wL, wN, wR, dx, dx, :minmod)
KB.reconstruct3!(randn(3, 2), randn(3, 2), randn(3, 2), randn(3, 2), dx, dx)
KB.reconstruct3!(
    randn(3, 2, 2),
    randn(3, 2, 2),
    randn(3, 2, 2),
    randn(3, 2, 2),
    dx,
    dx,
)
KB.reconstruct3!(
    randn(3, 2, 2, 2),
    randn(3, 2, 2, 2),
    randn(3, 2, 2, 2),
    randn(3, 2, 2, 2),
    dx,
    dx,
)

KB.weno5(-2.0, -1.0, 0.0, 1.0, 2.0)

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
