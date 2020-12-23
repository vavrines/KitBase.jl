inK = 2.
γ = 5 / 3
wL = [1.13, 0.13, 0.8]
wR = [0.3, -0.2, 1.5]
primL = KitBase.conserve_prim(wL, γ)
primR = KitBase.conserve_prim(wR, γ)

u = collect(-5.:0.2:5.)
ω = ones(length(u)) / length(u)

hL = KitBase.maxwellian(u, primL)
bL = hL .* inK ./ (2.0 * primL[end])
hR = KitBase.maxwellian(u, primR)
bR = hR .* inK ./ (2.0 * primR[end])

dt = 1e-3
dx = 1e-2

fw = similar(wL)
fh = similar(hL)
fb = similar(bL)

KitBase.flux_lax!(fw, wL, wR, γ, dt, dx)
KitBase.flux_hll!(fw, wL, wR, γ, dt)
KitBase.flux_roe!(fw, wL, wR, γ, dt)
KitBase.flux_roe!(zeros(4), [1., 0.3, 0., 1.], [0.3, -0.1, 0., 2.], γ, dt)

KitBase.flux_gks(0.3, 1e-3, dt, 1e-1, 0)
KitBase.flux_gks(0.3, 1e-3, dt, 1e-1, 1.)
KitBase.flux_gks(1., 0.125, 1e-3, dt, 1e-2, 1e-2, 0.1, 0.3)
KitBase.flux_gks!(fw, wL, wR, γ, inK, 1e-3, 0.72, dt, dx)
KitBase.flux_gks!(zeros(4), [1., 0.3, 0., 1.], [0.3, -0.1, 0., 2.], γ, inK, 1e-3, 0.72, dt, dx, dx)
KitBase.flux_gks!(fw, fh, fb, wL, wR, u, inK, γ, 1e-3, 0.72, dt, dx, dx)
KitBase.flux_gks!(
    zeros(4, 2),
    hcat([1., 0.3, 0., 1.], [1., 0.3, 0., 1.]),
    hcat([0.3, -0.1, 0., 2.], [0.3, -0.1, 0., 2.]),
    inK,
    γ,
    1.,
    1.,
    0.5,
    1.,
    1e-2,
    dt,
    dx,
    dx,
    dx,
)

KitBase.flux_ugks!(
    fw,
    fh,
    fb,
    wL,
    hL,
    bL,
    wR,
    hR,
    bR,
    u,
    ω,
    inK,
    γ,
    1e-3,
    0.81,
    0.72,
    dt,
    dx,
    dx,
)
KitBase.flux_ugks!(
    zeros(4),
    zeros(28, 28),
    zeros(28, 28),
    [1., 0.3, 0., 1.],
    rand(28, 28),
    rand(28, 28),
    [0.3, -0.1, 0., 2.],
    rand(28, 28),
    rand(28, 28),
    rand(28, 28),
    rand(28, 28),
    ones(28, 28),
    inK,
    γ,
    1e-3,
    0.81,
    0.72,
    dt,
    dx,
    dx,
    dx
)
KitBase.flux_ugks!(
    zeros(5, 2),
    zeros(16, 16, 2),
    zeros(16, 16, 2),
    zeros(16, 16, 2),
    hcat([1., 0.3, 0., 0., 1.], [1., 0.3, 0., 0., 1.]),
    rand(16, 16, 2),
    rand(16, 16, 2),
    rand(16, 16, 2),
    hcat([0.3, -0.1, 0., 0., 2.], [0.3, -0.1, 0., 0., 2.]),
    rand(16, 16, 2),
    rand(16, 16, 2),
    rand(16, 16, 2),
    rand(16, 16, 2),
    rand(16, 16, 2),
    ones(16, 16, 2),
    inK,
    γ,
    1.,
    1.,
    0.5,
    1.,
    0.1,
    dt,
    dx,
    dx,
    dx,
)