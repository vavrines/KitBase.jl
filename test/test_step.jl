fwL = zeros(3)
fwR = zeros(3)
γ = 5 / 3
prim = [1.0, 0.0, 1.0]
w = KitBase.prim_conserve(prim, γ)
dx = 1e-2
dt = 1e-3
res = zeros(3)
avg = zeros(3)

KitBase.step!(fwL, w, prim, fwR, γ, dx, res, avg)

ffL = zeros(16)
ffR = zeros(16)
f = rand(16)
u = randn(16)
ω = ones(16)

KitBase.step!(fwL, ffL, w, prim, f, fwR, ffR, u, ω, γ, 1e-3, 0.72, 1.0, dx, dt, res, avg)

KitBase.step!(
    zeros(5),
    zeros(16, 16, 16),
    KitBase.prim_conserve([1.0, 0.0, 0.0, 0.0, 1.0], γ),
    [1.0, 0.0, 0.0, 0.0, 1.0],
    rand(16, 16, 16),
    zeros(5),
    zeros(16, 16, 16),
    randn(16, 16, 16),
    randn(16, 16, 16),
    randn(16, 16, 16),
    ones(16, 16, 16),
    γ,
    1e-3,
    0.72,
    1.0,
    dx,
    dt,
    zeros(5),
    zeros(5),
)

KitBase.step!(
    fwL,
    ffL,
    ffL,
    w,
    prim,
    f,
    f,
    fwR,
    ffR,
    ffR,
    u,
    ω,
    2,
    γ,
    1e-3,
    0.72,
    1.0,
    dx,
    dt,
    res,
    avg,
)
KitBase.step!(
    zeros(3, 2),
    zeros(16, 2),
    zeros(16, 2),
    hcat(w, w),
    hcat(prim, prim),
    hcat(f, f),
    hcat(f, f),
    zeros(3, 2),
    zeros(16, 2),
    zeros(16, 2),
    hcat(u, u),
    hcat(ω, ω),
    2,
    γ,
    1.0,
    1.0,
    0.5,
    1.0,
    1e-2,
    1.0,
    dx,
    dt,
    hcat(res, res),
    hcat(avg, avg),
)
