begin
    fwL = zeros(3)
    fwR = zeros(3)
    γ = 5 / 3
    prim = [1.0, 0.0, 1.0]
    w = KitBase.prim_conserve(prim, γ)
    dx = 1e-2
    dt = 1e-3
    res = zeros(3)
    avg = zeros(3)
end

KitBase.step!(fwL, w, prim, fwR, γ, dx, res, avg)

begin
    ffL = zeros(16)
    ffR = zeros(16)
    f = rand(16)
    u = randn(16)
    ω = ones(16)
end

KitBase.step!(fwL, ffL, w, prim, f, fwR, ffR, u, ω, γ, 1e-3, 0.72, 1.0, dx, dt, res, avg)

KitBase.step!(
    KitBase.prim_conserve([1.0, 0.0, 0.0, 0.0, 1.0], γ),
    [1.0, 0.0, 0.0, 0.0, 1.0],
    rand(16, 16, 16),
    zeros(5),
    zeros(16, 16, 16),
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
    w,
    prim,
    f,
    f,
    fwL,
    ffL,
    ffL,
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
    hcat(w, w),
    hcat(prim, prim),
    hcat(f, f),
    hcat(f, f),
    zeros(3, 2),
    zeros(16, 2),
    zeros(16, 2),
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

# fsm
KitBase.step!(
    [1.0, 0.0, 0.0, 0.0, 1.0],
    KitBase.conserve_prim([1.0, 0.0, 0.0, 0.0, 1.0], 5 / 3),
    rand(16, 16, 16),
    zeros(5),
    zeros(16, 16, 16),
    zeros(5),
    zeros(16, 16, 16),
    5 / 3,
    1.0,
    5,
    ones(16, 16, 16, 20),
    ones(16, 16, 16, 20),
    ones(16, 16, 16),
    0.1,
    0.01,
    zeros(5),
    zeros(5),
    :fsm,
)
KitBase.step!(
    [1.0, 0.0, 0.0, 0.0, 1.0],
    KitBase.conserve_prim([1.0, 0.0, 0.0, 0.0, 1.0], 5 / 3),
    rand(16, 16, 16),
    zeros(5),
    zeros(16, 16, 16),
    zeros(5),
    zeros(16, 16, 16),
    zeros(5),
    zeros(16, 16, 16),
    zeros(5),
    zeros(16, 16, 16),
    5 / 3,
    1.0,
    5,
    ones(16, 16, 16, 20),
    ones(16, 16, 16, 20),
    ones(16, 16, 16),
    0.1,
    0.01,
    zeros(5),
    zeros(5),
    :fsm,
)

# Rykov
KitBase.step!(
    [1.0, 0.0, 1.0, 0.1],
    KitBase.conserve_prim([1.0, 0.0, 1.0, 0.1], 2, 2),
    f,
    f,
    f,
    zeros(4),
    ffL,
    ffL,
    ffL,
    zeros(4),
    ffR,
    ffR,
    ffR,
    u,
    ω,
    2,
    2,
    1e-3,
    0.81,
    0.72,
    91.5 / 273,
    18.1,
    1 / 1.55,
    0.2354,
    0.3049,
    dx,
    dt,
    zeros(4),
    zeros(4),
    :rykov,
)
