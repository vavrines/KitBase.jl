cd(@__DIR__)
D = KitBase.read_dict("config.txt")
for key in keys(D)
    s = key
    @eval $s = $(D[key])
end

#--- settings ---#
KitBase.Setup() |> show
KitBase.Scalar(1.0, 1e-3)
KitBase.Radiation(1e-2, 1.0, 0.0, 1e-3, 1000)
KitBase.Gas(knudsen, mach, prandtl, inK, 3.0, omega, alphaRef, omegaRef, 0.01) |> show
KitBase.Gas(knudsen, mach, prandtl, inK, 3.0, omega, alphaRef, omegaRef, 0.01, 1e-4, 10000) |> show
KitBase.Mixture([0.1, 0.5], mach, prandtl, inK, 3.0, 1.0, 0.5, 0.5, 0.5) |> show
KitBase.Plasma1D(
    [0.1, 0.5],
    mach,
    prandtl,
    inK,
    3.0,
    1.0,
    0.5,
    0.5,
    0.5,
    0.01,
    0.01,
    100.0,
    1.0,
    1.0,
) |> show
KitBase.Plasma2D(
    [0.1, 0.5],
    mach,
    prandtl,
    inK,
    3.0,
    1.0,
    0.5,
    0.5,
    0.5,
    0.01,
    0.01,
    100.0,
    1.0,
    1.0,
) |> show
KitBase.DiatomicGas(
    knudsen,
    mach,
    prandtl,
    inK,
    inK,
    1.4,
    0.81,
    1.0,
    0.5,
    1e-3,
    89.1 / 273,
    18.1,
    1 / 1.55,
    0.2354,
    0.3049,
) |> show

begin
    prim = [1.0, 0.0, 1.0]
    w = KitBase.prim_conserve(prim, 1.4)
    u = Float64.(collect(umin:nu:umax)) # u should be float
    h = KitBase.maxwellian(u, prim)
    b = h .* inK ./ (2.0 * prim[end])
    x = Float64(x0) # x & dx should be of same type
    dx = x0 / nx
end

KitBase.IB(w, prim, prim, w, prim, prim) |> show
KitBase.IB(
    hcat(w, w),
    hcat(prim, prim),
    hcat(prim, prim),
    hcat(w, w),
    hcat(prim, prim),
    hcat(prim, prim),
) |> show
KitBase.IB1F(w, prim, h, prim, w, prim, h, prim) |> show
KitBase.IB2F(w, prim, h, h, prim, w, prim, h, h, prim) |> show
KitBase.IB3F(
    w,
    prim,
    h,
    h,
    h,
    prim,
    zeros(3),
    zeros(3),
    zeros(3, 2),
    w,
    prim,
    h,
    h,
    h,
    prim,
    zeros(3),
    zeros(3),
    zeros(3, 2),
) |> show
KitBase.IB4F(
    w,
    prim,
    h,
    h,
    h,
    h,
    prim,
    zeros(3),
    zeros(3),
    zeros(3, 2),
    w,
    prim,
    h,
    h,
    h,
    h,
    prim,
    zeros(3),
    zeros(3),
    zeros(3, 2),
) |> show

#--- control volume ---#
KitBase.ControlVolume1D(x, dx, w, prim) |> show
KitBase.ControlVolume1D1F(x, dx, w, prim, h) |> show
KitBase.ControlVolume1D2F(x, dx, w, prim, h, b) |> show
KitBase.ControlVolume1D3F(
    x,
    dx,
    hcat(w, w),
    hcat(prim, prim),
    zeros(nu, nu, 2),
    zeros(nu, nu, 2),
    zeros(nu, nu, 2),
    zeros(3),
    zeros(3),
    zeros(3, 2),
) |> show
KitBase.ControlVolume1D3F(
    x,
    dx,
    zeros(5, 7, 2), # indexed with [flow entry, uq, species]
    zeros(5, 7, 2),
    zeros(nu, nu, 7, 2),
    zeros(nu, nu, 7, 2),
    zeros(nu, nu, 7, 2),
    zeros(3, 7),
    zeros(3, 7),
    zeros(3, 7, 2),
) |> show
# Rykov
KitBase.ControlVolume1D3F(x, dx, rand(4), rand(5), rand(nu), rand(nu), rand(nu)) |> show
KitBase.ControlVolume1D4F(
    x,
    dx,
    hcat(w, w),
    hcat(prim, prim),
    hcat(h, h),
    hcat(h, h),
    hcat(h, h),
    hcat(h, h),
    zeros(3),
    zeros(3),
    zeros(3, 2),
) |> show
KitBase.ControlVolume1D4F(
    x,
    dx,
    zeros(5, 7, 2), # indexed with [flow entry, uq, species]
    zeros(5, 7, 2),
    zeros(nu, 7, 2),
    zeros(nu, 7, 2),
    zeros(nu, 7, 2),
    zeros(nu, 7, 2),
    zeros(3, 7),
    zeros(3, 7),
    zeros(3, 7, 2),
) |> show
KitBase.ControlVolume2D(x, dx, x, dx, w, prim) |> show
KitBase.ControlVolume2D1F(x, dx, x, dx, w, prim, h) |> show
KitBase.ControlVolume2D2F(x, dx, x, dx, w, prim, h, b) |> show
KitBase.ControlVolume2D3F(x, dx, x, dx, w, prim, h, b, b, zeros(3), zeros(3), zeros(3, 2)) |> show
# Rykov
KitBase.ControlVolume2D3F(
    x,
    dx,
    x,
    dx,
    rand(4),
    rand(4),
    rand(8, 8),
    rand(8, 8),
    rand(8, 8),
) |> show
# unstructured
KitBase.ControlVolumeUS([1, 0], x, dx, w, prim)
KitBase.ControlVolumeUS1F([1,0], x, dx, w, prim, h)
KitBase.ControlVolumeUS2F([1,0], x, dx, w, prim, h, b)

#--- interface ---#
KitBase.Interface1D(w) |> show
KitBase.Interface1D1F(w, h) |> show
KitBase.Interface1D2F(w, h) |> show
KitBase.Interface1D3F(w, h, zeros(3)) |> show
KitBase.Interface1D3F(zeros(5, 7, 2), zeros(21, 21, 7, 2), zeros(3, 7)) |> show
KitBase.Interface1D3F(w, h) |> show # Rykov
KitBase.Interface1D4F(w, h, zeros(3)) |> show
KitBase.Interface1D4F(zeros(5, 7, 2), zeros(21, 7, 2), zeros(3, 7)) |> show
cosa = 1 / √2
sina = 1 / √2
KitBase.Interface2D(dx, cosa, sina, w) |> show
KitBase.Interface2D1F(dx, cosa, sina, w, h) |> show
KitBase.Interface2D2F(dx, cosa, sina, w, h) |> show
KitBase.Interface2D3F(dx, cosa, sina, w, h, zeros(3)) |> show
KitBase.Interface2D3F(dx, cosa, sina, w, h, zeros(3, 7)) |> show
KitBase.Interface2D3F(dx, cosa, sina, w, h) |> show

#--- solution ---#
sol_w = [w for i = 1:2]
sol_prim = [prim for i = 1:2]
sol_h = [h for i = 1:2]
KitBase.Solution1D(sol_w, sol_prim)
KitBase.Solution1D1F(sol_w, sol_prim, sol_h)
KitBase.Solution1D1F(sol_w, sol_prim, similar(sol_w), sol_h, similar(sol_h))
KitBase.Solution1D2F(sol_w, sol_prim, sol_h, sol_h)
KitBase.Solution1D2F(
    sol_w,
    sol_prim,
    similar(sol_w),
    sol_h,
    sol_h,
    similar(sol_h),
    similar(sol_h),
)

sol_w = [w for i = 1:2, j = 1:2]
sol_prim = [prim for i = 1:2, j = 1:2]
sol_h = [h for i = 1:2, j = 1:2]
KitBase.Solution2D(sol_w, sol_prim)
KitBase.Solution2D1F(sol_w, sol_prim, sol_h)
KitBase.Solution2D1F(sol_w, sol_prim, similar(sol_w), sol_h, similar(sol_h))
KitBase.Solution2D2F(sol_w, sol_prim, sol_h, sol_h)
KitBase.Solution2D2F(
    sol_w,
    sol_prim,
    similar(sol_w),
    sol_h,
    sol_h,
    similar(sol_h),
    similar(sol_h),
)

#--- flux ---#
KitBase.Flux1D(w, w)
KitBase.Flux1D1F(w, w, h)
KitBase.Flux1D2F(w, w, h, h)
KitBase.Flux2D(zeros(2), w, w, zeros(2), w, w)
KitBase.Flux2D1F(zeros(2), w, w, h, zeros(2), w, w, h)
KitBase.Flux2D2F(zeros(2), w, w, h, h, zeros(2), w, w, h, h)
