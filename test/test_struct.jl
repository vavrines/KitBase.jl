cd(@__DIR__)
D = KitBase.read_dict("config.txt")
for key in keys(D)
    s = key
    @eval $s = $(D[key])
end

#--- settings ---#
Setup(
    "gas",
    "cylinder",
    "2d2f",
    "kfvs",
    "bgk",
    1, # species
    2, # order of accuracy
    "vanleer", # limiter
    "maxwell",
    0.5, # cfl
    10.0, # time
) |> show
KitBase.Scalar(1.0, 1e-3)
KitBase.Radiation(1e-2, 1.0, 0.0, 1e-3, 1000)
gas = KitBase.Gas()
show(gas)
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

fw = (args...) -> [1.0, 0.0, rand()]
KitBase.IB(fw, gas) |> show

#--- control volume ---#
KitBase.ControlVolume1D(w, prim) |> show
KitBase.ControlVolume1D1F(w, prim, h) |> show
KitBase.ControlVolume1D2F(w, prim, h, b) |> show
KitBase.ControlVolume1D3F(
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
KitBase.ControlVolume1D3F(rand(4), rand(5), rand(nu), rand(nu), rand(nu)) |> show
KitBase.ControlVolume1D4F(
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
KitBase.ControlVolume2D(w, prim) |> show
KitBase.ControlVolume2D1F(w, prim, h) |> show
KitBase.ControlVolume2D2F(w, prim, h, b) |> show
KitBase.ControlVolume2D3F(w, prim, h, b, b, zeros(3), zeros(3), zeros(3, 2)) |> show
# Rykov
KitBase.ControlVolume2D3F(rand(4), rand(4), rand(8, 8), rand(8, 8), rand(8, 8)) |> show
# unstructured
KitBase.ControlVolumeUS([1, 0], x, dx, w, prim)
KitBase.ControlVolumeUS1F([1, 0], x, dx, w, prim, h)
KitBase.ControlVolumeUS2F([1, 0], x, dx, w, prim, h, b)

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
sol_w = rand(4, 3)
sol_prim = rand(4, 3)
sol_h = rand(6, 3)
KitBase.Solution1D(sol_w, sol_prim)
KitBase.Solution1D(sol_w, sol_prim, sol_h)
KitBase.Solution1D(sol_w, sol_prim, sol_h, sol_h)

sol_w = rand(4, 3, 2)
sol_prim = rand(4, 3, 2)
sol_h = rand(6, 3, 2)
KitBase.Solution2D(sol_w, sol_prim)
KitBase.Solution2D(sol_w, sol_prim, sol_h)
KitBase.Solution2D(sol_w, sol_prim, sol_h, sol_h)

#--- flux ---#
sol_w = rand(4, 3)
sol_h = rand(6, 3)
KitBase.Flux1D(sol_w)
KitBase.Flux1D(sol_w, sol_h)
KitBase.Flux1D(sol_w, sol_h, sol_h)

sol_w = rand(4, 3, 2)
sol_h = rand(6, 3, 2)
n = [rand(2) for i = 1:3, j = 1:2]
KitBase.Flux2D(n, sol_w)
KitBase.Flux2D(n, sol_w, sol_h)
KitBase.Flux2D(n, sol_w, sol_h, sol_h)
