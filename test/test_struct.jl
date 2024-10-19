cd(@__DIR__)
D = KB.read_dict("config.txt")
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
    false,
) |> show
config_ntuple(; u0=-8, c=1)

KB.Scalar(1.0, 1e-3)
KB.Radiation(1e-2, 1.0, 0.0, 1e-3, 1000)
gas = KB.Gas()
show(gas)
KB.Mixture([0.1, 0.5], mach, prandtl, inK, 3.0, 1.0, 0.5, 0.5, 0.5) |> show
KB.Plasma1D(
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
KB.Plasma2D(
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
KB.DiatomicGas(
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
    w = KB.prim_conserve(prim, 1.4)
    u = Float64.(collect(umin:nu:umax)) # u should be float
    h = KB.maxwellian(u, prim)
    b = h .* inK ./ (2.0 * prim[end])
    x = Float64(x0) # x & dx should be of same type
    dx = (x1 - x0) / nx
end

#fw = (args...) -> [1.0, 0.0, rand()]
#KB.IB(fw, gas) |> show

#--- control volume ---#
KB.ControlVolume(w, prim, 1)
KB.ControlVolume(w, prim, 2)
KB.ControlVolume(w, prim, h, 1)
KB.ControlVolume(w, prim, h, 2)
KB.ControlVolume(w, prim, h, h, 1)
KB.ControlVolume(w, prim, h, h, 2)

KB.ControlVolume1D(w, prim) |> show
KB.ControlVolume1D1F(w, prim, h) |> show
KB.ControlVolume1D2F(w, prim, h, b) |> show
KB.ControlVolume1D3F(
    hcat(w, w),
    hcat(prim, prim),
    zeros(nu, nu, 2),
    zeros(nu, nu, 2),
    zeros(nu, nu, 2),
    zeros(3),
    zeros(3),
    zeros(3, 2),
) |> show
KB.ControlVolume1D3F(
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
KB.ControlVolume1D3F(rand(4), rand(5), rand(nu), rand(nu), rand(nu)) |> show
KB.ControlVolume1D4F(
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
KB.ControlVolume1D4F(
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
KB.ControlVolume2D(w, prim) |> show
KB.ControlVolume2D1F(w, prim, h) |> show
KB.ControlVolume2D2F(w, prim, h, b) |> show
KB.ControlVolume2D3F(w, prim, h, b, b, zeros(3), zeros(3), zeros(3, 2)) |> show
# Rykov
KB.ControlVolume2D3F(rand(4), rand(4), rand(8, 8), rand(8, 8), rand(8, 8)) |> show
# unstructured
KB.ControlVolumeUS([1, 0], x, dx, w, prim)
KB.ControlVolumeUS1F([1, 0], x, dx, w, prim, h)
KB.ControlVolumeUS2F([1, 0], x, dx, w, prim, h, b)

#--- interface ---#
KB.Interface1D(w) |> show
KB.Interface1D1F(w, h) |> show
KB.Interface1D2F(w, h) |> show
KB.Interface1D3F(w, h, zeros(3)) |> show
KB.Interface1D3F(zeros(5, 7, 2), zeros(21, 21, 7, 2), zeros(3, 7)) |> show
KB.Interface1D3F(w, h) |> show # Rykov
KB.Interface1D4F(w, h, zeros(3)) |> show
KB.Interface1D4F(zeros(5, 7, 2), zeros(21, 7, 2), zeros(3, 7)) |> show
cosa = 1 / √2
sina = 1 / √2
KB.Interface2D(dx, cosa, sina, w) |> show
KB.Interface2D1F(dx, cosa, sina, w, h) |> show
KB.Interface2D2F(dx, cosa, sina, w, h) |> show
KB.Interface2D3F(dx, cosa, sina, w, h, zeros(3)) |> show
KB.Interface2D3F(dx, cosa, sina, w, h, zeros(3, 7)) |> show
KB.Interface2D3F(dx, cosa, sina, w, h) |> show

#--- solution ---#
sol_w = rand(4, 3)
sol_prim = rand(4, 3)
sol_h = rand(6, 3)

x = KB.Solution1D(sol_w, sol_prim)
KB.Solution1D(sol_w, sol_prim, sol_h)
KB.Solution1D(sol_w, sol_prim, sol_h, sol_h)

sol_w = rand(4, 3, 2)
sol_prim = rand(4, 3, 2)
sol_h = rand(6, 3, 2)

sol_w = [rand(4) for i in 1:4, j in 1:3]
sol_prim = [rand(4) for i in 1:4, j in 1:3]
sol_h = [rand(16, 16) for i in 1:4, j in 1:3]

KB.Solution2D(sol_w, sol_prim)
KB.Solution2D(sol_w, sol_prim, sol_h)
KB.Solution2D(sol_w, sol_prim, sol_h, sol_h)

#--- flux ---#
sol_w = rand(4, 3)
sol_h = rand(6, 3)
KB.Flux1D(sol_w)
KB.Flux1D(sol_w, sol_h)
KB.Flux1D(sol_w, sol_h, sol_h)

sol_w = rand(4, 3, 2)
sol_h = rand(6, 3, 2)
n = [rand(2) for i in 1:3, j in 1:2]
KB.Flux2D(n, sol_w)
KB.Flux2D(n, sol_w, sol_h)
KB.Flux2D(n, sol_w, sol_h, sol_h)
