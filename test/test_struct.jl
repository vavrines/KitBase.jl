cd(@__DIR__)
D = read_dict("config.txt")
for key in keys(D)
    s = Symbol(key)
    @eval $s = $(D[key])
end

#--- settings ---#
Setup()
Gas(knudsen, mach, prandtl, inK, 3.0, omega, alphaRef, omegaRef, 0.01)
Particle(knudsen, mach, prandtl, inK, 3.0, omega, alphaRef, omegaRef, 0.01, 1e-4, 10000)
Mixture([0.1, 0.5], mach, prandtl, inK, 3.0, 1.0, 0.5, 0.5, 0.5)
Plasma1D([0.1, 0.5], mach, prandtl, inK, 3.0, 1.0, 0.5, 0.5, 0.5, 0.01, 0.01, 100.0, 1.0, 1.0)
Plasma2D([0.1, 0.5], mach, prandtl, inK, 3.0, 1.0, 0.5, 0.5, 0.5, 0.01, 0.01, 100.0, 1.0, 1.0)

prim = [1.0, 0.0, 1.0]
w = prim_conserve(prim, 1.4)
u = Float64.(collect(umin:nu:umax)) # u should be float
h = maxwellian(u, prim)
b = h .* inK ./ (2.0 * prim[end])
x = Float64(x0) # x & dx should be of same type
dx = x0 / nx

#--- control volume ---#
ControlVolume1D(x, dx, w, prim)
ControlVolume1D1F(x, dx, w, prim, h)
ControlVolume1D2F(x, dx, w, prim, h, b)
ControlVolume1D3F(
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
)
ControlVolume1D4F(
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
)
ControlVolume2D(x, dx, x, dx, w, prim)
ControlVolume2D1F(x, dx, x, dx, w, prim, h)
ControlVolume2D2F(x, dx, x, dx, w, prim, h, b)
ControlVolume2D3F(x, dx, x, dx, w, prim, h, b, b, zeros(3), zeros(3), zeros(3, 2))

#--- interface ---#
Interface1D(w)
Interface1D1F(w, h)
Interface1D2F(w, h)
Interface1D3F(w, h, zeros(3))
Interface1D4F(w, h, zeros(3))
cosa = 1 / √2
sina = 1 / √2
Interface2D(dx, cosa, sina, w)
Interface2D1F(dx, cosa, sina, w, h)
Interface2D2F(dx, cosa, sina, w, h)

#--- solution ---#
sol_w = [w for i in 1:2]
sol_prim = [prim for i in 1:2]
sol_h = [h for i in 1:2]
Solution1D(sol_w, sol_prim)
Solution1D1F(sol_w, sol_prim, sol_h)
Solution1D2F(sol_w, sol_prim, sol_h, sol_h)

sol_w = [w for i in 1:2, j in 1:2]
sol_prim = [prim for i in 1:2, j in 1:2]
sol_h = [h for i in 1:2, j in 1:2]
Solution2D(sol_w, sol_prim)
Solution2D1F(sol_w, sol_prim, sol_h)
Solution2D2F(sol_w, sol_prim, sol_h, sol_h)

#--- flux ---#
Flux1D(w, w)
Flux1D1F(w, w, h)
Flux1D2F(w, w, h, h)
Flux2D(zeros(2), w, w, zeros(2), w, w)
Flux2D1F(zeros(2), w, w, h, zeros(2), w, w, h)
Flux2D2F(zeros(2), w, w, h, h, zeros(2), w, w, h, h)

#--- particle ---#
Particle1D(1e-4, 0.1, randn(3), 34, 0.2)
Particle2D(1e-4, 0.1, 0.3, randn(3), 34, 21, 0.2)
ControlVolumeParticle1D(x, dx, w, prim)
ControlVolumeParticle2D(x, dx, x, dx, w, prim)