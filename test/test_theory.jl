#--- continuum ---#
prim = [1.0, 0.0, 1.0]
prim_conserve(prim, 3.0)
prim_conserve(prim[1], prim[2], prim[3], 3.0)

mprim = hcat(prim, prim)
mixture_prim_conserve(mprim, 3.0)

conserve_prim(prim[1])
conserve_prim(prim[1], 1.0)
conserve_prim(prim, 3.0)
conserve_prim(prim[1], prim[2], prim[3], 3.0)
mixture_conserve_prim(mprim, 3.0)

prim = [1.0, 0.2, 0.3, -0.1, 1.0]
mprim = hcat(prim, prim)
em_coefficients(mprim, randn(3), randn(3), 100, 0.01, 0.01, 0.001)

advection_flux(1.0, -0.1)
burgers_flux(1.0)
euler_flux(prim, 3.0)
euler_jacobi(prim, 3.0)

#--- atom ---#
pdf_slope(1.0, 0.1)
pdf_slope(prim, randn(5), 0.0)
mixture_pdf_slope(mprim, randn(5, 2), 0.0)

u = collect(-5:0.2:5)
ω = ones(51) .* 0.2
prim = [1.0, 0.0, 1.0]
M = maxwellian(u, prim)
mprim = hcat(prim, prim)
mixture_maxwellian(hcat(u, u), mprim)

shakhov(u, M, 0.01, prim, 1.0)
shakhov(u, M, M, 0.01, prim, 1.0, 2.0)

reduce_distribution(randn(16, 51), ω, 1)
reduce_distribution(randn(16, 24, 24), ones(24, 24), 1)
full_distribution(M, M, u, ω, ones(51, 24, 24), ones(51, 24, 24), 1.0, 3.0)

ref_vhs_vis(1.0, 1.0, 0.5)
vhs_collision_time(prim, 1e-3, 0.81)
hs_boltz_kn(1e-3, 1.0)
kernel_mode(5, 5.0, 5.0, 5.0, 0.1, 0.1, 0.1, 16, 16, 16, 1.0, quad_num=16)
 
