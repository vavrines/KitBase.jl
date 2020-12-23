#--- continuum ---#
prim = [1.0, 0.0, 1.0]
KitBase.prim_conserve(prim, 3.0)
KitBase.prim_conserve(prim[1], prim[2], prim[3], 3.0)

mprim = hcat(prim, prim)
KitBase.mixture_prim_conserve(mprim, 3.0)

KitBase.conserve_prim(prim[1])
KitBase.conserve_prim(prim[1], 1.0)
KitBase.conserve_prim(prim, 3.0)
KitBase.conserve_prim(prim[1], prim[2], prim[3], 3.0)
KitBase.mixture_conserve_prim(mprim, 3.0)

prim = [1.0, 0.2, 0.3, -0.1, 1.0]
mprim = hcat(prim, prim)
KitBase.em_coefficients(mprim, randn(3), randn(3), 100, 0.01, 0.01, 0.001)

KitBase.advection_flux(1.0, -0.1)
KitBase.burgers_flux(1.0)
KitBase.euler_flux(prim, 3.0)
KitBase.euler_jacobi(prim, 3.0)

#--- atom ---#
KitBase.pdf_slope(1.0, 0.1)
KitBase.pdf_slope(prim, randn(5), 0.0)
KitBase.mixture_pdf_slope(mprim, randn(5, 2), 0.0)

u = collect(-5:0.2:5)
ω = ones(51) .* 0.2
prim = [1.0, 0.0, 1.0]
M = KitBase.maxwellian(u, prim)
mprim = hcat(prim, prim)
KitBase.mixture_maxwellian(hcat(u, u), mprim)

KitBase.shakhov(u, M, 0.01, prim, 1.0)
KitBase.shakhov(u, M, M, 0.01, prim, 1.0, 2.0)

KitBase.reduce_distribution(randn(16, 51), ω, 1)
KitBase.reduce_distribution(randn(16, 24, 24), ones(24, 24), 1)
KitBase.full_distribution(M, M, u, ω, ones(51, 24, 24), ones(51, 24, 24), 1.0, 3.0)

KitBase.ref_vhs_vis(1.0, 1.0, 0.5)
KitBase.vhs_collision_time(prim, 1e-3, 0.81)
KitBase.hs_boltz_kn(1e-3, 1.0)
phi, psi, phipsi =
    KitBase.kernel_mode(5, 5.0, 5.0, 5.0, 0.1, 0.1, 0.1, 16, 16, 16, 1.0, quad_num = 16)
KitBase.boltzmann_fft(rand(16, 16, 16), 1.0, 5, phi, psi, phipsi)
KitBase.boltzmann_fft!(rand(16, 16, 16), rand(16, 16, 16), 1.0, 5, phi, psi, phipsi)

τ = KitBase.aap_hs_collision_time(mprim, 1.0, 0.5, 0.5, 0.5, 1.0)
KitBase.aap_hs_prim(mprim, τ, 1.0, 0.5, 0.5, 0.5, 1.0)

KitBase.aap_hs_diffeq!(
    similar(mprim),
    mprim,
    [τ[1], τ[2], 1.0, 0.5, 0.5, 0.5, 1.0, 3.0],
    0.0,
)
KitBase.shift_pdf!(M, 1.0, 1e-4, 1e-4)

#--- particle ---#
KitBase.sample_maxwell([1.0, 0.0, 1.0])
KitBase.next_collision_time(1.0)

ptc = KitBase.Particle1D(1e-3, 0.0, randn(3), 0.1, 1)
KitBase.sample_particle!(ptc, 1e-4, rand(), randn(3), rand(), 2, 0, 0.1)
