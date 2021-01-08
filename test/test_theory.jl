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
KitBase.maxwellian(randn(16, 16), randn(16, 16), [1.0, 0.0, 0.0, 1.0])
KitBase.maxwellian(
    randn(8, 8, 8),
    randn(8, 8, 8),
    randn(8, 8, 8),
    [1.0, 0.0, 0.0, 0.0, 1.0],
)

KitBase.maxwellian!(rand(16), rand(16), prim)
KitBase.maxwellian!(rand(16, 16), randn(16, 16), randn(16, 16), [1.0, 0.0, 0.0, 1.0])
KitBase.maxwellian!(
    rand(8, 8, 8),
    randn(8, 8, 8),
    randn(8, 8, 8),
    randn(8, 8, 8),
    [1.0, 0.0, 0.0, 0.0, 1.0],
)

mprim = hcat(prim, prim)
KitBase.mixture_maxwellian(hcat(u, u), mprim)
KitBase.mixture_maxwellian(randn(8, 8, 2), randn(8, 8, 2), rand(4, 2))
KitBase.mixture_maxwellian(
    randn(8, 8, 8, 2),
    randn(8, 8, 8, 2),
    randn(8, 8, 8, 2),
    rand(5, 2),
)

KitBase.mixture_maxwellian!(randn(8, 2), randn(8, 2), rand(3, 2))
KitBase.mixture_maxwellian!(randn(8, 8, 2), randn(8, 8, 2), randn(8, 8, 2), rand(4, 2))
KitBase.mixture_maxwellian!(
    randn(8, 8, 8, 2),
    randn(8, 8, 8, 2),
    randn(8, 8, 8, 2),
    randn(8, 8, 8, 2),
    rand(5, 2),
)

KitBase.shakhov(u, M, 0.01, prim, 1.0)
KitBase.shakhov(u, M, M, 0.01, prim, 1.0, 2.0)
KitBase.shakhov(randn(16, 16), randn(16, 16), rand(16, 16), rand(2), [1.0, 0.0, 1.0], 1.0)
KitBase.shakhov(
    randn(16, 16),
    randn(16, 16),
    rand(16, 16),
    rand(16, 16),
    rand(2),
    [1.0, 0.0, 0.0, 1.0],
    1.0,
    1.0,
)
KitBase.shakhov(
    randn(8, 8, 8),
    randn(8, 8, 8),
    rand(8, 8, 8),
    rand(8, 8, 8),
    rand(3),
    [1.0, 0.0, 0.0, 0.0, 1.0],
    1.0,
)

KitBase.shakhov!(randn(16), randn(16), rand(16), 0.01, prim, 1.0)
KitBase.shakhov!(randn(16), randn(16), randn(16), rand(16), rand(16), 0.01, prim, 1.0, 2.0)
KitBase.shakhov!(
    randn(16, 16),
    randn(16, 16),
    randn(16, 16),
    rand(16, 16),
    rand(2),
    [1.0, 0.0, 1.0],
    1.0,
)
KitBase.shakhov!(
    randn(16, 16),
    randn(16, 16),
    randn(16, 16),
    randn(16, 16),
    rand(16, 16),
    rand(16, 16),
    rand(2),
    [1.0, 0.0, 0.0, 1.0],
    1.0,
    1.0,
)
KitBase.shakhov!(
    randn(8, 8, 8),
    randn(8, 8, 8),
    randn(8, 8, 8),
    rand(8, 8, 8),
    rand(8, 8, 8),
    rand(3),
    [1.0, 0.0, 0.0, 0.0, 1.0],
    1.0,
)

# Rykov
KitBase.maxwellian!(
    zeros(16),
    zeros(16),
    zeros(16),
    zeros(16),
    zeros(16),
    zeros(16),
    randn(16),
    rand(5),
    4,
    2,
)
KitBase.maxwellian!(
    zeros(8, 8),
    zeros(8, 8),
    zeros(8, 8),
    zeros(8, 8),
    zeros(8, 8),
    zeros(8, 8),
    randn(8, 8),
    randn(8, 8),
    rand(6),
    4,
    2,
)

KitBase.rykov!(
    zeros(16),
    zeros(16),
    zeros(16),
    zeros(16),
    zeros(16),
    zeros(16),
    randn(16),
    rand(16),
    rand(16),
    rand(16),
    rand(16),
    rand(16),
    rand(16),
    rand(2),
    rand(5),
    0.72,
    4,
    1/1.55,
    0.2354,
    0.3049,
)
KitBase.rykov!(
    zeros(8, 8),
    zeros(8, 8),
    zeros(8, 8),
    zeros(8, 8),
    zeros(8, 8),
    zeros(8, 8),
    randn(8, 8),
    randn(8, 8),
    rand(8, 8),
    rand(8, 8),
    rand(8, 8),
    rand(8, 8),
    rand(8, 8),
    rand(8, 8),
    rand(4),
    rand(6),
    0.72,
    4,
    1/1.55,
    0.2354,
    0.3049,
)

KitBase.reduce_distribution(randn(16, 51), ω, 1)
KitBase.reduce_distribution(randn(16, 24, 24), ones(24, 24), 1)
KitBase.reduce_distribution(
    randn(16, 24, 24),
    randn(16, 24, 24),
    randn(16, 24, 24),
    ones(24, 24),
    1,
)
KitBase.full_distribution(M, M, u, ω, ones(51, 24, 24), ones(51, 24, 24), 1.0, 3.0)
KitBase.full_distribution(M, M, u, ω, ones(51, 24, 24), ones(51, 24, 24), prim, 3.0)

KitBase.ref_vhs_vis(1.0, 1.0, 0.5)
KitBase.vhs_collision_time(prim, 1e-3, 0.81)
KitBase.hs_boltz_kn(1e-3, 1.0)
phi, psi, phipsi =
    KitBase.kernel_mode(5, 5.0, 5.0, 5.0, 0.1, 0.1, 0.1, 16, 16, 16, 1.0, quad_num = 16)
KitBase.boltzmann_fft(rand(16, 16, 16), 1.0, 5, phi, psi, phipsi)
KitBase.boltzmann_fft!(rand(16, 16, 16), rand(16, 16, 16), 1.0, 5, phi, psi, phipsi)

KitBase.boltzmann_ode!(
    zeros(16, 16, 16),
    rand(16, 16, 16),
    (1.0, 5, phi, psi, phipsi),
    0.0,
)
KitBase.bgk_ode!(
    zeros(16, 16, 16),
    rand(16, 16, 16),
    (rand(16, 16, 16), 1e-2),
    0.0,
)

τ = KitBase.aap_hs_collision_time(mprim, 1.0, 0.5, 0.5, 0.5, 1.0)
KitBase.aap_hs_prim(mprim, τ, 1.0, 0.5, 0.5, 0.5, 1.0)
KitBase.aap_hs_prim(rand(4, 2), rand(2), 1.0, 0.5, 0.5, 0.5, 1e-2)
KitBase.aap_hs_prim(rand(5, 2), rand(2), 1.0, 0.5, 0.5, 0.5, 1e-2)

KitBase.aap_hs_diffeq!(
    similar(mprim),
    mprim,
    [τ[1], τ[2], 1.0, 0.5, 0.5, 0.5, 1.0, 3.0],
    0.0,
)
KitBase.shift_pdf!(M, 1.0, 1e-4, 1e-4)
KitBase.shift_pdf!(rand(16, 2), randn(2), rand(2), 1e-4)
