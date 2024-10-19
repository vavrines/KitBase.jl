#--- continuum ---#
prim = [1.0, 0.0, 1.0]
KB.prim_conserve(prim, 3.0)
KB.prim_conserve(prim[1], prim[2], prim[3], 3.0)

mprim = hcat(prim, prim)
KB.mixture_prim_conserve(mprim, 3.0)

KB.conserve_prim(prim[1])
KB.conserve_prim(prim[1], 1.0)
KB.conserve_prim(prim, 3.0)
KB.conserve_prim(prim[1], prim[2], prim[3], 3.0)
KB.mixture_conserve_prim(mprim, 3.0)

# polyatomic
KB.prim_conserve([1.0, 0.0, 1.0, 1.0, 1.0], 5 / 3, 2)
KB.prim_conserve([1.0, 0.0, 0.0, 1.0, 1.0, 1.0], 5 / 3, 2)
KB.conserve_prim([1.0, 0.0, 1.0, 0.1], 5 / 3, 2)
KB.conserve_prim([1.0, 0.0, 0.0, 1.0, 0.1], 5 / 3, 2)
# multi-species polyatomic
KB.mixture_prim_conserve(rand(5, 2), 5 / 3, 2)
KB.mixture_prim_conserve(rand(6, 2), 5 / 3, 2)
KB.mixture_prim_conserve(rand(7, 2), 5 / 3, 2)
KB.mixture_conserve_prim(rand(4, 2), 2, 2)
KB.mixture_conserve_prim(rand(5, 2), 1, 2)
KB.mixture_conserve_prim(rand(6, 2), 0, 2)

prim = [1.0, 0.2, 0.3, -0.1, 1.0]
mprim = hcat(prim, prim)
KB.em_coefficients(mprim, randn(3), randn(3), 100, 0.01, 0.01, 0.001)

KB.advection_flux(1.0, -0.1)
KB.burgers_flux(1.0)
KB.euler_flux(prim, 3.0)
KB.euler_jacobi(prim, 3.0)

#--- thermo ---#
KB.heat_capacity_ratio(2.0, 1)
KB.heat_capacity_ratio(2.0, 2)
KB.heat_capacity_ratio(2.0, 3)
KB.heat_capacity_ratio(2.0, 2, 1)
KB.heat_capacity_ratio(2.0, 2, 2)
KB.heat_capacity_ratio(2.0, 2, 3)

KB.sound_speed(1.0, 5 / 3)
KB.sound_speed([1.0, 0.0, 1.0], 5 / 3)
KB.sound_speed(rand(3, 2), 5 / 3)

#--- atom ---#
KB.pdf_slope(1.0, 0.1)
KB.pdf_slope(prim, randn(5), 0.0)
KB.mixture_pdf_slope(mprim, randn(5, 2), 0.0)

u = collect(-5:0.2:5)
ω = ones(51) .* 0.2
prim = [1.0, 0.0, 1.0]

M = KB.maxwellian(u, prim)
KB.maxwellian(randn(16, 16), randn(16, 16), [1.0, 0.0, 0.0, 1.0])
KB.maxwellian(randn(8, 8, 8), randn(8, 8, 8), randn(8, 8, 8), [1.0, 0.0, 0.0, 0.0, 1.0])

KB.energy_maxwellian(M, prim, 2)
KB.mixture_energy_maxwellian(hcat(M, M), hcat(prim, prim), 2)
KB.mixture_energy_maxwellian(rand(16, 16, 2), hcat(prim, prim), 2)

KB.maxwellian!(rand(16), rand(16), prim)
KB.maxwellian!(rand(16, 16), randn(16, 16), randn(16, 16), [1.0, 0.0, 0.0, 1.0])
KB.maxwellian!(
    rand(8, 8, 8),
    randn(8, 8, 8),
    randn(8, 8, 8),
    randn(8, 8, 8),
    [1.0, 0.0, 0.0, 0.0, 1.0],
)

mprim = hcat(prim, prim)
KB.mixture_maxwellian(hcat(u, u), mprim)
KB.mixture_maxwellian(randn(8, 8, 2), randn(8, 8, 2), rand(4, 2))
KB.mixture_maxwellian(randn(8, 8, 8, 2), randn(8, 8, 8, 2), randn(8, 8, 8, 2), rand(5, 2))

KB.mixture_maxwellian!(randn(8, 2), randn(8, 2), rand(3, 2))
KB.mixture_maxwellian!(randn(8, 8, 2), randn(8, 8, 2), randn(8, 8, 2), rand(4, 2))
KB.mixture_maxwellian!(
    randn(8, 8, 8, 2),
    randn(8, 8, 8, 2),
    randn(8, 8, 8, 2),
    randn(8, 8, 8, 2),
    rand(5, 2),
)

KB.shakhov(u, M, 0.01, prim, 1.0)
KB.shakhov(u, M, M, 0.01, prim, 1.0, 2.0)
KB.shakhov(randn(16, 16), randn(16, 16), rand(16, 16), rand(2), [1.0, 0.0, 1.0], 1.0)
KB.shakhov(
    randn(16, 16),
    randn(16, 16),
    rand(16, 16),
    rand(16, 16),
    rand(2),
    [1.0, 0.0, 0.0, 1.0],
    1.0,
    1.0,
)
KB.shakhov(
    randn(8, 8, 8),
    randn(8, 8, 8),
    rand(8, 8, 8),
    rand(8, 8, 8),
    rand(3),
    [1.0, 0.0, 0.0, 0.0, 1.0],
    1.0,
)

KB.shakhov!(randn(16), randn(16), rand(16), 0.01, prim, 1.0)
KB.shakhov!(randn(16), randn(16), randn(16), rand(16), rand(16), 0.01, prim, 1.0, 2.0)
KB.shakhov!(
    randn(16, 16),
    randn(16, 16),
    randn(16, 16),
    rand(16, 16),
    rand(2),
    [1.0, 0.0, 1.0],
    1.0,
)
KB.shakhov!(
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
KB.shakhov!(
    randn(8, 8, 8),
    randn(8, 8, 8),
    randn(8, 8, 8),
    rand(8, 8, 8),
    rand(8, 8, 8),
    rand(3),
    [1.0, 0.0, 0.0, 0.0, 1.0],
    1.0,
)

# ES-BGK
vs = VSpace1D()
prim = [1.0, 0, 1]
f = maxwellian(vs.u, prim)
KB.esbgk_ode!(zero(f), f, (vs.u, vs.weights, prim, 2 / 3, 1), 0.0)

vs = VSpace2D()
prim = [1.0, 0, 0, 1]
f = maxwellian(vs.u, vs.v, prim)
KB.esbgk_ode!(zero(f), f, (vs.u, vs.v, vs.weights, prim, 2 / 3, 1), 0.0)

vs = VSpace3D()
prim = [1.0, 0, 0, 0, 1]
f = maxwellian(vs.u, vs.v, vs.w, prim)
KB.esbgk_ode!(zero(f), f, (vs.u, vs.v, vs.w, vs.weights, prim, 2 / 3, 1), 0.0)

# Rykov
KB.polyatomic_maxwellian!(
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
KB.polyatomic_maxwellian!(
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

KB.rykov!(
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
    1 / 1.55,
    0.2354,
    0.3049,
)
KB.rykov!(
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
    1 / 1.55,
    0.2354,
    0.3049,
)

# BIP
KB.polyatomic_maxwellian!(zeros(16), zeros(16), zeros(16), randn(16), rand(5), 2, 2)
KB.polyatomic_maxwellian!(
    zeros(16),
    zeros(16),
    zeros(16),
    randn(16),
    randn(16),
    rand(6),
    1,
    2,
)
KB.polyatomic_maxwellian!(
    zeros(16),
    zeros(16),
    zeros(16),
    randn(16),
    randn(16),
    randn(16),
    rand(7),
    1,
    2,
)

# multi-species BIP
KB.mixture_polyatomic_maxwellian!(
    zeros(16, 2),
    zeros(16, 2),
    zeros(16, 2),
    randn(16, 2),
    rand(5, 2),
    2,
    2,
)
KB.mixture_polyatomic_maxwellian!(
    zeros(16, 2),
    zeros(16, 2),
    zeros(16, 2),
    randn(16, 2),
    randn(16, 2),
    rand(6, 2),
    1,
    2,
)
KB.mixture_polyatomic_maxwellian!(
    zeros(16, 2),
    zeros(16, 2),
    zeros(16, 2),
    randn(16, 2),
    randn(16, 2),
    randn(16, 2),
    rand(7, 2),
    0,
    2,
)

KB.f_maxwellian(rand(16, 2))
KB.f_maxwellian(rand(16, 2), rand(16, 2))

KB.reduce_distribution(randn(16, 51), ω, 1)
KB.reduce_distribution(randn(16, 24, 24), ones(24, 24), 1)
KB.reduce_distribution(
    randn(16, 24, 24),
    randn(16, 24, 24),
    randn(16, 24, 24),
    ones(24, 24),
    1,
)
KB.full_distribution(M, M, u, ω, ones(51, 24, 24), ones(51, 24, 24), 1.0, 3.0)
KB.full_distribution(M, M, u, ω, ones(51, 24, 24), ones(51, 24, 24), prim, 3.0)

KB.ref_vhs_vis(1.0, 1.0, 0.5)
KB.vhs_collision_time(prim, 1e-3, 0.81)

KB.νbgk_relaxation_time(0.1, 1, rand(3), Class{1})
KB.νbgk_relaxation_time(0.1, 1, rand(3), Class{2})
KB.νbgk_relaxation_time(0.1, 1, 1, rand(4), Class{3})
KB.νbgk_relaxation_time(0.1, 1, 1, 1, rand(5), Class{4})

KB.νshakhov_relaxation_time(0.1, 1, rand(3))

KB.rykov_zr(100, 91.5, 18.1)

KB.hs_boltz_kn(1e-3, 1.0)
vs = VSpace3D(-5, 5, 16, -5, 5, 16, -5, 5, 16)
fsm = KB.fsm_kernel(vs, 1e-3)
phi, psi, phipsi =
    KB.kernel_mode(5, 5.0, 5.0, 5.0, 0.1, 0.1, 0.1, 16, 16, 16, 1.0; quad_num=16)
KB.kernel_mode(5, 5.0, 5.0, 5.0, 16, 16, 16, 1.0; quad_num=16)
KB.kernel_mode(5, 5.0, 5.0, 0.1, 0.1, 16, 16; quad_num=16)
KB.boltzmann_fft(rand(16, 16, 16), fsm)
KB.boltzmann_fft!(rand(16, 16, 16), rand(16, 16, 16), fsm)

KB.boltzmann_ode!(zeros(16, 16, 16), rand(16, 16, 16), (1.0, 5, phi, psi, phipsi), 0.0)
KB.bgk_ode!(zeros(16, 16, 16), rand(16, 16, 16), (rand(16, 16, 16), 1e-2), 0.0)

vs = KB.VSpace3D(-5.0, 5.0, 16, -5.0, 5.0, 16, -5.0, 5.0, 16; type="algebra")
u, v, w = vs.u[:, 1, 1], vs.v[1, :, 1], vs.w[1, 1, :]
vnu = hcat(vs.u[:], vs.v[:], vs.w[:])
uuni1d = linspace(vs.u[1, 1, 1], vs.u[end, 1, 1], 16)
vuni1d = linspace(vs.v[1, 1, 1], vs.v[1, end, 1], 16)
wuni1d = linspace(vs.w[1, 1, 1], vs.w[1, 1, end], 16)
u13d = [uuni1d[i] for i in 1:16, j in 1:16, k in 1:16]
v13d = [vuni1d[j] for i in 1:16, j in 1:16, k in 1:16]
w13d = [wuni1d[k] for i in 1:16, j in 1:16, k in 1:16]
vuni = hcat(u13d[:], v13d[:], w13d[:])

τ = KB.aap_hs_collision_time(mprim, 1.0, 0.5, 0.5, 0.5, 1.0)
KB.aap_hs_prim(mprim, τ, 1.0, 0.5, 0.5, 0.5, 1.0)
KB.aap_hs_prim(rand(4, 2), rand(2), 1.0, 0.5, 0.5, 0.5, 1e-2)
KB.aap_hs_prim(rand(5, 2), rand(2), 1.0, 0.5, 0.5, 0.5, 1e-2)

KB.aap_hs_diffeq!(similar(mprim), mprim, [τ[1], τ[2], 1.0, 0.5, 0.5, 0.5, 1.0, 3.0], 0.0)
KB.shift_pdf!(M, 1.0, 1e-4, 1e-4)
KB.shift_pdf!(rand(16, 2), randn(2), rand(2), 1e-4)

KB.chapman_enskog(rand(16), [1.0, 0.0, 1.0], rand(3), rand(3), 0.1)
KB.chapman_enskog(rand(16), [1.0, 0.0, 1.0], zeros(3), 0, 0.1)
KB.chapman_enskog(
    rand(16, 16),
    rand(16, 16),
    [1.0, 0.0, 0.0, 1.0],
    rand(4),
    rand(4),
    rand(4),
    0.1,
)
KB.chapman_enskog(
    rand(16, 16),
    rand(16, 16),
    [1.0, 0.0, 0.0, 1.0],
    zeros(4),
    zeros(4),
    0.0,
    0.1,
)
KB.chapman_enskog(
    rand(16, 16, 16),
    rand(16, 16, 16),
    rand(16, 16, 16),
    [1.0, 0.0, 0.0, 0.0, 1.0],
    rand(5),
    rand(5),
    rand(5),
    rand(5),
    0.1,
)
KB.chapman_enskog(
    rand(16, 16, 16),
    rand(16, 16, 16),
    rand(16, 16, 16),
    [1.0, 0.0, 0.0, 0.0, 1.0],
    zeros(5),
    zeros(5),
    zeros(5),
    0.0,
    0.1,
)

vs = KB.VSpace1D()
KB.collision_invariant(rand(3), vs)
KB.collision_invariant(rand(3, 2), vs)
vs2 = KB.VSpace2D()
KB.collision_invariant(rand(4), vs2)
vs3 = KB.VSpace3D()
KB.collision_invariant(rand(5), vs3)

#--- quantum ---#
f0 = 0.5 * (1 / π)^0.5 .* (exp.(-(vs.u .- 0.99) .^ 2) .+ exp.(-(vs.u .+ 0.99) .^ 2))
w0 = moments_conserve(f0, vs.u, vs.weights)
prim0 = quantum_conserve_prim(w0, 2, :fd)
quantum_prim_conserve(prim0, 2, :fd)
prim0 = quantum_conserve_prim(w0, 2, :be)
quantum_prim_conserve(prim0, 2, :be)

f0 =
    0.5 * (1 / π)^0.5 .* (exp.(-(vs2.u .- 0.99) .^ 2) .+ exp.(-(vs2.u .+ 0.99) .^ 2)) .*
    exp.(-vs2.v .^ 2)
w0 = moments_conserve(f0, vs2.u, vs2.v, vs2.weights)
prim0 = quantum_conserve_prim(w0, 2, :fd)
quantum_prim_conserve(prim0, 2, :fd)
prim0 = quantum_conserve_prim(w0, 2, :be)
quantum_prim_conserve(prim0, 2, :be)

#--- Hermite ---#
vs = VSpace1D(-5, 5, 36)
prim = [2.0, 0.5, 0.6]
f = maxwellian(vs.u, prim)
df = hermite_force(f, vs.u, vs.weights, prim, 11, 1.0)

#--- Riemann solution ---#
KB.sample_riemann_solution(
    [-0.5, -0.2, 0.1, 0.3, 0.5],
    0.2,
    KB.HydroStatus(1.0, 0.0, 1.0, 1.4),
    KB.HydroStatus(0.125, 0.0, 0.1, 1.4),
)
