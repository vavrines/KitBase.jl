x0 = 0
x1 = 1
nx = 5

KitBase.linspace(x0, x1, nx)
KitBase.heaviside(x0)
KitBase.fortsign(x0, x1)

KitBase.mat_split(randn(2, 2))
KitBase.mat_split(randn(3, 3))

KB.convergence_order(1e-2, 1e-3)
KB.L1_error(rand(3), rand(3), 1e-2)
KB.L2_error(rand(3), rand(3), 1e-2)
KB.L∞_error(rand(3), rand(3), 1e-2)

u = Float64[0, 1, 2, 3, 2, 1, 0]
dx = 1.0
KitBase.finite_difference(u, dx; method = :central, bc = :period)
KitBase.finite_difference(u, dx; method = :central, bc = :none)
KitBase.finite_difference(u, dx, -1; method = :upwind, bc = :none)

x = randn(16)
y = randn(16)
KitBase.@nametuple x y # NamedTuple constructor
KitBase.central_diff(y, x)
KitBase.central_diff(y, x0)
KitBase.central_diff2(y, x)
KitBase.central_diff2(y, x0)

res = similar(x)
KitBase.central_diff!(res, y, x)
KitBase.central_diff!(res, y, x0)
KitBase.central_diff2!(res, y, x)
KitBase.central_diff2!(res, y, x0)

KitBase.upwind_diff(y, x; stream = :right)
KitBase.upwind_diff(y, x; stream = :left)
KitBase.upwind_diff(y, x0)

KitBase.upwind_diff!(res, y, x)
KitBase.upwind_diff!(res, y, x0)

KitBase.unstruct_diff(y, x, 4; mode = :central)
KitBase.unstruct_diff(y, x, 4; mode = :upwind)
KitBase.unstruct_diff(sin, randn(12), 4, 1)
KitBase.unstruct_diff(sin, randn(12), 4, 2)

KitBase.lgwt(12, -1, 1)

KitBase.extract_last(randn(2, 3), 2; mode = :view)
KitBase.extract_last(randn(2, 3), 2; mode = :copy)

#--- entropy closure ---#
KitBase.maxwell_boltzmann(rand())
KitBase.maxwell_boltzmann_prime(rand())
KitBase.maxwell_boltzmann_dual(rand())
KitBase.maxwell_boltzmann_dual_prime(rand())

quadratureorder = 2
points, weights = KitBase.octa_quadrature(quadratureorder)
nq = size(points, 1)
L = 1
ne = (L + 1)^2

α = zeros(ne)
u = [2.0, 0.0, 0.0, 0.0]
m = KitBase.eval_spherharmonic(points, L)

KitBase.eval_sphermonomial(rand(6), L)
KitBase.eval_sphermonomial(points, L)

KitBase.kinetic_entropy(α, m, weights)

using KitBase.TypedPolynomials
@polyvar _x _y _z
KitBase.rlylm(2, 2, _x, _y, _z)

res = KitBase.optimize_closure(α, m, weights, u, KitBase.maxwell_boltzmann_dual)
u1 = KitBase.realizable_reconstruct(
    res.minimizer,
    m,
    weights,
    KitBase.maxwell_boltzmann_dual_prime,
)

using Distributions
u = KitBase.linspace(-5, 5, 100)
m = moment_basis(u, 4)
mm = moment_basis(u, u, 4)
mmm = moment_basis(u, u, u, 4)

pdf = Normal(0, 0.01)
prim = [1.0, 0.0, 0.8]
sample_pdf(m, 4, prim, pdf)
sample_pdf(mm, 4, [1, 0, 0, 1], pdf)
sample_pdf(mmm, 4, [1, 0, 0, 0, 1], pdf)

###
# Polylogarithms
###

function get_μt(z)
    μ = log(convert(Complex{Float64}, z))
    t = abs(μ / 2 / π)
    return μ, t
end

# abs(μ) < 1.0e-14
s, z = 1.5, 1.0
μ, t = get_μt(z)
polylog(s, z)

# abs(z) <= 0.5 && abs(z) < t
z = 0.1
μ, t = get_μt(z)
polylog(s, z)

# t <= T && (abs(round(real(s)) - s) > tau_threshold || real(s) <= 0)
s, z = -1.0, 1.0
μ, t = get_μt(z)
polylog(s, z)

# t <= T = 0.512
s, z = 1.5, 0.5
μ, t = get_μt(z)
polylog(s, z)

# else
s, z = 1.0, 100.0
μ, t = get_μt(z)
polylog(s, z)

KitBase.harmonic(9)
KitBase.harmonic(9.0)
KitBase.harmonic(0, 9.0)
KitBase.harmonic(1, 9.0)
KitBase.harmonic(2, 9.0)
KitBase.harmonic(2, 1)

KitBase.f_crandall(0, 0)
KitBase.f_crandall(0, 1)
KitBase.f_crandall(1, 0)
KitBase.f_crandall(1, 1)

KitBase.g_crandall(1)
KitBase.b_crandall(1, 1, 1)
KitBase.c_crandall(1, 1, 1)

KitBase.Q(1, 0, 1)
KitBase.Q(1, 1, 1)
KitBase.Q_closed(1, 1, 1)

KitBase.c_closed(1, 0, 1)
KitBase.c_closed(1, 1, 1)
KitBase.c_closed(1, 2, 1)

KitBase.dirac_delta(1)
KitBase.dirac_delta(1, KitBase.Class{2})
KitBase.dirac_delta(1, KitBase.Class{3})
KitBase.dirac_delta(1, 1, 0.1, 0.1)
KitBase.dirac_delta(1, 1, 0.1, 0.1, KitBase.Class{3})
