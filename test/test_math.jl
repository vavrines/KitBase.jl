x0 = 0
x1 = 1
nx = 5

KB.linspace(x0, x1, nx)
KB.heaviside(x0)
KB.fortsign(x0, x1)

KB.mat_split(randn(2, 2))
KB.mat_split(randn(3, 3))

KB.convergence_order(1e-2, 1e-3)
KB.L1_error(rand(3), rand(3), 1e-2)
KB.L2_error(rand(3), rand(3), 1e-2)
KB.L∞_error(rand(3), rand(3), 1e-2)

x = randn(16)
y = randn(16)
KB.@nametuple x y # NamedTuple constructor
KB.central_diff(y, x)
KB.central_diff(y, x0)
KB.central_diff2(y, x)
KB.central_diff2(y, x0)

res = similar(x)
KB.central_diff!(res, y, x)
KB.central_diff!(res, y, x0)
KB.central_diff2!(res, y, x)
KB.central_diff2!(res, y, x0)

KB.upwind_diff(y, x; stream=:right)
KB.upwind_diff(y, x; stream=:left)
KB.upwind_diff(y, x0)

KB.upwind_diff!(res, y, x)
KB.upwind_diff!(res, y, x0)

KB.unstruct_diff(y, x, 4; mode=:central)
KB.unstruct_diff(y, x, 4; mode=:upwind)
KB.unstruct_diff(sin, randn(12), 4, 1)
KB.unstruct_diff(sin, randn(12), 4, 2)

KB.lgwt(12, -1, 1)

#--- entropy closure ---#
KB.maxwell_boltzmann(rand())
KB.maxwell_boltzmann_prime(rand())
KB.maxwell_boltzmann_dual(rand())
KB.maxwell_boltzmann_dual_prime(rand())

quadratureorder = 2
points, weights = KB.octa_quadrature(quadratureorder)
nq = size(points, 1)
L = 1
ne = (L + 1)^2

α = zeros(ne)
u = [2.0, 0.0, 0.0, 0.0]
m = KB.eval_spherharmonic(points, L)

KB.eval_sphermonomial(rand(6), L)
KB.eval_sphermonomial(points, L)

KB.kinetic_entropy(α, m, weights)

using KitBase.TypedPolynomials
@polyvar _x _y _z
KB.rlylm(2, 2, _x, _y, _z)

res = KB.optimize_closure(α, m, weights, u, KB.maxwell_boltzmann_dual)
u1 = KB.realizable_reconstruct(res.minimizer, m, weights, KB.maxwell_boltzmann_dual_prime)

using KitBase.Distributions
u = KB.linspace(-5, 5, 100)
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
KB.polylog(s, z)

# abs(z) <= 0.5 && abs(z) < t
z = 0.1
μ, t = get_μt(z)
KB.polylog(s, z)

# t <= T && (abs(round(real(s)) - s) > tau_threshold || real(s) <= 0)
s, z = -1.0, 1.0
μ, t = get_μt(z)
KB.polylog(s, z)

# t <= T = 0.512
s, z = 1.5, 0.5
μ, t = get_μt(z)
KB.polylog(s, z)

# else
s, z = 1.0, 100.0
μ, t = get_μt(z)
KB.polylog(s, z)

KB.harmonic(9)
KB.harmonic(9.0)
KB.harmonic(0, 9.0)
KB.harmonic(1, 9.0)
KB.harmonic(2, 9.0)
KB.harmonic(2, 1)

KB.f_crandall(0, 0)
KB.f_crandall(0, 1)
KB.f_crandall(1, 0)
KB.f_crandall(1, 1)

KB.g_crandall(1)
KB.b_crandall(1, 1, 1)
KB.c_crandall(1, 1, 1)

KB.Q(1, 0, 1)
KB.Q(1, 1, 1)
KB.Q_closed(1, 1, 1)

KB.c_closed(1, 0, 1)
KB.c_closed(1, 1, 1)
KB.c_closed(1, 2, 1)

KB.dirac_delta(1)
KB.dirac_delta(1, KB.Class{2})
KB.dirac_delta(1, KB.Class{3})
KB.dirac_delta(1, 1, 0.1, 0.1)
KB.dirac_delta(1, 1, 0.1, 0.1, KB.Class{3})
