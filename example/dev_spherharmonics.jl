# https://github.com/hofmannmartin/SphericalHarmonicExpansions.jl
using SphericalHarmonicExpansions

#--- Yₗᵐ in Cartesian---#
@polyvar x y z

ylm(0, 0, x, y, z) # 1 / 2√π = 0.282
ylm(1, 0, x, y, z) # √(3/π) / 2 = 0.4886
ylm(1, 1, x, y, z)
ylm(2, 2, x, y, z)
p = ylm(3, 3, x, y, z)

# spherical harmonics expansion
L = 1
coeff = SphericalHarmonicCoefficients(L) # length(coeff.c) = (L+1)²
for i in eachindex(coeff.c)
    coeff.c[i] = rand()
end
f = sphericalHarmonicsExpansion(coeff, x, y, z)

# evaluate at quadrature point
f(x=>0.5, y=>-1.0, z=>0.25) # MultivariatePolynomials.jl

g = @fastfunc f # we can create fast evaluating function
g(0.5, -1.0, 0.25)
