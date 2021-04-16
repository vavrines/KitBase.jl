# Ref: https://github.com/hofmannmartin/SphericalHarmonicExpansions.jl
using MultivariatePolynomials, TypedPolynomials
import KitBase

#--- Yₗᵐ in Cartesian---#
@polyvar x y z

KitBase.ylm(0, 0, x, y, z) # 1 / 2√π = 0.282
KitBase.ylm(1, 0, x, y, z) # √(3/π) / 2 = 0.4886
KitBase.ylm(1, 1, x, y, z)
KitBase.ylm(2, 2, x, y, z)
f = KitBase.ylm(3, 3, x, y, z)

# evaluate at quadrature point
f(x => 0.5, y => -1.0, z => 0.25) # MultivariatePolynomials.jl
