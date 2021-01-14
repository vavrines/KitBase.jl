using ProgressMeter, LinearAlgebra, Optim#, SphericalHarmonicExpansions
import KitBase

# one-cell test
begin
    quadratureorder = 5
    points, triangulation = KitBase.octa_quadrature(quadratureorder)
    weights = KitBase.quadrature_weights(points, triangulation)
    nq = size(points, 1)
    L = 1
    ne = (L + 1)^2

    α = zeros(ne)
    u = [2., 0., 0., 0.]
    m = KitBase.eval_spherharmonic(points, L)

    res = KitBase.optimize_closure(α, m, weights, u, KitBase.maxwell_boltzmann_dual)
    u1 = KitBase.realizable_reconstruct(res.minimizer, m, weights, KitBase.maxwell_boltzmann_dual_prime)
end

# let's start multi-cell case
begin
    # space
    x0 = -1.5
    x1 = 1.5
    y0 = -1.5
    y1 = 1.5
    nx = 100
    ny = 100
    dx = (x1 - x0) / nx
    dy = (y1 - y0) / ny

    pspace = PSpace2D(x0, x1, nx, y0, y1, ny)

    # time
    tEnd = 1.0
    cfl = 0.95
   
    # quadrature
    quadratureorder = 5
    points, triangulation = octa_quadrature(quadratureorder)
    weights = KitBase.quadrature_weights(points, triangulation)
    nq = size(points, 1)

    # particle
    SigmaS = 1 * ones(ny + 4, nx + 4)
    SigmaA = 0 * ones(ny + 4, nx + 4)
    SigmaT = SigmaS + SigmaA
end

# initial distribution
# maximal moment degree = 1 -> 4 entries
ne = 4 # number of entries
phi = zeros(ne, nx, ny)
s2 = 0.03^2
flr = 1e-4

init_field(x, y) = max(flr, 1.0 / (4.0 * pi * s2) * exp(-(x^2 + y^2) / 4.0 / s2))

for j = 1:nx
    for i = 1:ny
        y = y0 + dy / 2 + (i - 3) * dy
        x = x0 + dx / 2 + (j - 3) * dx
        # only zeroth order moment is non-zero
        phi[1, i, j] = init_field(x, y)
    end
end

### to be done
