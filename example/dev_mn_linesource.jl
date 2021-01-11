using KitBase, ProgressMeter, LinearAlgebra, Optim

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
    weights = quadrature_weights(points, triangulation)
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

# Lagrange multiplier
α = zeros(ne, nx, ny)
# moment basis
m = zeros(ne, nq)
# work together with solution phi & quadrature weights

# toy example ↓
###
nquad = 16
nentry = 4
m = randn(nentry, nquad)
α = randn(nentry)
ω = ones(nquad) ./ nquad
u = rand(nentry)

function η(f)
    exp(-f)
end

"""
Optimizer for the entropy closure problem
    
argmin(<η(α*m)> - α*u)

"""
function optimize_closure(_α, _m, _ω, _u, _η::Function)
    loss(x) = sum(_η.(x' * _m) .* _ω) - dot(x, _u)
    res = Optim.optimize(loss, _α, Newton()) # Optim.jl
    return res
end

res = optimize_closure(α, m, ω, u, η)
res.minimizer # optimized parameters
###
