using KitBase, ProgressMeter, LinearAlgebra, Optim, SphericalHarmonicExpansions

function entropy_prime_dual(y)
    exp(y)
end

function upwind_flux(ψL::Real, ψR::Real, Ω::AbstractVector, n::AbstractVector, dt::Real)
    if dot(Ω, n) > 0
        return dt * ψL
    else
        return dt * ψR
    end
end

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

    nl = 2
    nm = 2
end


@polyvar x y z
L = 1

f = sphericalHarmonicsExpansion(coeff, x, y, z)

m = zeros(ne, nq)
for i in axes(m, 1)
    for j in axes(m, 2)
        coeff = SphericalHarmonicCoefficients(L) # length(coeff.c) = (L+1)²
        for k in eachindex(coeff.c)
            coeff.c[k] = rand()
        end


        m[i, j] = f(x=>0.5, y=>-1.0, z=>0.25)
    end
end



@polyvar x y z

p = ylm(0, 0, x, y, z) 

f = sphericalHarmonicsExpansion(coeff, x, y, z)




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

dt = cfl / 2 * (dx * dy) / (dx + dy)
global t = 0.0

flux1 = zeros(ne, nx + 1, ny)
flux2 = zeros(ne, nx, ny + 1)

@showprogress for iter = 1:10
    for i = 2:nx, j = 1:ny
        flux1[:, i, j] .= calc_flux(phi[:, i-1, j], phi[:, i, j], dt)
    end
    for i = 1:nx, j = 2:ny
        flux2[:, i, j] .= calc_flux(phi[:, i, j-1], phi[:, i, j], dt)
    end

    for j = 1:ny, i = 1:nx
        #integral = discrete_moments(phi[:, i, j], weights)
        #integral *= 1.0 / 4.0 / pi

        for q = 1:ne
        #    phi[q, i, j] =
        #        phi[q, i, j] #+
                #(flux1[q, i, j] - flux1[q, i+1, j]) / dx +
                #(flux2[q, i, j] - flux2[q, i, j+1]) / dy #+
                #(integral - phi[q, i, j]) * dt
        end
    end

    global t += dt
end

using Plots
contourf(pspace.x[1:nx, 1], pspace.y[1, 1:ny], phi[1, :, :])

α = zeros(ne, nx, ny)
m = zeros(ne, nq)

u = zeros(ne, nx, ny)

points ~ (nq, 3)
weights ~ (nq)
