using KitBase, Plots
using ProgressMeter: @showprogress

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
    quadratureorder = 8
    #points, weights = legendre_quadrature(quadratureorder)
    points, triangulation = octa_quadrature(quadratureorder)
    weights = quadrature_weights(points, triangulation)
    nq = size(points, 1)

    # particle
    SigmaS = 1 .* ones(ny, nx)
    SigmaA = 0 .* ones(ny, nx)
    SigmaT = SigmaS + SigmaA
end

# initial distribution
phi = ones(nq, nx, ny) .* 1e-4
s2 = 0.03^2
flr = 1e-4
init_field(x, y) = max(flr, 1.0 / (4.0 * pi * s2) * exp(-(x^2 + y^2) / 4.0 / s2))
for j = 1:nx
    for i = 1:ny
        y = y0 + (i - 0.5) * dy
        x = x0 + (j - 0.5) * dx
        for q = 1:nq
            phi[q, i, j] = init_field(x, y)

            #if -0.01 <= pspace.x[i, j] <= 0.01 && -0.01 <= pspace.y[i, j] <= 0.01
            #    phi[q, i, j] = 1.0
            #end
        end
    end
end

dt = cfl / 2 * (dx * dy) / (dx + dy)
global t = 0.0

flux1 = zeros(nq, nx + 1, ny)
flux2 = zeros(nq, nx, ny + 1)

@showprogress for iter = 1:20
    for i = 2:nx, j = 1:ny
        tmp = @view flux1[:, i, j]
        flux_kfvs!(tmp, phi[:, i-1, j], phi[:, i, j], points[:, 1], dt)
    end
    for i = 1:nx, j = 2:ny
        tmp = @view flux2[:, i, j]
        flux_kfvs!(tmp, phi[:, i, j-1], phi[:, i, j], points[:, 2], dt)
    end

    for j = 1:ny, i = 1:nx
        integral = discrete_moments(phi[:, i, j], weights)
        integral /= 4.0 * π

        for q = 1:nq
            phi[q, i, j] =
                phi[q, i, j] +
                (flux1[q, i, j] - flux1[q, i+1, j]) / dx +
                (flux2[q, i, j] - flux2[q, i, j+1]) / dy +
                (integral - phi[q, i, j]) * dt
        end
    end

    global t += dt
end

ρ = zeros(nx, ny)
for i = 1:nx, j = 1:ny
    ρ[i, j] = discrete_moments(phi[:, i, j], weights)
end

contourf(pspace.x[1:nx, 1], pspace.y[1, 1:ny], ρ', color = :PiYG_3)
