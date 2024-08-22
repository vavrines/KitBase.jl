# ============================================================
# Quadrature Methods
# ============================================================

"""
$(SIGNATURES)

Maxwell quadrature

## Arguments
* `N`: quadrature order (MUST less than 33)
"""
function maxwell_quadrature(N2, C = 1.0)
    N = N2 ÷ 2

    a = zeros(N)
    b = zeros(N)
    a[1] = 1.0 / sqrt(pi)
    a[2] = 2.0 / sqrt(pi) / (pi - 2.0)
    b[2] = a[1] / (a[1] + a[2]) / 2.0

    for i = 3:N
        b[i] = (i - 2) + 1.0 / 2.0 - b[i-1] - a[i-1]^2
        a[i] = ((i - 1)^2 / (4.0 * b[i]) - b[i-1] - 1.0 / 2) / a[i-1] - a[i-1]
    end

    J = Diagonal(a) + diagm(1 => sqrt.(b[2:N])) + diagm(-1 => sqrt.(b[2:N]))

    v, V = eigen(J)

    w = V[1, :] .^ 2 .* sqrt(pi) / 2.0

    vw = [v w]
    vw = sortslices(vw, dims = 1, by = x -> x[1])
    v = vw[:, 1]
    w = vw[:, 2]

    Xis = vcat(-reverse(v), v)
    weights = vcat(reverse(w), w)
    weights .*= exp.(Xis .^ 2) .* C
    Xis .*= C

    return Xis, weights
end


"""
$(SIGNATURES)

Gauss-Legendre quadrature

## Arguments
* `n`: quadrature order (MUST be even)

## Outputs
* `points`: quadrature points in 3D coordinate
* `weights`: quadrature weights
"""
function legendre_quadrature(n::Integer)
    pointsxyz = zeros(n * n, 3)
    weights = zeros(n * n)

    # construct Gauss quadrature
    mu, gaussweights = gausslegendre(n)

    # transform between (mu,phi) and (x,y,z)
    phi = [(k + 0.5) * pi / n for k = 0:2*n-1] # equidistance in z axis
    range = 1:n÷2 # only use upper half of the sphere as quadrature point due to pseudo 3D

    x = sqrt.(1.0 .- mu[range] .^ 2) .* cos.(phi)'
    y = sqrt.(1.0 .- mu[range] .^ 2) .* sin.(phi)'
    z = mu[range] .* ones(size(phi))'
    weights = 2.0 * pi / n * repeat(gaussweights[range], 1, 2 * n)

    # assign
    pointsxyz[:, 1] .= x[:]
    pointsxyz[:, 2] .= y[:]
    pointsxyz[:, 3] .= z[:]
    weights = weights[:]

    return pointsxyz, weights
end


"""
$(SIGNATURES)

Octaeder quadrature

## Arguments
* `n`: quadrature order
* `slerpflag`: flag of spherical linear interpolation

## Outputs
* points and triangulation
"""
function octa_quadrature(n::Integer, slerpflag = true::Bool)
    points, triangulation = octa_triangle(n, slerpflag)
    weights = triangle_weights(points, triangulation)

    return points, weights
end

function octa_triangle(n::Integer, slerpflag = true::Bool)
    # integral range
    pt0 = [0.0, 0.0, 1.0]
    pt1 = [0.0, 1.0, 0.0]
    pt2 = [1.0, 0.0, 0.0]

    # slerp / linspace
    if slerpflag
        pts01 = slerp(pt0, pt1, n)
        pts02 = slerp(pt0, pt2, n)
    else
        pts01 = linspace(pt0, pt1, n)
        pts02 = linspace(pt0, pt2, n)
    end

    nptsoctant = Int64(n * (n + 1) / 2)
    pts = zeros(3, nptsoctant)

    # generate points in planar geometry
    counter = 0
    for i = 1:n
        if slerpflag
            if i == 1
                tmp = pts01[:, 1]
            else
                tmp = slerp(pts01[:, i], pts02[:, i], i)
            end
        else
            tmp = linspace(pts01[i], pts02[i], i)
        end
        for j = 1:i
            counter += 1
            if slerpflag
                pts[:, counter] = tmp[:, j]
            else
                pts[:, counter] = tmp[j]
            end
        end

    end

    # project points onto sphere
    for i = 1:nptsoctant
        pts[:, i] = pts[:, i] / norm(pts[:, i])
    end

    # enumerate over points and write their connectivity
    ids = zeros(Int64, n, n) # matrix that assigns an ID to the points
    nTrianglesOctant = Int64(n * n - 2 * n + 1)
    triangles = zeros(Int64, 3, nTrianglesOctant) # matrix that contains all triangles

    counter = 0
    for i = 1:n
        for j = 1:i
            counter += 1
            ids[i, j] = counter
        end
    end

    # create triangles
    counter = 0
    tmp = zeros(Int64, 1)
    for i = 1:n
        for j = 1:i-1
            tmp = [ids[i, j], ids[i, j+1], ids[i-1, j]]
            counter += 1
            triangles[:, counter] = tmp
        end
        if i < n
            for j = 1:i-1
                tmp = [ids[i, j], ids[i, j+1], ids[i+1, j+1]]
                counter += 1
                triangles[:, counter] = tmp
            end
        end
    end

    # now we have the quadrature points and triangles for a single octant
    ptsAll = deepcopy(pts)

    tmp = deepcopy(pts)
    tmp[1, :] *= -1.0
    ptsAll = deepcopy(hcat(ptsAll, tmp))

    tmp = deepcopy(pts)
    tmp[2, :] *= -1.0
    ptsAll = deepcopy(hcat(ptsAll, tmp))

    tmp = deepcopy(pts)
    tmp[1, :] *= -1.0
    tmp[2, :] *= -1.0
    ptsAll = deepcopy(hcat(ptsAll, tmp))

    tmp = deepcopy(ptsAll)
    tmp[3, :] *= -1.0
    ptsAll = deepcopy(hcat(ptsAll, tmp))

    trianglesAll = deepcopy(triangles)
    for i = 1:7
        trianglesAll = (hcat(trianglesAll, triangles .+ i .* nptsoctant))
    end
    ptsAll, triangulation = unique(ptsAll, trianglesAll)

    points = permutedims(ptsAll)
    triangulation = permutedims(triangulation)

    return points, triangulation
end


"""
$(SIGNATURES)

Create quadrature weights from points and triangulation

## Arguments
* `xyz`: quadrature points
* `triangles`: triangulation

## Outputs
* quadrature weights
"""
function triangle_weights(xyz::AM, triangles::AM)
    weights = zeros(axes(xyz, 1))
    nTriangles = size(triangles, 1)
    xy = zeros(3)
    yz = zeros(3)
    zx = zeros(3)
    mid = zeros(3)

    for n = 1:nTriangles
        # get three nodes of a triangle
        i, j, k = triangles[n, :]

        # get the corresponding points
        x = xyz[i, :]
        y = xyz[j, :]
        z = xyz[k, :]

        # Now get the midpoint of the triangle and the midpoints along the lines
        mid = (x + y + z) / 3.0
        xy = (x + y) / 2.0
        yz = (y + z) / 2.0
        zx = (z + x) / 2.0

        # These points still have to be projected onto the sphere
        mid = mid / norm(mid, 2)
        xy = xy / norm(xy, 2)
        yz = yz / norm(yz, 2)
        zx = zx / norm(zx, 2)

        # By these four points, plus x,y,z we can span 6 triangles
        # Each of these triangles is assigned to one of the three nodes of the triangle x,y,z.
        # the area = the weight

        # i
        weights[i] += area(x, mid, xy, :sphere)
        weights[i] += area(x, mid, zx, :sphere)

        # j
        weights[j] += area(y, mid, xy, :sphere)
        weights[j] += area(y, mid, yz, :sphere)

        # k
        weights[k] += area(z, mid, yz, :sphere)
        weights[k] += area(z, mid, zx, :sphere)
    end

    return weights
end


"""
$(SIGNATURES)

Evaluate quadrature weight from Newton-Cotes rule
"""
function newton_cotes(idx::T, num::T) where {T<:Integer}
    if idx == 1 || idx == num
        nc_coeff = 14.0 / 45.0
    elseif (idx - 5) % 4 == 0
        nc_coeff = 28.0 / 45.0
    elseif (idx - 3) % 4 == 0
        nc_coeff = 24.0 / 45.0
    else
        nc_coeff = 64.0 / 45.0
    end

    return nc_coeff
end
