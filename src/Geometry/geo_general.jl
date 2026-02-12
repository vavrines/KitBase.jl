"""
$(SIGNATURES)

Transform local flow variables to global frame
"""
function global_frame(w::AV, cosa, sina)
    if eltype(w) <: Int
        G = similar(w, Float64)
    else
        G = similar(w)
    end

    if length(w) == 2
        G[1] = w[1] * cosa - w[2] * sina
        G[2] = w[1] * sina + w[2] * cosa
    elseif length(w) == 4
        G[1] = w[1]
        G[2] = w[2] * cosa - w[3] * sina
        G[3] = w[2] * sina + w[3] * cosa
        G[4] = w[4]
    else
        throw("local -> global: dimension dismatch")
    end

    return G
end

# global_frame(w::AV, n::AV) = global_frame(w, n[1], n[2])

"""
$(SIGNATURES)

Build direction cosine matrix from unit normal vector (3D)

Returns a 3×3 matrix where:
- Row 1: normal vector n
- Row 2: tangent vector t1
- Row 3: tangent vector t2
"""
function direction_cosine(nx, ny, nz)
    if abs(nx) < 0.9
        t1x = 0.0
        t1y = -nz
        t1z = ny
    else
        t1x = -nz
        t1y = 0.0
        t1z = nx
    end
    
    t1_norm = sqrt(t1x^2 + t1y^2 + t1z^2)
    t1x /= t1_norm
    t1y /= t1_norm
    t1z /= t1_norm
    
    t2x = ny * t1z - nz * t1y
    t2y = nz * t1x - nx * t1z
    t2z = nx * t1y - ny * t1x

    dirccos = [
        nx ny nz
        t1x t1y t1z
        t2x t2y t2z
    ]

    return dirccos
end

"""
$(SIGNATURES)
"""
function global_frame(w::AV, dirccos::AM)
    if eltype(w) <: Int
        G = similar(w, Float64)
    else
        G = similar(w)
    end

    if length(w) == 3
        G[1] = w[1] * dirccos[1, 1] + w[2] * dirccos[2, 1] + w[3] * dirccos[3, 1]
        G[2] = w[1] * dirccos[1, 2] + w[2] * dirccos[2, 2] + w[3] * dirccos[3, 2]
        G[3] = w[1] * dirccos[1, 3] + w[2] * dirccos[2, 3] + w[3] * dirccos[3, 3]
    elseif length(w) == 5
        G[1] = w[1]
        G[2] = w[2] * dirccos[1, 1] + w[3] * dirccos[2, 1] + w[4] * dirccos[3, 1]
        G[3] = w[2] * dirccos[1, 2] + w[3] * dirccos[2, 2] + w[4] * dirccos[3, 2]
        G[4] = w[2] * dirccos[1, 3] + w[3] * dirccos[2, 3] + w[4] * dirccos[3, 3]
        G[5] = w[5]
    else
        throw("local -> global: dimension dismatch")
    end

    return G
end

"""
$(SIGNATURES)

Transform from local frame to global frame (3D)
"""
function global_frame(w::AV, nx, ny, nz)
    dirccos = direction_cosine(nx, ny, nz)
    
    return global_frame(w, dirccos)
end

"""
$(SIGNATURES)

Transform from local frame to global frame

For 2D: n = [nx, ny]
For 3D: n = [nx, ny, nz]
"""
function global_frame(w::AV, n)
    if length(n) == 2
        return global_frame(w, n[1], n[2])
    elseif length(n) == 3
        return global_frame(w, n[1], n[2], n[3])
    else
        throw("global_frame: unsupported normal vector dimension")
    end
end

"""
$(SIGNATURES)

Transform global flow variables to local frame
"""
function local_frame(w::AV, cosa, sina)
    if eltype(w) <: Int
        L = similar(w, Float64)
    else
        L = similar(w)
    end

    if length(w) == 2
        L[1] = w[1] * cosa + w[2] * sina
        L[2] = w[2] * cosa - w[1] * sina
    elseif length(w) == 4
        L[1] = w[1]
        L[2] = w[2] * cosa + w[3] * sina
        L[3] = w[3] * cosa - w[2] * sina
        L[4] = w[4]
    else
        throw("global -> local: dimension dismatch")
    end

    return L
end

# """
# $(SIGNATURES)
# """
# local_frame(w::AV, n) = local_frame(w::AV, n[1], n[2])

"""
$(SIGNATURES)
"""
function local_frame(w::AV, dirccos::AM)
    if eltype(w) <: Int
        L = similar(w, Float64)
    else
        L = similar(w)
    end

    if length(w) == 3
        L[1] = w[1] * dirccos[1, 1] + w[2] * dirccos[1, 2] + w[3] * dirccos[1, 3]
        L[2] = w[1] * dirccos[2, 1] + w[2] * dirccos[2, 2] + w[3] * dirccos[2, 3]
        L[3] = w[1] * dirccos[3, 1] + w[2] * dirccos[3, 2] + w[3] * dirccos[3, 3]
    elseif length(w) == 5
        L[1] = w[1]
        L[2] = w[2] * dirccos[1, 1] + w[3] * dirccos[1, 2] + w[4] * dirccos[1, 3]
        L[3] = w[2] * dirccos[2, 1] + w[3] * dirccos[2, 2] + w[4] * dirccos[2, 3]
        L[4] = w[2] * dirccos[3, 1] + w[3] * dirccos[3, 2] + w[4] * dirccos[3, 3]
        L[5] = w[5]
    else
        throw("global -> local: dimension dismatch")
    end

    return L
end

"""
$(SIGNATURES)

Transform from global frame to local frame (3D)
"""
function local_frame(w::AV, nx, ny, nz)
    dirccos = direction_cosine(nx, ny, nz)
    
    return local_frame(w, dirccos)
end

"""
$(SIGNATURES)

Transform from global frame to local frame

For 2D: n = [nx, ny]
For 3D: n = [nx, ny, nz]
"""
function local_frame(w::AV, n)
    if length(n) == 2
        return local_frame(w, n[1], n[2])
    elseif length(n) == 3
        return local_frame(w, n[1], n[2], n[3])
    else
        throw("local_frame: unsupported normal vector dimension")
    end
end

"""
$(SIGNATURES)
"""
function local_velocity(u, v, cosa, sina)
    vn = @. u * cosa + v * sina
    vt = @. v * cosa - u * sina

    return vn, vt
end

local_velocity(u, v, n) = local_velocity(u, v, n[1], n[2])

"""
$(SIGNATURES)

Calculate unit normal vector
"""
function unit_normal(p1::T, p2::T) where {T<:AV}
    Δ = p2 .- p1
    l = norm(Δ) + 1e-6

    return [-Δ[2], Δ[1]] ./ l
end

"""
$(SIGNATURES)
"""
function unit_normal(p1::T, p2::T, p3::T) where {T<:AV}
    v1 = p2 .- p1
    v2 = p3 .- p1

    n = cross(v1, v2)
    l = norm(n) + 1e-6

    return n ./ l
end

"""
$(SIGNATURES)

Calculate point-point/line/surface distance
"""
point_distance(p1::T, p2::T) where {T<:AV} = norm(p1 .- p2)

"""
$(SIGNATURES)
"""
function point_distance(p::T, p1::T, p2::T) where {T<:AV}
    x0, y0 = p
    x1, y1 = p1
    x2, y2 = p2

    return abs((x2 - x1) * (y1 - y0) - (x1 - x0) * (y2 - y1)) /
           sqrt((x2 - x1)^2 + (y2 - y1)^2 + 1e-6)
end
