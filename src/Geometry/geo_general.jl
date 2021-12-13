"""
    2D: global_frame(w::AA{<:Real,1}, cosa, sina)
    3D: global_frame(w::AA{<:Real,1}, dirccos::AA{<:Real,2})

Transform local flow variables to global frame
"""
function global_frame(w::T, cosa, sina) where {T<:AA{<:Real,1}}

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

function global_frame(w::T, dirccos::X) where {T<:AA{<:Real,1},X<:AA{<:Real,2}}
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
    2D: local_frame(w::AA{<:Real,1}, cosa, sina)
    3D: local_frame(w::AA{<:Real,1}, dirccos::AA{<:Real,2})

Transform global flow variables to local frame
"""
function local_frame(w::T, cosa, sina) where {T<:AA{<:Real,1}}

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

function local_frame(w::T, dirccos::X) where {T<:AA{<:Real,1},X<:AA{<:Real,2}}
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
    2D: unit_normal(p1::T, p2::T) where {T<:AV}
    3D: unit_normal(p1::T, p2::T, p3::T) where {T<:AV}

Calculate unit normal vector
"""
function unit_normal(p1::T, p2::T) where {T<:AV}
    Δ = p2 .- p1
    l = norm(Δ) + 1e-6

    return [-Δ[2], Δ[1]] ./ l
end

function unit_normal(p1::T, p2::T, p3::T) where {T<:AV}
    v1 = p2 .- p1
    v2 = p3 .- p1

    n = cross(v1, v2)
    l = norm(n) + 1e-6

    return n ./ l
end


"""
    point_distance(p1::T, p2::T) where {T<:AV}
    point_distance(p::T, p1::T, p2::T) where {T<:AV}

Calculate point-point/line/surface distance
"""
point_distance(p1::T, p2::T) where {T<:AV} = norm(p1 .- p2)

function point_distance(p::T, p1::T, p2::T) where {T<:AV}
    x0, y0 = p
    x1, y1 = p1
    x2, y2 = p2

    return abs((x2 - x1) * (y1 - y0) - (x1 - x0) * (y2 - y1)) /
           sqrt((x2 - x1)^2 + (y2 - y1)^2 + 1e-6)
end
