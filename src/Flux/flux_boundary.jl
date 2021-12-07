"""
    flux_boundary_maxwell!(fw, bc, w, inK, γ, dt, rot)
    flux_boundary_maxwell!(fw, bc, w, inK, γ, dt, len, rot)
    flux_boundary_maxwell!(fw, fh, fb, bc, h, b, u, ω, inK, dt, rot)
    flux_boundary_maxwell!(fw, ff, bc, f, u, v, ω, dt, len, rot)
    flux_boundary_maxwell!(fw, fh, fb, bc, h, b, u, v, ω, inK, dt, len, rot)
    flux_boundary_maxwell!(fw, ff, bc, f, u, v, w, ω, dt, area, rot)

Maxwell's diffusive boundary flux

- @args: distribution functions and their slopes at left/right sides of interface
- @args: particle velocity quadrature points and weights
- @args: time step

"""
function flux_boundary_maxwell!(
    fw::AbstractVector{T1},
    bc::AbstractVector{T2},
    w::AbstractVector{T3},
    inK,
    γ,
    dt,
    rot, # 1 / -1
) where {T1<:AbstractFloat,T2<:Real,T3<:Real} # 1D continuum

    @assert length(bc) == 3

    primL, primR = ifelse(rot == 1, (bc, conserve_prim(w, γ)), (conserve_prim(w, γ), bc))

    Mu1, Mxi1, MuL1, MuR1 = gauss_moments(primL, inK)
    Mu2, Mxi2, MuL2, MuR2 = gauss_moments(primR, inK)

    Muv_L = moments_conserve(MuL1, Mxi1, 1, 0)
    Muv_R = moments_conserve(MuR2, Mxi2, 1, 0)

    ρ = ifelse(rot == 1, -primR[1] * Muv_R[1] / Muv_L[1], -primL[1] * Muv_L[1] / Muv_R[1])

    @. fw = ρ * (Muv_L + Muv_R) * dt

    return nothing

end

#--- 1D2F1V ---#
function flux_boundary_maxwell!(
    fw::AbstractVector{T1},
    fh::T2,
    fb::T2,
    bc::AbstractVector{T3},
    h::T4,
    b::T4,
    u::T5,
    ω::T5,
    inK,
    dt,
    rot = 1,
) where {
    T1<:AbstractFloat,
    T2<:AbstractVector{<:AbstractFloat},
    T3<:Real,
    T4<:AbstractVector{<:AbstractFloat},
    T5<:AbstractVector{<:AbstractFloat},
} # 1D2F1V

    @assert length(bc) == 3

    δ = heaviside.(u .* rot)
    SF = sum(ω .* u .* h .* (1.0 .- δ))
    SG = (bc[end] / π)^0.5 * sum(ω .* u .* exp.(-bc[end] .* (u .- bc[2]) .^ 2) .* δ)
    prim = [-SF / SG; bc[2:end]]

    H = maxwellian(u, prim)
    B = H .* inK ./ (2.0 * prim[end])

    hWall = H .* δ .+ h .* (1.0 .- δ)
    bWall = B .* δ .+ b .* (1.0 .- δ)

    fw[1] = discrete_moments(hWall, u, ω, 1) * dt
    fw[2] = discrete_moments(hWall, u, ω, 2) * dt
    fw[3] =
        (
            0.5 * discrete_moments(hWall .* u .^ 2, u, ω, 1) +
            0.5 * discrete_moments(bWall, u, ω, 1)
        ) * dt

    @. fh = u * hWall * dt
    @. fb = u * bWall * dt

    return nothing

end

#--- 2D Continuum ---#
function flux_boundary_maxwell!(
    fw::AbstractVector{T1},
    bc::AbstractVector{T2},
    w::AbstractVector{T3},
    inK,
    γ,
    dt,
    len,
    rot,
) where {T1<:Real,T2<:Real,T3<:Real}

    @assert length(bc) == 4

    primL, primR = ifelse(rot == 1, (bc, conserve_prim(w, γ)), (conserve_prim(w, γ), bc))

    Mu1, Mv1, Mxi1, MuL1, MuR1 = gauss_moments(primL, inK)
    Mu2, Mv2, Mxi2, MuL2, MuR2 = gauss_moments(primR, inK)

    Muv_L = moments_conserve(MuL1, Mv1, Mxi1, 1, 0, 0)
    Muv_R = moments_conserve(MuR2, Mv2, Mxi2, 1, 0, 0)

    ρ = ifelse(rot == 1, -primR[1] * Muv_R[1] / Muv_L[1], -primL[1] * Muv_L[1] / Muv_R[1])

    @. fw = ρ * (Muv_L + Muv_R) * dt * len

    return nothing

end

# ------------------------------------------------------------
# 2D1F2V
# ------------------------------------------------------------
function flux_boundary_maxwell!(
    fw::AbstractVector{T1},
    fh::AbstractMatrix{T2},
    bc::AbstractVector{T3},
    h::AbstractMatrix{T4},
    u::T5,
    v::T5,
    ω::T5,
    dt,
    len,
    rot = 1,
) where {
    T1<:AbstractFloat,
    T2<:AbstractFloat,
    T3<:Real,
    T4<:AbstractFloat,
    T5<:AbstractMatrix{<:AbstractFloat},
}

    @assert length(bc) == 4

    δ = heaviside.(u .* rot)

    SF = sum(ω .* u .* h .* (1.0 .- δ))
    SG =
        (bc[end] / π) *
        sum(ω .* u .* exp.(-bc[end] .* ((u .- bc[2]) .^ 2 .+ (v .- bc[3]) .^ 2)) .* δ)
    prim = [-SF / SG; bc[2:end]]

    H = maxwellian(u, v, prim)

    hWall = H .* δ .+ h .* (1.0 .- δ)

    fw[1] = discrete_moments(hWall, u, ω, 1) * len * dt
    fw[2] = discrete_moments(hWall, u, ω, 2) * len * dt
    fw[3] = discrete_moments(hWall .* u, v, ω, 1) * len * dt
    fw[4] = 0.5 * discrete_moments(hWall .* (u .^ 2 .+ v .^ 2), u, ω, 1) * len * dt

    @. fh = u * hWall * len * dt

    return nothing

end

# ------------------------------------------------------------
# 2D2F2V
# ------------------------------------------------------------
function flux_boundary_maxwell!(
    fw::AbstractVector{T1},
    fh::T2,
    fb::T2,
    bc::AbstractVector{T3},
    h::T4,
    b::T4,
    u::T5,
    v::T5,
    ω::T5,
    inK,
    dt,
    len,
    rot = 1,
) where {
    T1<:AbstractFloat,
    T2<:AbstractMatrix{<:AbstractFloat},
    T3<:Real,
    T4<:AbstractMatrix{<:AbstractFloat},
    T5<:AbstractMatrix{<:AbstractFloat},
}

    @assert length(bc) == 4

    δ = heaviside.(u .* rot)

    SF = sum(ω .* u .* h .* (1.0 .- δ))
    SG =
        (bc[end] / π) *
        sum(ω .* u .* exp.(-bc[end] .* ((u .- bc[2]) .^ 2 .+ (v .- bc[3]) .^ 2)) .* δ)
    prim = [-SF / SG; bc[2:end]]

    H = maxwellian(u, v, prim)
    B = H .* inK ./ (2.0 * prim[end])

    hWall = H .* δ .+ h .* (1.0 .- δ)
    bWall = B .* δ .+ b .* (1.0 .- δ)

    fw[1] = discrete_moments(hWall, u, ω, 1) * len * dt
    fw[2] = discrete_moments(hWall, u, ω, 2) * len * dt
    fw[3] = discrete_moments(hWall .* u, v, ω, 1) * len * dt
    fw[4] =
        (
            0.5 * discrete_moments(hWall .* (u .^ 2 .+ v .^ 2), u, ω, 1) +
            0.5 * discrete_moments(bWall, u, ω, 1)
        ) *
        len *
        dt

    @. fh = u * hWall * len * dt
    @. fb = u * bWall * len * dt

    return nothing

end

# ------------------------------------------------------------
# 1F3V
# ------------------------------------------------------------
function flux_boundary_maxwell!(
    fw::AbstractVector{T1},
    ff::AA{T2,3},
    bc::AbstractVector{T3},
    f::AA{T4,3},
    u::T5,
    v::T5,
    w::T5,
    ω::T5,
    dt,
    area,
    rot = 1,
) where {
    T1<:AbstractFloat,
    T2<:AbstractFloat,
    T3<:Real,
    T4<:AbstractFloat,
    T5<:AA{<:AbstractFloat,3},
}

    @assert length(bc) == 5

    δ = heaviside.(u .* rot)
    SF = sum(ω .* u .* f .* (1.0 .- δ))
    SG =
        (bc[end] / π)^1.5 * sum(
            ω .* u .*
            exp.(
                -bc[end] .* ((u .- bc[2]) .^ 2 .+ (v .- bc[3]) .^ 2 .+ (w .- bc[4]) .^ 2),
            ) .* δ,
        )
    prim = [-SF / SG; bc[2:end]]

    M = maxwellian(u, v, w, prim)
    fWall = M .* δ .+ f .* (1.0 .- δ)

    fw[1] = discrete_moments(fWall, u, ω, 1) * area * dt
    fw[2] = discrete_moments(fWall, u, ω, 2) * area * dt
    fw[3] = discrete_moments(fWall .* u, v, ω, 1) * area * dt
    fw[4] = discrete_moments(fWall .* u, w, ω, 1) * area * dt
    fw[5] =
        (0.5 * discrete_moments(fWall .* (u .^ 2 .+ v .^ 2 + w .^ 2), u, ω, 1)) * area * dt

    @. ff = u * fWall * area * dt

    return nothing

end


"""
    flux_boundary_specular!(fw, ff, f, u, ω, dt)
    flux_boundary_specular!(fw, fh, fb, h, b, u, ω, dt)

Specular reflection boundary flux

- @args: fluxes of conservative variables and distribution functions
- @args: distribution functions
- @args: velocity quadrature points and weights
- @args: time step

"""
function flux_boundary_specular!(
    fw::AbstractVector{T1},
    ff::AbstractVector{T2},
    f::AbstractVector{T3},
    u::T4,
    ω::T4,
    dt,
) where {T1<:Real,T2<:AbstractFloat,T3<:AbstractFloat,T4<:AbstractVector{<:AbstractFloat}} # 1D1F1V

    fWall = similar(f)
    for i in eachindex(f)
        fWall[i] = f[end-i+1]
    end

    fw[1] = discrete_moments(fWall, u, ω, 1) * dt
    fw[2] = discrete_moments(fWall, u, ω, 2) * dt
    fw[3] = 0.5 * discrete_moments(fWall .* u .^ 2, u, ω, 1) * dt

    @. ff = u * fWall * dt

    return nothing

end

# ------------------------------------------------------------
# 1D2F1V
# ------------------------------------------------------------
function flux_boundary_specular!(
    fw::AbstractVector{T1},
    fh::T2,
    fb::T2,
    h::T3,
    b::T3,
    u::T4,
    ω::T4,
    dt,
) where {
    T1<:Real,
    T2<:AbstractVector{<:AbstractFloat},
    T3<:AbstractVector{<:AbstractFloat},
    T4<:AbstractVector{<:AbstractFloat},
}

    hWall = similar(h)
    bWall = similar(b)
    for i in eachindex(h)
        hWall[i] = h[end-i+1]
        bWall[i] = b[end-i+1]
    end

    fw[1] = discrete_moments(hWall, u, ω, 1) * dt
    fw[2] = discrete_moments(hWall, u, ω, 2) * dt
    fw[3] =
        (
            0.5 * discrete_moments(hWall .* u .^ 2, u, ω, 1) +
            0.5 * discrete_moments(bWall, u, ω, 1)
        ) * dt

    @. fh = u * hWall * dt
    @. fb = u * bWall * dt

    return nothing

end
