"""
$(SIGNATURES)

Maxwell's diffusive boundary flux

1D continuum

# Arguments

* ``fw``: fluxes for conservative variables
* ``bc``: boundary condition for primitive variables
* ``w``: conservative variables
* ``inK``: internal degrees of freedom
* ``dt``: time step
* ``rot``: rotation indicator (1/-1)
"""
function flux_boundary_maxwell!(
    fw::AV,
    bc::AV,
    w::AV,
    inK,
    γ,
    dt,
    rot, # 1 / -1
)

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

"""
$(SIGNATURES)

1D2F1V
"""
function flux_boundary_maxwell!(
    fw::AV,
    fh::T1,
    fb::T1,
    bc::AV,
    h::T2,
    b::T2,
    u::T3,
    ω::T3,
    inK,
    dt,
    rot = 1,
) where {T1<:AV,T2<:AV,T3<:AV}

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

"""
$(SIGNATURES)

Mixture
"""
function flux_boundary_maxwell!(
    fw::AM,
    fh::T1,
    fb::T1,
    bc::AM,
    h::T2,
    b::T2,
    u::T3,
    ω::T3,
    inK,
    dt,
    rot = 1,
) where {T1<:AM,T2<:AM,T3<:AM}

    @assert size(bc, 1) == 3

    δ = heaviside.(u .* rot)
    SF = [sum(ω[:, j] .* u[:, j] .* h[:, j] .* (1.0 .- δ[:, j])) for j in axes(h, 2)]
    SG = [(bc[end, j] / π)^0.5 * sum(ω[:, j] .* u[:, j] .* exp.(-bc[end, j] .* (u[:, j] .- bc[2, j]) .^ 2) .* δ[:, j]) for j in axes(h, 2)]
    
    prim = zero(bc)
    for j in axes(prim, 2)
        prim[:, j] .= [-SF[j] / SG[j]; bc[2:end, j]]
    end

    H = mixture_maxwellian(u, prim)
    B = mixture_energy_maxwellian(H, prim, inK)

    hWall = H .* δ .+ h .* (1.0 .- δ)
    bWall = B .* δ .+ b .* (1.0 .- δ)

    for j in axes(fw, 2)
        fw[1, j] = discrete_moments(hWall[:, j], u[:, j], ω[:, j], 1) * dt
        fw[2, j] = discrete_moments(hWall[:, j], u[:, j], ω[:, j], 2) * dt
        fw[3, j] =
            (
                0.5 * discrete_moments(hWall[:, j] .* u[:, j] .^ 2, u[:, j], ω[:, j], 1) +
                0.5 * discrete_moments(bWall[:, j], u[:, j], ω[:, j], 1)
            ) * dt
    end

    @. fh = u * hWall * dt
    @. fb = u * bWall * dt

    return nothing

end

"""
$(SIGNATURES)

2D Continuum
"""
function flux_boundary_maxwell!(
    fw::AV,
    bc::AV,
    w::AV,
    inK,
    γ,
    dt,
    len,
    rot,
)

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

"""
$(SIGNATURES)

2D1F2V
"""
function flux_boundary_maxwell!(
    fw::AV,
    fh::AM,
    bc::AV,
    h::AM,
    u::T,
    v::T,
    ω::T,
    dt,
    len,
    rot = 1,
) where {T<:AM}

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

"""
$(SIGNATURES)

2D2F2V
"""
function flux_boundary_maxwell!(
    fw::AV,
    fh::T2,
    fb::T2,
    bc::AV,
    h::T4,
    b::T4,
    u::T5,
    v::T5,
    ω::T5,
    inK,
    dt,
    len,
    rot = 1,
) where {T2<:AM,T4<:AM,T5<:AM}

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

"""
$(SIGNATURES)

1F3V
"""
function flux_boundary_maxwell!(
    fw::AV,
    ff::AA{T1,3},
    bc::AV,
    f::AA{T2,3},
    u::T3,
    v::T3,
    w::T3,
    ω::T3,
    dt,
    area,
    rot = 1,
) where {T1,T2,T3<:AA{T4,3} where {T4}}

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
$(SIGNATURES)

Specular reflection boundary flux

1D1F1V
"""
function flux_boundary_specular!(
    fw::AV,
    ff::AV,
    f::AV,
    u::T,
    ω::T,
    dt,
) where {T<:AV}

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

"""
$(SIGNATURES)

1D2F1V
"""
function flux_boundary_specular!(
    fw::AV,
    fh::T2,
    fb::T2,
    h::T3,
    b::T3,
    u::T4,
    ω::T4,
    dt,
) where {T2<:AV,T3<:AV,T4<:AV}

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
