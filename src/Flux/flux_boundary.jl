"""
$(SIGNATURES)

Maxwell's diffusive boundary flux

1D continuum

## Arguments
* `fw`: fluxes for conservative variables
* `bc`: boundary condition for primitive variables
* `w`: conservative variables
* `inK`: internal degrees of freedom
* `dt`: time step
* `rot`: rotation indicator (1/-1)
"""
function flux_boundary_maxwell!(fw, bc::AV, w, inK, γ, dt, rot)

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

1D2F1V & 1F3V
"""
function flux_boundary_maxwell!(fw::AV, a1, a2, a3, a4, a5, a6, a7, a8, a9, rot)

    if length(fw) == 3
        fh, fb, bc, h, b, u, ω, inK, dt = a1, a2, a3, a4, a5, a6, a7, a8, a9

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
    elseif length(fw) == 5
        ff, bc, f, u, v, w, ω, dt, area = a1, a2, a3, a4, a5, a6, a7, a8, a9

        δ = heaviside.(u .* rot)
        SF = sum(ω .* u .* f .* (1.0 .- δ))
        SG =
            (bc[end] / π)^1.5 * sum(
                ω .* u .*
                exp.(
                    -bc[end] .*
                    ((u .- bc[2]) .^ 2 .+ (v .- bc[3]) .^ 2 .+ (w .- bc[4]) .^ 2),
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
            (0.5 * discrete_moments(fWall .* (u .^ 2 .+ v .^ 2 + w .^ 2), u, ω, 1)) *
            area *
            dt

        @. ff = u * fWall * area * dt
    end


    return nothing

end

"""
$(SIGNATURES)

1D2F1V mixture
"""
function flux_boundary_maxwell!(fw::AM, fh, fb, bc, h, b, u, ω, inK, dt, rot)

    @assert size(bc, 1) == 3

    δ = heaviside.(u .* rot)
    SF = [sum(ω[:, j] .* u[:, j] .* h[:, j] .* (1.0 .- δ[:, j])) for j in axes(h, 2)]
    SG = [
        (bc[end, j] / π)^0.5 * sum(
            ω[:, j] .* u[:, j] .* exp.(-bc[end, j] .* (u[:, j] .- bc[2, j]) .^ 2) .*
            δ[:, j],
        ) for j in axes(h, 2)
    ]

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
function flux_boundary_maxwell!(fw, bc::AV, w, inK, γ, dt, len, rot)

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
function flux_boundary_maxwell!(fw, fh, bc::AV, h, u, v, ω, dt, len, rot)

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
function flux_boundary_maxwell!(fw, fh, fb, bc::AV, h, b, u, v, ω, inK, dt, len, rot)

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

Specular reflection boundary flux

1D1F1V
"""
function flux_boundary_specular!(fw, ff, f, u, ω, dt)
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
function flux_boundary_specular!(fw, fh, fb, h, b, u, ω, dt)
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
