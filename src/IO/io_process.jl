extract_sol(ks::AbstractSolverSet, ctr) = extract_sol(ks.ps, ctr)

function extract_sol(ps::AbstractPhysicalSpace1D, ctr)
    sol = zeros(ps.nx, axes(ctr[1].prim)...)
    if ndims(sol) == 2
        for i in axes(sol, 1)
            sol[i, :] .= ctr[i].prim
        end
    elseif ndims(sol) == 3
        for i in axes(sol, 1)
            sol[i, :, :] .= ctr[i].prim
        end
    end

    return sol
end

function extract_sol(ps::AbstractPhysicalSpace2D, ctr)
    nx, ny = begin
        if ps isa CSpace2D
            ps.nr, ps.nθ
        else
            ps.nx, ps.ny
        end
    end

    sol = zeros(nx, ny, axes(ctr[1].prim)...)
    if ndims(sol) == 3
        for i in axes(sol, 1), j in axes(sol, 2)
            sol[i, j, :] .= ctr[i, j].prim
        end
    elseif ndims(sol) == 4
        for i in axes(sol, 1), j in axes(sol, 2)
            sol[i, j, :, :] .= ctr[i, j].prim
        end
    end

    return sol
end

"""
$(SIGNATURES)

Extract wall quantities
"""
function extract_wall(bc, h, b, u, v, ω, inK, γ, rot)
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

    w = moments_conserve(hWall, bWall, u, v, ω, VDF{2,2})
    prim = conserve_prim(w, γ)

    wp = pressure(hWall, bWall, prim, u, v, ω, inK, VDF{2,2})
    ws = sum(ω .* u * v * hWall)
    wq = heat_flux(hWall, bWall, prim, u, v, ω) |> norm

    return wp, ws, wq
end

function extract_wall(ctr, bc, vs, inK, γ, n, rot)
    vn, vt = local_velocity(vs.u, vs.v, n)
    bcL = local_frame(bc, n)

    return extract_wall(bcL, ctr.h, ctr.b, vn, vt, vs.weights, inK, γ, rot)
end
