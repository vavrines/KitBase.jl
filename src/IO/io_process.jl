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
    sol = zeros(ps.nx, ps.ny, axes(ctr[1].prim)...)
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
