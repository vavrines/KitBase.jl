# ============================================================
# Extrapolation Functions
# ============================================================

"""
$(SIGNATURES)

Extrapolation boundary condition
"""
function bc_extra!(ctr::AV, ng=1::Integer; dirc)
    if Symbol(dirc) in (:xl, :xL)
        for i in 1:ng
            copy_ctr!(ctr[1-i], ctr[1])
        end
    elseif Symbol(dirc) in (:xr, :xR)
        nx = size(ctr, 1) - 2 * ng
        for i in 1:ng
            copy_ctr!(ctr[nx+i], ctr[nx])
        end
    end

    return nothing
end

"""
$(SIGNATURES)
"""
function bc_extra!(ctr::AM, ng=1::Integer; dirc)
    if Symbol(dirc) in (:xl, :xL)
        for j in axes(ctr, 2), i in 1:ng
            copy_ctr!(ctr[1-i, j], ctr[1, j])
        end
    elseif Symbol(dirc) in (:xr, :xR)
        nx = size(ctr, 1) - 2 * ng
        for j in axes(ctr, 2), i in 1:ng
            copy_ctr!(ctr[nx+i, j], ctr[nx, j])
        end
    elseif Symbol(dirc) in (:yl, :yL)
        for j in 1:ng, i in axes(ctr, 1)
            copy_ctr!(ctr[i, 1-j], ctr[i, 1])
        end
    elseif Symbol(dirc) in (:yr, :yR)
        ny = size(ctr, 2) - 2 * ng
        for j in 1:ng, i in axes(ctr, 1)
            copy_ctr!(ctr[i, ny+j], ctr[i, ny])
        end
    end

    return nothing
end

"""
$(SIGNATURES)
"""
function bc_extra!(ctr::AA{T,3}, ng=1::Integer; dirc) where {T}
    if Symbol(dirc) in (:xl, :xL)
        for k in axes(ctr, 3), j in axes(ctr, 2), i in 1:ng
            copy_ctr!(ctr[1-i, j, k], ctr[1, j, k])
        end
    elseif Symbol(dirc) in (:xr, :xR)
        nx = size(ctr, 1) - 2 * ng
        for k in axes(ctr, 3), j in axes(ctr, 2), i in 1:ng
            copy_ctr!(ctr[nx+i, j, k], ctr[nx, j, k])
        end
    elseif Symbol(dirc) in (:yl, :yL)
        for k in axes(ctr, 3), j in 1:ng, i in axes(ctr, 1)
            copy_ctr!(ctr[i, 1-j, k], ctr[i, 1, k])
        end
    elseif Symbol(dirc) in (:yr, :yR)
        ny = size(ctr, 2) - 2 * ng
        for k in axes(ctr, 3), j in 1:ng, i in axes(ctr, 1)
            copy_ctr!(ctr[i, ny+j, k], ctr[i, ny, k])
        end
    elseif Symbol(dirc) in (:zl, :zL)
        for k in 1:ng, j in axes(ctr, 2), i in axes(ctr, 1)
            copy_ctr!(ctr[i, j, 1-k], ctr[i, j, 1])
        end
    elseif Symbol(dirc) in (:zr, :zR)
        nz = size(ctr, 3) - 2 * ng
        for k in 1:ng, j in axes(ctr, 2), i in axes(ctr, 1)
            copy_ctr!(ctr[i, j, nz+k], ctr[i, j, nz])
        end
    end

    return nothing
end