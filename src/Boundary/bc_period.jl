# ============================================================
# Period Functions
# ============================================================

"""
$(SIGNATURES)

Periodic boundary condition
"""
function bc_period!(ctr::AV, ng=1)
    nx = length(ctr) - 2 * ng

    for i in 1:ng
        copy_ctr!(ctr[1-i], ctr[nx+1-i])
        copy_ctr!(ctr[nx+i], ctr[i])
    end

    return nothing
end

"""
$(SIGNATURES)
"""
function bc_period!(ctr::AM, ng; dirc)
    if dirc == :x
        nx = size(ctr, 1) - 2 * ng
        for j in axes(ctr, 2), i in 1:ng
            copy_ctr!(ctr[1-i, j], ctr[nx-i+1, j])
        end
    elseif dirc == :y
        ny = size(ctr, 2) - 2 * ng
        for j in 1:ng, i in axes(ctr, 1)
            copy_ctr!(ctr[i, 1-j], ctr[i, ny-j+1])
        end
    end

    return nothing
end

"""
$(SIGNATURES)

Periodic boundary condition for 3D
"""
function bc_period!(ctr::AA{T,3}, ng; dirc) where {T}
    if Symbol(dirc) in (:x, :X)
        nx = size(ctr, 1) - 2 * ng
        for k in axes(ctr, 3), j in axes(ctr, 2), i in 1:ng
            copy_ctr!(ctr[1-i, j, k], ctr[nx-i+1, j, k])
            copy_ctr!(ctr[nx+i, j, k], ctr[i, j, k])
        end
    elseif Symbol(dirc) in (:y, :Y)
        ny = size(ctr, 2) - 2 * ng
        for k in axes(ctr, 3), j in 1:ng, i in axes(ctr, 1)
            copy_ctr!(ctr[i, 1-j, k], ctr[i, ny-j+1, k])
            copy_ctr!(ctr[i, ny+j, k], ctr[i, j, k])
        end
    elseif Symbol(dirc) in (:z, :Z)
        nz = size(ctr, 3) - 2 * ng
        for k in 1:ng, j in axes(ctr, 2), i in axes(ctr, 1)
            copy_ctr!(ctr[i, j, 1-k], ctr[i, j, nz-k+1])
            copy_ctr!(ctr[i, j, nz+k], ctr[i, j, k])
        end
    end

    return nothing
end
