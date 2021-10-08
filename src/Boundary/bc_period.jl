# ============================================================
# Period Functions
# ============================================================

function bc_period!(ctr::AbstractVector{T}, ng = 1) where {T<:AbstractControlVolume}
    nx = length(ctr) - 2 * ng

    for i = 1:ng
        copy_ctr!(ctr[1-i], ctr[nx+1-i])
        copy_ctr!(ctr[nx+i], ctr[i])
    end

    return nothing
end

function bc_period!(ctr::AbstractMatrix{T}, ng; dirc) where {T<:AbstractControlVolume}
    if dirc == :x
        nx = size(ctr, 1) - 2 * ng
        for j in axes(ctr, 2), i = 1:ng
            copy_ctr!(ctr[1-i, j], ctr[nx-i+1, j])
        end
    elseif dirc == :y
        ny = size(ctr, 2) - 2 * ng
        for j = 1:ng, i in axes(ctr, 1)
            copy_ctr!(ctr[i, 1-j], ctr[i, ny-j+1])
        end
    end

    return nothing
end
