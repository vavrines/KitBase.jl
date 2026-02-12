# ============================================================
# Mirroring Functions
# ============================================================

"""
$(SIGNATURES)

Mirror boundary condition
"""
function bc_mirror!(ctr::AM, ng=1::Integer; dirc)
    if Symbol(dirc) in (:xl, :xL)
        for j in axes(ctr, 2), i in 1:ng
            bc_mirror!(ctr[1-i, j], ctr[i, j], :x)
        end
    elseif Symbol(dirc) in (:xr, :xR)
        nx = size(ctr, 1) - 2 * ng
        for j in axes(ctr, 2), i in 1:ng
            bc_mirror!(ctr[nx+i, j], ctr[nx-i+1, j], :x)
        end
    elseif Symbol(dirc) in (:yl, :yL)
        for j in 1:ng, i in axes(ctr, 1)
            bc_mirror!(ctr[i, 1-j], ctr[i, j], :y)
        end
    elseif Symbol(dirc) in (:yr, :yR)
        ny = size(ctr, 2) - 2 * ng
        for j in 1:ng, i in axes(ctr, 1)
            bc_mirror!(ctr[i, ny+j], ctr[i, ny-j+1], :y)
        end
    end

    return nothing
end

"""
$(SIGNATURES)

Mirror boundary condition for 3D
"""
function bc_mirror!(ctr::AA{T,3}, ng=1::Integer; dirc) where {T}
    if Symbol(dirc) in (:xl, :xL)
        for k in axes(ctr, 3), j in axes(ctr, 2), i in 1:ng
            bc_mirror!(ctr[1-i, j, k], ctr[i, j, k], :x)
        end
    elseif Symbol(dirc) in (:xr, :xR)
        nx = size(ctr, 1) - 2 * ng
        for k in axes(ctr, 3), j in axes(ctr, 2), i in 1:ng
            bc_mirror!(ctr[nx+i, j, k], ctr[nx-i+1, j, k], :x)
        end
    elseif Symbol(dirc) in (:yl, :yL)
        for k in axes(ctr, 3), j in 1:ng, i in axes(ctr, 1)
            bc_mirror!(ctr[i, 1-j, k], ctr[i, j, k], :y)
        end
    elseif Symbol(dirc) in (:yr, :yR)
        ny = size(ctr, 2) - 2 * ng
        for k in axes(ctr, 3), j in 1:ng, i in axes(ctr, 1)
            bc_mirror!(ctr[i, ny+j, k], ctr[i, ny-j+1, k], :y)
        end
    elseif Symbol(dirc) in (:zl, :zL)
        for k in 1:ng, j in axes(ctr, 2), i in axes(ctr, 1)
            bc_mirror!(ctr[i, j, 1-k], ctr[i, j, k], :z)
        end
    elseif Symbol(dirc) in (:zr, :zR)
        nz = size(ctr, 3) - 2 * ng
        for k in 1:ng, j in axes(ctr, 2), i in axes(ctr, 1)
            bc_mirror!(ctr[i, j, nz+k], ctr[i, j, nz-k+1], :z)
        end
    end

    return nothing
end

"""
$(SIGNATURES)
"""
function bc_mirror!(ctr::T, ctr0::T, dirc) where {T<:Union{ControlVolume,ControlVolume2D,ControlVolume3D}}
    copy_ctr!(ctr, ctr0)

    if Symbol(dirc) == :x
        ctr.w[2] *= -1.0
        ctr.prim[2] *= -1.0
        ctr.sw[2, 1] *= -1.0
    elseif Symbol(dirc) == :y
        ctr.w[3] *= -1.0
        ctr.prim[3] *= -1.0
        ctr.sw[3, 2] *= -1.0
    elseif Symbol(dirc) == :z
        ctr.w[4] *= -1.0
        ctr.prim[4] *= -1.0
        ctr.sw[4, 3] *= -1.0
    end

    return nothing
end

"""
$(SIGNATURES)
"""
function bc_mirror!(
    ctr::T,
    ctr0::T,
    dirc,
) where {T<:Union{ControlVolume1F,ControlVolume2D1F}}
    copy_ctr!(ctr, ctr0)

    nu = size(ctr.f, 1)
    nv = size(ctr.f, 2)

    if Symbol(dirc) == :x
        ctr.w[2] *= -1.0
        ctr.prim[2] *= -1.0
        ctr.sw[2, 1] *= -1.0

        for j in 1:nv, i in 1:nu
            ctr.f[i, j] = ctr0.f[nu+1-i, j]
        end
    elseif Symbol(dirc) == :y
        ctr.w[3] *= -1.0
        ctr.prim[3] *= -1.0
        ctr.sw[3, 2] *= -1.0

        for j in 1:nv, i in 1:nu
            ctr.f[i, j] = ctr0.f[i, nv+1-j]
        end
    end

    return nothing
end

"""
$(SIGNATURES)
"""
function bc_mirror!(
    ctr::T,
    ctr0::T,
    dirc,
) where {T<:Union{ControlVolume2F,ControlVolume2D2F}}
    copy_ctr!(ctr, ctr0)

    nu = size(ctr.h, 1)
    nv = size(ctr.h, 2)

    if Symbol(dirc) == :x
        ctr.w[2] *= -1.0
        ctr.prim[2] *= -1.0
        ctr.sw[2, 1] *= -1.0

        for j in 1:nv, i in 1:nu
            ctr.h[i, j] = ctr0.h[nu+1-i, j]
            ctr.b[i, j] = ctr0.b[nu+1-i, j]
        end
    elseif Symbol(dirc) == :y
        ctr.w[3] *= -1.0
        ctr.prim[3] *= -1.0
        ctr.sw[3, 2] *= -1.0

        for j in 1:nv, i in 1:nu
            ctr.h[i, j] = ctr0.h[i, nv+1-j]
            ctr.b[i, j] = ctr0.b[i, nv+1-j]
        end
    end

    return nothing
end
