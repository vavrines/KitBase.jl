"""
$(SIGNATURES)

Update solver for boundary cells
"""
function update_boundary!(
    KS,
    ctr::AV{TC},
    face,
    dt,
    residual;
    bc,
    fn = step!,
    kwargs...,
) where {
    TC<:Union{
        ControlVolume,
        ControlVolume1D,
        ControlVolume1F,
        ControlVolume1D1F,
        ControlVolume2F,
        ControlVolume1D2F,
    },
}

    resL = zero(ctr[1].w)
    avgL = zero(ctr[1].w)
    resR = zero(ctr[1].w)
    avgR = zero(ctr[1].w)

    if bc[1] != :fix
        i = 1
        fn(KS, ctr[i], face[i], face[i+1], (dt, KS.ps.dx[i], resL, avgL))
    end

    if bc[2] != :fix
        j = KS.ps.nx
        fn(KS, ctr[j], face[j], face[j+1], (dt, KS.ps.dx[j], resR, avgR))
    end

    #@. residual += sqrt((resL + resR) * 2) / (avgL + avgR + 1.e-7)
    for i in eachindex(residual)
        residual[i] += sqrt((resL[i] + resR[i]) * 2) / (avgL[i] + avgR[i] + 1.e-7) # residual for first species only
    end

    ng = 1 - first(eachindex(KS.ps.x))
    if bc[1] == :period
        bc_period!(ctr, ng)
    elseif bc[1] == :extra
        bc_extra!(ctr, ng; dirc = :xl)
    elseif bc[1] == :balance
        bc_balance!(ctr[0], ctr[1], ctr[2])
    end
    if bc[2] == :extra
        bc_extra!(ctr, ng; dirc = :xr)
    elseif bc[2] == :balance
        bc_balance!(ctr[KS.ps.nx+1], ctr[KS.ps.nx], ctr[KS.ps.nx-1])
    end

    return nothing

end

function update_boundary!(
    KS,
    ctr::AV{TC},
    face,
    dt,
    residual;
    coll = symbolize(KS.set.collision),
    bc,
    fn = step!,
    isMHD = false,
) where {TC<:Union{ControlVolume3F,ControlVolume1D3F}}

    resL = zero(ctr[1].w)
    avgL = zero(ctr[1].w)
    resR = zero(ctr[1].w)
    avgR = zero(ctr[1].w)

    if bc[1] != :fix
        i = 1
        fn(KS, ctr[i], face[i], face[i+1], KS.ps.dx[i], dt, resL, avgL, coll, isMHD)
    end

    if bc[2] != :fix
        j = KS.ps.nx
        fn(KS, ctr[j], face[j], face[j+1], KS.ps.dx[j], dt, resR, avgR, coll, isMHD)
    end

    for i in eachindex(residual)
        residual[i] += sqrt((resL[i] + resR[i]) * 2) / (avgL[i] + avgR[i] + 1.e-7)
    end

    ng = 1 - first(eachindex(KS.ps.x))
    if bc[1] == :period
        bc_period!(ctr, ng)
    elseif bc[1] == :extra
        bc_extra!(ctr, ng; dirc = :xl)
    elseif bc[1] == :balance
        bc_balance!(ctr[0], ctr[1], ctr[2])
    end
    if bc[2] == :extra
        bc_extra!(ctr, ng; dirc = :xr)
    elseif bc[2] == :balance
        bc_balance!(ctr[KS.ps.nx+1], ctr[KS.ps.nx], ctr[KS.ps.nx-1])
    end

    return nothing

end

function update_boundary!(
    KS,
    ctr::AV{TC},
    face,
    dt,
    residual;
    coll = symbolize(KS.set.collision),
    bc,
    fn = step!,
    isMHD = false::Bool,
) where {TC<:Union{ControlVolume4F,ControlVolume1D4F}}

    resL = zero(ctr[1].w)
    avgL = zero(ctr[1].w)
    resR = zero(ctr[1].w)
    avgR = zero(ctr[1].w)

    if bc[1] != :fix
        i = 1
        fn(KS, ctr[i], face[i], face[i+1], KS.ps.dx[i], dt, resL, avgL, coll, isMHD)
    end

    if bc[2] != :fix
        j = KS.ps.nx
        fn(KS, ctr[j], face[j], face[j+1], KS.ps.dx[j], dt, resR, avgR, coll, isMHD)
    end

    for i in eachindex(residual)
        residual[i] += sqrt((resL[i] + resR[i]) * 2) / (avgL[i] + avgR[i] + 1.e-7)
    end

    ng = 1 - first(eachindex(KS.ps.x))
    if bc[1] == :period
        bc_period!(ctr, ng)
    elseif bc[1] == :extra
        bc_extra!(ctr, ng; dirc = :xl)
    elseif bc[1] == :balance
        bc_balance!(ctr[0], ctr[1], ctr[2])
    end
    if bc[2] == :extra
        bc_extra!(ctr, ng; dirc = :xr)
    elseif bc[2] == :balance
        bc_balance!(ctr[KS.ps.nx+1], ctr[KS.ps.nx], ctr[KS.ps.nx-1])
    end

    return nothing

end

function update_boundary!(
    KS,
    ctr::AM{TC},
    a1face,
    a2face,
    dt,
    residual;
    coll = symbolize(KS.set.collision),
    bc,
    fn = step!,
) where {
    TC<:Union{
        ControlVolume,
        ControlVolume2D,
        ControlVolume1F,
        ControlVolume2D1F,
        ControlVolume2F,
        ControlVolume2D2F,
    },
}

    nx, ny, dx, dy = begin
        if KS.ps isa CSpace2D
            KS.ps.nr, KS.ps.nθ, KS.ps.dr, KS.ps.darc
        else
            KS.ps.nx, KS.ps.ny, KS.ps.dx, KS.ps.dy
        end
    end

    resL = zero(ctr[1].w)
    avgL = zero(ctr[1].w)
    resR = zero(ctr[1].w)
    avgR = zero(ctr[1].w)
    resU = zero(ctr[1].w)
    avgU = zero(ctr[1].w)
    resD = zero(ctr[1].w)
    avgD = zero(ctr[1].w)

    if bc[1] != :fix
        @inbounds for j = 1:ny
            fn(
                KS,
                ctr[1, j],
                a1face[1, j],
                a1face[2, j],
                a2face[1, j],
                a2face[1, j+1],
                (dt, dx[1, j] * dy[1, j], resL, avgL),
                coll,
            )
        end
    end

    if bc[2] != :fix
        @inbounds for j = 1:ny
            fn(
                KS,
                ctr[nx, j],
                a1face[nx, j],
                a1face[nx+1, j],
                a2face[nx, j],
                a2face[nx, j+1],
                (dt, dx[nx, j] * dy[nx, j], resR, avgR),
                coll,
            )
        end
    end

    if bc[3] != :fix
        @inbounds for i = 2:nx-1 # skip overlap
            fn(
                KS,
                ctr[i, 1],
                a1face[i, 1],
                a1face[i+1, 1],
                a2face[i, 1],
                a2face[i, 2],
                (dt, dx[i, 1] * dy[i, 1], resD, avgD),
                coll,
            )
        end
    end

    if bc[4] != :fix
        @inbounds for i = 2:nx-1 # skip overlap
            fn(
                KS,
                ctr[i, ny],
                a1face[i, ny],
                a1face[i+1, ny],
                a2face[i, ny],
                a2face[i, ny+1],
                (dt, dx[i, ny] * dy[i, ny], resU, avgU),
                coll,
            )
        end
    end

    for i in eachindex(residual)
        residual[i] +=
            sqrt((resL[i] + resR[i] + resU[i] + resD[i]) * 2) /
            (avgL[i] + avgR[i] + avgU[i] + avgD[i] + 1.e-7)
    end

    ngx = 1 - first(eachindex(KS.ps.x[:, 1]))
    if bc[1] == :period
        bc_period!(ctr, ngx; dirc = :x)
    elseif bc[1] in (:extra, :mirror)
        bcfun = eval(Symbol("bc_" * string(bc[1]) * "!"))
        bcfun(ctr, ngx; dirc = :xl)
    end
    if bc[2] in (:extra, :mirror)
        bcfun = eval(Symbol("bc_" * string(bc[2]) * "!"))
        bcfun(ctr, ngx; dirc = :xr)
    end

    ngy = 1 - first(eachindex(KS.ps.y[1, :]))
    if bc[3] == :period
        bc_period!(ctr, ngy; dirc = :y)
    elseif bc[3] in (:extra, :mirror)
        bcfun = eval(Symbol("bc_" * string(bc[3]) * "!"))
        bcfun(ctr, ngy; dirc = :yl)
    end
    if bc[4] in (:extra, :mirror)
        bcfun = eval(Symbol("bc_" * string(bc[4]) * "!"))
        bcfun(ctr, ngy; dirc = :yr)
    end

    return nothing

end

function update_boundary!(
    KS,
    ctr::AV{TC},
    face,
    dt,
    residual;
    coll,
    bc,
    fn = step!,
) where {TC<:ControlVolumeUS}

    for i in eachindex(KS.ps.cellType)
        if KS.ps.cellType[i] == 3
            ids = KS.ps.cellNeighbors[i, :]
            deleteat!(ids, findall(x -> x == -1, ids))
            id1, id2 = ids
            ctr[i].w .= 0.5 .* (ctr[id1].w .+ ctr[id2].w)
            ctr[i].prim .= conserve_prim(ctr[i].w, KS.gas.γ)
        end
    end

    return nothing

end

function update_boundary!(
    KS,
    ctr::AV{TC},
    face,
    dt,
    residual;
    coll,
    bc,
    fn = step!,
) where {TC<:ControlVolumeUS1F}

    for i in eachindex(KS.ps.cellType)
        if KS.ps.cellType[i] == 3
            ids = KS.ps.cellNeighbors[i, :]
            deleteat!(ids, findall(x -> x == -1, ids))
            id1, id2 = ids
            ctr[i].w .= 0.5 .* (ctr[id1].w .+ ctr[id2].w)
            ctr[i].f .= 0.5 .* (ctr[id1].f .+ ctr[id2].f)
            ctr[i].prim .= conserve_prim(ctr[i].w, KS.gas.γ)
        end
    end

    return nothing

end

function update_boundary!(
    KS,
    ctr::AV{TC},
    face,
    dt,
    residual;
    coll,
    bc,
    fn = step!,
) where {TC<:ControlVolumeUS2F}

    for i in eachindex(KS.ps.cellType)
        if KS.ps.cellType[i] == 3
            ids = KS.ps.cellNeighbors[i, :]
            deleteat!(ids, findall(x -> x == -1, ids))
            id1, id2 = ids
            ctr[i].w .= 0.5 .* (ctr[id1].w .+ ctr[id2].w)
            ctr[i].h .= 0.5 .* (ctr[id1].h .+ ctr[id2].h)
            ctr[i].b .= 0.5 .* (ctr[id1].b .+ ctr[id2].b)
            ctr[i].prim .= conserve_prim(ctr[i].w, KS.gas.γ)
        end
    end

    return nothing

end
