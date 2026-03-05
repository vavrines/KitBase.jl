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
    fn=step!,
    st=fn,
    kwargs...,
) where {TC<:Union{
    ControlVolume,
    ControlVolume1D,
    ControlVolume1F,
    ControlVolume1D1F,
    ControlVolume2F,
    ControlVolume1D2F,
},}
    bcs = ifelse(bc isa Symbol, [bc, bc], bc)

    resL = zero(ctr[1].w)
    avgL = zero(ctr[1].w)
    resR = zero(ctr[1].w)
    avgR = zero(ctr[1].w)

    if bcs[1] != :fix
        i = 1
        fn(KS, ctr[i], face[i], face[i+1], (dt, KS.ps.dx[i], resL, avgL); st=st)
    end

    if bcs[2] != :fix
        j = KS.ps.nx
        fn(KS, ctr[j], face[j], face[j+1], (dt, KS.ps.dx[j], resR, avgR); st=st)
    end

    #@. residual += sqrt((resL + resR) * 2) / (avgL + avgR + 1.e-7)
    for i in eachindex(residual)
        residual[i] += sqrt((resL[i] + resR[i]) * 2) / (avgL[i] + avgR[i] + 1.e-7) # residual for first species only
    end

    ng = 1 - first(eachindex(KS.ps.x))
    if bcs[1] == :period
        bc_period!(ctr, ng)
    elseif bcs[1] == :extra
        bc_extra!(ctr, ng; dirc=:xl)
    elseif bcs[1] == :balance
        bc_balance!(ctr[0], ctr[1], ctr[2])
    end
    if bcs[2] == :extra
        bc_extra!(ctr, ng; dirc=:xr)
    elseif bcs[2] == :balance
        bc_balance!(ctr[KS.ps.nx+1], ctr[KS.ps.nx], ctr[KS.ps.nx-1])
    end

    return nothing
end

"""
$(SIGNATURES)
"""
function update_boundary!(
    KS,
    ctr::AV{TC},
    face,
    dt,
    residual;
    coll=symbolize(KS.set.collision),
    bc,
    fn=step!,
    st=fn,
    isMHD=false,
) where {TC<:Union{ControlVolume3F,ControlVolume1D3F}}
    bcs = ifelse(bc isa Symbol, [bc, bc], bc)

    resL = zero(ctr[1].w)
    avgL = zero(ctr[1].w)
    resR = zero(ctr[1].w)
    avgR = zero(ctr[1].w)

    if bcs[1] != :fix
        i = 1
        fn(KS, ctr[i], face[i], face[i+1], KS.ps.dx[i], dt, resL, avgL, coll, isMHD)
    end

    if bcs[2] != :fix
        j = KS.ps.nx
        fn(KS, ctr[j], face[j], face[j+1], KS.ps.dx[j], dt, resR, avgR, coll, isMHD)
    end

    for i in eachindex(residual)
        residual[i] += sqrt((resL[i] + resR[i]) * 2) / (avgL[i] + avgR[i] + 1.e-7)
    end

    ng = 1 - first(eachindex(KS.ps.x))
    if bcs[1] == :period
        bc_period!(ctr, ng)
    elseif bcs[1] == :extra
        bc_extra!(ctr, ng; dirc=:xl)
    elseif bcs[1] == :balance
        bc_balance!(ctr[0], ctr[1], ctr[2])
    end
    if bcs[2] == :extra
        bc_extra!(ctr, ng; dirc=:xr)
    elseif bcs[2] == :balance
        bc_balance!(ctr[KS.ps.nx+1], ctr[KS.ps.nx], ctr[KS.ps.nx-1])
    end

    return nothing
end

"""
$(SIGNATURES)
"""
function update_boundary!(
    KS,
    ctr::AV{TC},
    face,
    dt,
    residual;
    coll=symbolize(KS.set.collision),
    bc,
    fn=step!,
    st=fn,
    isMHD=false::Bool,
) where {TC<:Union{ControlVolume4F,ControlVolume1D4F}}
    bcs = ifelse(bc isa Symbol, [bc, bc], bc)

    resL = zero(ctr[1].w)
    avgL = zero(ctr[1].w)
    resR = zero(ctr[1].w)
    avgR = zero(ctr[1].w)

    if bcs[1] != :fix
        i = 1
        fn(KS, ctr[i], face[i], face[i+1], KS.ps.dx[i], dt, resL, avgL, coll, isMHD)
    end

    if bcs[2] != :fix
        j = KS.ps.nx
        fn(KS, ctr[j], face[j], face[j+1], KS.ps.dx[j], dt, resR, avgR, coll, isMHD)
    end

    for i in eachindex(residual)
        residual[i] += sqrt((resL[i] + resR[i]) * 2) / (avgL[i] + avgR[i] + 1.e-7)
    end

    ng = 1 - first(eachindex(KS.ps.x))
    if bcs[1] == :period
        bc_period!(ctr, ng)
    elseif bcs[1] == :extra
        bc_extra!(ctr, ng; dirc=:xl)
    elseif bcs[1] == :balance
        bc_balance!(ctr[0], ctr[1], ctr[2])
    end
    if bcs[2] == :extra
        bc_extra!(ctr, ng; dirc=:xr)
    elseif bcs[2] == :balance
        bc_balance!(ctr[KS.ps.nx+1], ctr[KS.ps.nx], ctr[KS.ps.nx-1])
    end

    return nothing
end

"""
$(SIGNATURES)
"""
function update_boundary!(
    KS,
    ctr::AM{TC},
    a1face,
    a2face,
    dt,
    residual;
    coll=symbolize(KS.set.collision),
    bc,
    fn=step!,
    st=fn,
    kwargs...,
) where {TC<:Union{
    ControlVolume,
    ControlVolume2D,
    ControlVolume1F,
    ControlVolume2D1F,
    ControlVolume2F,
    ControlVolume2D2F,
},}
    bcs = ifelse(bc isa Symbol, [bc, bc, bc, bc], bc)

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

    if bcs[1] != :fix
        @inbounds for j in 1:ny
            fn(
                KS,
                ctr[1, j],
                a1face[1, j],
                a1face[2, j],
                a2face[1, j],
                a2face[1, j+1],
                (dt, dx[1, j] * dy[1, j], resL, avgL),
                coll;
                st=st,
            )
        end
    end

    if bcs[2] != :fix
        @inbounds for j in 1:ny
            fn(
                KS,
                ctr[nx, j],
                a1face[nx, j],
                a1face[nx+1, j],
                a2face[nx, j],
                a2face[nx, j+1],
                (dt, dx[nx, j] * dy[nx, j], resR, avgR),
                coll;
                st=st,
            )
        end
    end

    if bcs[3] != :fix
        @inbounds for i in 2:nx-1 # skip overlap
            fn(
                KS,
                ctr[i, 1],
                a1face[i, 1],
                a1face[i+1, 1],
                a2face[i, 1],
                a2face[i, 2],
                (dt, dx[i, 1] * dy[i, 1], resD, avgD),
                coll;
                st=st,
            )
        end
    end

    if bcs[4] != :fix
        @inbounds for i in 2:nx-1 # skip overlap
            fn(
                KS,
                ctr[i, ny],
                a1face[i, ny],
                a1face[i+1, ny],
                a2face[i, ny],
                a2face[i, ny+1],
                (dt, dx[i, ny] * dy[i, ny], resU, avgU),
                coll;
                st=st,
            )
        end
    end

    for i in eachindex(residual)
        residual[i] +=
            sqrt((resL[i] + resR[i] + resU[i] + resD[i]) * 2) /
            (avgL[i] + avgR[i] + avgU[i] + avgD[i] + 1.e-7)
    end

    ngx = 1 - first(eachindex(KS.ps.x[:, 1]))
    if bcs[1] == :period
        bc_period!(ctr, ngx; dirc=:x)
    elseif bcs[1] in (:extra, :mirror)
        bcfun = eval(Symbol("bc_" * string(bcs[1]) * "!"))
        bcfun(ctr, ngx; dirc=:xl)
    end
    if bcs[2] in (:extra, :mirror)
        bcfun = eval(Symbol("bc_" * string(bcs[2]) * "!"))
        bcfun(ctr, ngx; dirc=:xr)
    end

    ngy = 1 - first(eachindex(KS.ps.y[1, :]))
    if bcs[3] == :period
        bc_period!(ctr, ngy; dirc=:y)
    elseif bcs[3] in (:extra, :mirror)
        bcfun = eval(Symbol("bc_" * string(bcs[3]) * "!"))
        bcfun(ctr, ngy; dirc=:yl)
    end
    if bcs[4] in (:extra, :mirror)
        bcfun = eval(Symbol("bc_" * string(bcs[4]) * "!"))
        bcfun(ctr, ngy; dirc=:yr)
    end

    return nothing
end

"""
$(SIGNATURES)
"""
function update_boundary!(
    KS,
    ctr::AA{TC,3},
    a1face,
    a2face,
    a3face,
    dt,
    residual;
    coll=symbolize(KS.set.collision),
    bc,
    fn=step!,
    st=fn,
    kwargs...,
) where {TC<:Union{
    ControlVolume,
    ControlVolume3D,
},}
    bcs = ifelse(bc isa Symbol, [bc, bc, bc, bc, bc, bc], bc)

    nx, ny, nz = KS.ps.nx, KS.ps.ny, KS.ps.nz
    dx, dy, dz = KS.ps.dx, KS.ps.dy, KS.ps.dz

    resXL = zero(ctr[1, 1, 1].w)
    avgXL = zero(ctr[1, 1, 1].w)
    resXR = zero(ctr[1, 1, 1].w)
    avgXR = zero(ctr[1, 1, 1].w)
    resYL = zero(ctr[1, 1, 1].w)
    avgYL = zero(ctr[1, 1, 1].w)
    resYR = zero(ctr[1, 1, 1].w)
    avgYR = zero(ctr[1, 1, 1].w)
    resZL = zero(ctr[1, 1, 1].w)
    avgZL = zero(ctr[1, 1, 1].w)
    resZR = zero(ctr[1, 1, 1].w)
    avgZR = zero(ctr[1, 1, 1].w)

    if bcs[1] != :fix
        @inbounds for k in 1:nz
            for j in 1:ny
                fn(
                    KS,
                    ctr[1, j, k],
                    a1face[1, j, k],
                    a1face[2, j, k],
                    a2face[1, j, k],
                    a2face[1, j+1, k],
                    a3face[1, j, k],
                    a3face[1, j, k+1],
                    (dt, dx[1, j, k] * dy[1, j, k] * dz[1, j, k], resXL, avgXL),
                    coll;
                    st=st,
                )
            end
        end
    end

    if bcs[2] != :fix
        @inbounds for k in 1:nz
            for j in 1:ny
                fn(
                    KS,
                    ctr[nx, j, k],
                    a1face[nx, j, k],
                    a1face[nx+1, j, k],
                    a2face[nx, j, k],
                    a2face[nx, j+1, k],
                    a3face[nx, j, k],
                    a3face[nx, j, k+1],
                    (dt, dx[nx, j, k] * dy[nx, j, k] * dz[nx, j, k], resXR, avgXR),
                    coll;
                    st=st,
                )
            end
        end
    end

    if bcs[3] != :fix
        @inbounds for k in 1:nz
            for i in 2:nx-1  # skip overlap with x boundaries
                fn(
                    KS,
                    ctr[i, 1, k],
                    a1face[i, 1, k],
                    a1face[i+1, 1, k],
                    a2face[i, 1, k],
                    a2face[i, 2, k],
                    a3face[i, 1, k],
                    a3face[i, 1, k+1],
                    (dt, dx[i, 1, k] * dy[i, 1, k] * dz[i, 1, k], resYL, avgYL),
                    coll;
                    st=st,
                )
            end
        end
    end

    if bcs[4] != :fix
        @inbounds for k in 1:nz
            for i in 2:nx-1  # skip overlap with x boundaries
                fn(
                    KS,
                    ctr[i, ny, k],
                    a1face[i, ny, k],
                    a1face[i+1, ny, k],
                    a2face[i, ny, k],
                    a2face[i, ny+1, k],
                    a3face[i, ny, k],
                    a3face[i, ny, k+1],
                    (dt, dx[i, ny, k] * dy[i, ny, k] * dz[i, ny, k], resYR, avgYR),
                    coll;
                    st=st,
                )
            end
        end
    end

    if bcs[5] != :fix
        @inbounds for j in 2:ny-1  # skip overlap with x and y boundaries
            for i in 2:nx-1
                fn(
                    KS,
                    ctr[i, j, 1],
                    a1face[i, j, 1],
                    a1face[i+1, j, 1],
                    a2face[i, j, 1],
                    a2face[i, j+1, 1],
                    a3face[i, j, 1],
                    a3face[i, j, 2],
                    (dt, dx[i, j, 1] * dy[i, j, 1] * dz[i, j, 1], resZL, avgZL),
                    coll;
                    st=st,
                )
            end
        end
    end

    if bcs[6] != :fix
        @inbounds for j in 2:ny-1  # skip overlap with x and y boundaries
            for i in 2:nx-1
                fn(
                    KS,
                    ctr[i, j, nz],
                    a1face[i, j, nz],
                    a1face[i+1, j, nz],
                    a2face[i, j, nz],
                    a2face[i, j+1, nz],
                    a3face[i, j, nz],
                    a3face[i, j, nz+1],
                    (dt, dx[i, j, nz] * dy[i, j, nz] * dz[i, j, nz], resZR, avgZR),
                    coll;
                    st=st,
                )
            end
        end
    end

    for i in eachindex(residual)
        residual[i] +=
            sqrt((resXL[i] + resXR[i] + resYL[i] + resYR[i] + resZL[i] + resZR[i]) * 2) /
            (avgXL[i] + avgXR[i] + avgYL[i] + avgYR[i] + avgZL[i] + avgZR[i] + 1.e-7)
    end

    ngx = 1 - first(eachindex(KS.ps.x[:, 1, 1]))
    if bcs[1] == :period
        bc_period!(ctr, ngx; dirc=:x)
    elseif bcs[1] in (:extra, :mirror)
        bcfun = eval(Symbol("bc_" * string(bcs[1]) * "!"))
        bcfun(ctr, ngx; dirc=:xl)
    end
    if bcs[2] in (:extra, :mirror)
        bcfun = eval(Symbol("bc_" * string(bcs[2]) * "!"))
        bcfun(ctr, ngx; dirc=:xr)
    end

    ngy = 1 - first(eachindex(KS.ps.y[1, :, 1]))
    if bcs[3] == :period
        bc_period!(ctr, ngy; dirc=:y)
    elseif bcs[3] in (:extra, :mirror)
        bcfun = eval(Symbol("bc_" * string(bcs[3]) * "!"))
        bcfun(ctr, ngy; dirc=:yl)
    end
    if bcs[4] in (:extra, :mirror)
        bcfun = eval(Symbol("bc_" * string(bcs[4]) * "!"))
        bcfun(ctr, ngy; dirc=:yr)
    end

    ngz = 1 - first(eachindex(KS.ps.z[1, 1, :]))
    if bcs[5] == :period
        bc_period!(ctr, ngz; dirc=:z)
    elseif bcs[5] in (:extra, :mirror)
        bcfun = eval(Symbol("bc_" * string(bcs[5]) * "!"))
        bcfun(ctr, ngz; dirc=:zl)
    end
    if bcs[6] in (:extra, :mirror)
        bcfun = eval(Symbol("bc_" * string(bcs[6]) * "!"))
        bcfun(ctr, ngz; dirc=:zr)
    end

    return nothing
end

"""
$(SIGNATURES)
"""
function update_boundary!(
    KS,
    ctr::AV{TC},
    face,
    dt,
    residual;
    coll,
    bc,
    fn=step!,
    st=fn,
    kwargs...,
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

"""
$(SIGNATURES)
"""
function update_boundary!(
    KS,
    ctr::AV{TC},
    face,
    dt,
    residual;
    coll,
    bc,
    fn=step!,
    st=fn,
    kwargs...,
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

"""
$(SIGNATURES)
"""
function update_boundary!(
    KS,
    ctr::AV{TC},
    face,
    dt,
    residual;
    coll,
    bc,
    fn=step!,
    st=fn,
    kwargs...,
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
