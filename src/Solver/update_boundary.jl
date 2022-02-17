"""
$(SIGNATURES)

Update solver for boundary cells
"""
function update_boundary!(
    KS::AbstractSolverSet,
    ctr::AV{TC},
    face::AV{TF},
    dt,
    residual;
    bc,
) where {TC<:Union{ControlVolume,ControlVolume1D},TF<:Union{Interface,Interface1D}}

    resL = zero(ctr[1].w)
    avgL = zero(ctr[1].w)
    resR = zero(ctr[1].w)
    avgR = zero(ctr[1].w)

    if bc[1] != :fix
        i = 1
        j = KS.pSpace.nx
        if KS.set.nSpecies == 1
            step!(
                face[i].fw,
                ctr[i].w,
                ctr[i].prim,
                face[i+1].fw,
                KS.gas.γ,
                KS.ps.dx[i],
                resL,
                avgL,
            )
        elseif KS.set.nSpecies == 2
            step!(
                face[i].fw,
                ctr[i].w,
                ctr[i].prim,
                face[i+1].fw,
                KS.gas.γ,
                KS.gas.mi,
                KS.gas.ni,
                KS.gas.me,
                KS.gas.ne,
                KS.gas.Kn[1],
                KS.ps.dx[i],
                dt,
                resL,
                avgL,
            )
        end
    end

    if bc[2] != :fix
        j = KS.pSpace.nx
        if KS.set.nSpecies == 1
            step!(
                face[j].fw,
                ctr[j].w,
                ctr[j].prim,
                face[j+1].fw,
                KS.gas.γ,
                KS.ps.dx[j],
                resR,
                avgR,
            )
        elseif KS.set.nSpecies == 2
            step!(
                face[j].fw,
                ctr[j].w,
                ctr[j].prim,
                face[j+1].fw,
                KS.gas.γ,
                KS.gas.mi,
                KS.gas.ni,
                KS.gas.me,
                KS.gas.ne,
                KS.gas.Kn[1],
                KS.ps.dx[j],
                dt,
                resR,
                avgR,
            )
        end
    end

    @. residual += sqrt((resL + resR) * 2) / (avgL + avgR + 1.e-7)

    ng = 1 - first(eachindex(KS.pSpace.x))
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
    KS::AbstractSolverSet,
    ctr::AV{TC},
    face::AV{TF},
    dt,
    residual;
    coll = symbolize(KS.set.collision),
    bc,
) where {TC<:Union{ControlVolume1F,ControlVolume1D1F},TF<:Union{Interface1F,Interface1D1F}}

    resL = zero(ctr[1].w)
    avgL = zero(ctr[1].w)
    resR = zero(ctr[1].w)
    avgR = zero(ctr[1].w)

    if bc[1] != :fix
        i = 1

        if KS.set.space[5:6] == "1v"
            step!(
                face[i].fw,
                face[i].ff,
                ctr[i].w,
                ctr[i].prim,
                ctr[i].f,
                face[i+1].fw,
                face[i+1].ff,
                KS.vSpace.u,
                KS.vSpace.weights,
                KS.gas.γ,
                KS.gas.μᵣ,
                KS.gas.ω,
                KS.gas.Pr,
                KS.ps.dx[i],
                dt,
                resL,
                avgL,
                coll,
            )
        elseif KS.set.space[5:6] == "3v"
            step!(
                face[i].fw,
                face[i].ff,
                ctr[i].w,
                ctr[i].prim,
                ctr[i].f,
                face[i+1].fw,
                face[i+1].ff,
                KS.vSpace.u,
                KS.vSpace.v,
                KS.vSpace.w,
                KS.vSpace.weights,
                KS.gas.γ,
                KS.gas.μᵣ,
                KS.gas.ω,
                KS.gas.Pr,
                KS.ps.dx[i],
                dt,
                resL,
                avgL,
                coll,
            )
        end
    end

    if bc[2] != :fix
        j = KS.pSpace.nx

        if KS.set.space[5:6] == "1v"
            step!(
                face[j].fw,
                face[j].ff,
                ctr[j].w,
                ctr[j].prim,
                ctr[j].f,
                face[j+1].fw,
                face[j+1].ff,
                KS.vSpace.u,
                KS.vSpace.weights,
                KS.gas.γ,
                KS.gas.μᵣ,
                KS.gas.ω,
                KS.gas.Pr,
                KS.ps.dx[j],
                dt,
                resR,
                avgR,
                coll,
            )
        elseif KS.set.space[5:6] == "3v"
            step!(
                face[j].fw,
                face[j].ff,
                ctr[j].w,
                ctr[j].prim,
                ctr[j].f,
                face[j+1].fw,
                face[j+1].ff,
                KS.vSpace.u,
                KS.vSpace.v,
                KS.vSpace.w,
                KS.vSpace.weights,
                KS.gas.γ,
                KS.gas.μᵣ,
                KS.gas.ω,
                KS.gas.Pr,
                KS.ps.dx[j],
                dt,
                resL,
                avgL,
                coll,
            )
        end
    end

    for i in eachindex(residual)
        residual[i] += sqrt((resL[i] + resR[i]) * 2) / (avgL[i] + avgR[i] + 1.e-7)
    end

    ng = 1 - first(eachindex(KS.pSpace.x))
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
    KS::AbstractSolverSet,
    ctr::AV{TC},
    face::AV{TF},
    dt,
    residual;
    coll = symbolize(KS.set.collision),
    bc,
) where {TC<:Union{ControlVolume2F,ControlVolume1D2F},TF<:Union{Interface2F,Interface1D2F}}

    resL = zero(ctr[1].w)
    avgL = zero(ctr[1].w)
    resR = zero(ctr[1].w)
    avgR = zero(ctr[1].w)

    if bc[1] != :fix
        i = 1
        if KS.set.nSpecies == 1
            step!(
                face[i].fw,
                face[i].fh,
                face[i].fb,
                ctr[i].w,
                ctr[i].prim,
                ctr[i].h,
                ctr[i].b,
                face[i+1].fw,
                face[i+1].fh,
                face[i+1].fb,
                KS.vSpace.u,
                KS.vSpace.weights,
                KS.gas.K,
                KS.gas.γ,
                KS.gas.μᵣ,
                KS.gas.ω,
                KS.gas.Pr,
                KS.ps.dx[i],
                dt,
                resL,
                avgL,
                coll,
            )
        else
            step!(
                face[i].fw,
                face[i].fh,
                face[i].fb,
                ctr[i].w,
                ctr[i].prim,
                ctr[i].h,
                ctr[i].b,
                face[i+1].fw,
                face[i+1].fh,
                face[i+1].fb,
                KS.vSpace.u,
                KS.vSpace.weights,
                KS.gas.K,
                KS.gas.γ,
                KS.gas.mi,
                KS.gas.ni,
                KS.gas.me,
                KS.gas.ne,
                KS.gas.Kn[1],
                KS.gas.Pr,
                KS.ps.dx[i],
                dt,
                resL,
                avgL,
                coll,
            )
        end
    end

    if bc[2] != :fix
        j = KS.pSpace.nx
        if KS.set.nSpecies == 1
            step!(
                face[j].fw,
                face[j].fh,
                face[j].fb,
                ctr[j].w,
                ctr[j].prim,
                ctr[j].h,
                ctr[j].b,
                face[j+1].fw,
                face[j+1].fh,
                face[j+1].fb,
                KS.vSpace.u,
                KS.vSpace.weights,
                KS.gas.K,
                KS.gas.γ,
                KS.gas.μᵣ,
                KS.gas.ω,
                KS.gas.Pr,
                KS.ps.dx[j],
                dt,
                resL,
                avgL,
                coll,
            )
        else
            step!(
                face[j].fw,
                face[j].fh,
                face[j].fb,
                ctr[j].w,
                ctr[j].prim,
                ctr[j].h,
                ctr[j].b,
                face[j+1].fw,
                face[j+1].fh,
                face[j+1].fb,
                KS.vSpace.u,
                KS.vSpace.weights,
                KS.gas.K,
                KS.gas.γ,
                KS.gas.mi,
                KS.gas.ni,
                KS.gas.me,
                KS.gas.ne,
                KS.gas.Kn[1],
                KS.gas.Pr,
                KS.ps.dx[j],
                dt,
                resL,
                avgL,
                coll,
            )
        end
    end

    for i in eachindex(residual)
        residual[i] += sqrt((resL[i] + resR[i]) * 2) / (avgL[i] + avgR[i] + 1.e-7)
    end

    ng = 1 - first(eachindex(KS.pSpace.x))
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
    KS::AbstractSolverSet,
    ctr::AV{TC},
    face::AV{TF},
    dt,
    residual;
    coll = symbolize(KS.set.collision)::Symbol,
    bc,
    isMHD = false::Bool,
) where {TC<:Union{ControlVolume3F,ControlVolume1D3F},TF<:Union{Interface3F,Interface1D3F}}

    resL = zero(ctr[1].w)
    avgL = zero(ctr[1].w)
    resR = zero(ctr[1].w)
    avgR = zero(ctr[1].w)

    if bc[1] != :fix
        i = 1
        step!(KS, face[i], ctr[i], face[i+1], KS.ps.dx[i], dt, resL, avgL, coll, isMHD)
    end

    if bc[2] != :fix
        j = KS.pSpace.nx
        step!(KS, face[j], ctr[j], face[j+1], KS.ps.dx[j], dt, resR, avgR, coll, isMHD)
    end

    for i in eachindex(residual)
        residual[i] += sqrt((resL[i] + resR[i]) * 2) / (avgL[i] + avgR[i] + 1.e-7)
    end

    ng = 1 - first(eachindex(KS.pSpace.x))
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
    KS::AbstractSolverSet,
    ctr::AV{TC},
    face::AV{TF},
    dt,
    residual;
    coll = symbolize(KS.set.collision)::Symbol,
    bc,
    isMHD = false::Bool,
) where {TC<:Union{ControlVolume4F,ControlVolume1D4F},TF<:Union{Interface4F,Interface1D4F}}

    resL = zero(ctr[1].w)
    avgL = zero(ctr[1].w)
    resR = zero(ctr[1].w)
    avgR = zero(ctr[1].w)

    if bc[1] != :fix
        i = 1
        step!(KS, face[i], ctr[i], face[i+1], KS.ps.dx[i], dt, resL, avgL, coll, isMHD)
    end

    if bc[2] != :fix
        j = KS.pSpace.nx
        step!(KS, face[j], ctr[j], face[j+1], KS.ps.dx[j], dt, resR, avgR, coll, isMHD)
    end

    for i in eachindex(residual)
        residual[i] += sqrt((resL[i] + resR[i]) * 2) / (avgL[i] + avgR[i] + 1.e-7)
    end

    ng = 1 - first(eachindex(KS.pSpace.x))
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
    KS::AbstractSolverSet,
    ctr::AM{TC},
    a1face::AM{TF},
    a2face::AM{TF},
    dt,
    residual;
    coll = symbolize(KS.set.collision)::Symbol,
    bc,
) where {TC<:Union{ControlVolume,ControlVolume2D},TF<:Union{Interface,Interface2D}}

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
            step!(
                ctr[1, j].w,
                ctr[1, j].prim,
                a1face[1, j].fw,
                a1face[2, j].fw,
                a2face[1, j].fw,
                a2face[1, j+1].fw,
                KS.gas.γ,
                dx[1, j] * dy[1, j],
                resL,
                avgL,
                coll,
            )
        end
    end

    if bc[2] != :fix
        @inbounds for j = 1:ny
            step!(
                ctr[nx, j].w,
                ctr[nx, j].prim,
                a1face[nx, j].fw,
                a1face[nx+1, j].fw,
                a2face[nx, j].fw,
                a2face[nx, j+1].fw,
                KS.gas.γ,
                dx[nx, j] * dy[nx, j],
                resR,
                avgR,
                coll,
            )
        end
    end

    if bc[3] != :fix
        @inbounds for i = 2:nx-1 # skip overlap
            step!(
                ctr[i, 1].w,
                ctr[i, 1].prim,
                a1face[i, 1].fw,
                a1face[i+1, 1].fw,
                a2face[i, 1].fw,
                a2face[i, 2].fw,
                KS.gas.γ,
                dx[i, 1] * dy[i, 1],
                resD,
                avgD,
                coll,
            )
        end
    end

    if bc[4] != :fix
        @inbounds for i = 2:nx-1 # skip overlap
            step!(
                ctr[i, ny].w,
                ctr[i, ny].prim,
                a1face[i, ny].fw,
                a1face[i+1, ny].fw,
                a2face[i, ny].fw,
                a2face[i, ny+1].fw,
                KS.gas.γ,
                dx[i, ny] * dy[i, ny],
                resU,
                avgU,
                coll,
            )
        end
    end

    for i in eachindex(residual)
        residual[i] +=
            sqrt((resL[i] + resR[i] + resU[i] + resD[i]) * 2) /
            (avgL[i] + avgR[i] + avgU[i] + avgD[i] + 1.e-7)
    end

    ngx = 1 - first(eachindex(KS.pSpace.x[:, 1]))
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

    ngy = 1 - first(eachindex(KS.pSpace.y[1, :]))
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
    KS::AbstractSolverSet,
    ctr::AM{TC},
    a1face::AM{TF},
    a2face::AM{TF},
    dt,
    residual;
    coll = symbolize(KS.set.collision)::Symbol,
    bc,
) where {TC<:Union{ControlVolume1F,ControlVolume2D1F},TF<:Union{Interface1F,Interface2D1F}}

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
            step!(
                ctr[1, j].w,
                ctr[1, j].prim,
                ctr[1, j].f,
                a1face[1, j].fw,
                a1face[1, j].ff,
                a1face[2, j].fw,
                a1face[2, j].ff,
                a2face[1, j].fw,
                a2face[1, j].ff,
                a2face[1, j+1].fw,
                a2face[1, j+1].ff,
                KS.vSpace.u,
                KS.vSpace.v,
                KS.vSpace.weights,
                KS.gas.γ,
                KS.gas.μᵣ,
                KS.gas.ω,
                KS.gas.Pr,
                dx[1, j] * dy[1, j],
                dt,
                resL,
                avgL,
                coll,
            )
        end
    end
    if bc[2] != :fix
        @inbounds for j = 1:ny
            step!(
                ctr[nx, j].w,
                ctr[nx, j].prim,
                ctr[nx, j].f,
                a1face[nx, j].fw,
                a1face[nx, j].ff,
                a1face[nx+1, j].fw,
                a1face[nx+1, j].ff,
                a2face[nx, j].fw,
                a2face[nx, j].ff,
                a2face[nx, j+1].fw,
                a2face[nx, j+1].ff,
                KS.vSpace.u,
                KS.vSpace.v,
                KS.vSpace.weights,
                KS.gas.γ,
                KS.gas.μᵣ,
                KS.gas.ω,
                KS.gas.Pr,
                dx[nx, j] * dy[nx, j],
                dt,
                resR,
                avgR,
                coll,
            )
        end
    end
    if bc[3] != :fix
        @inbounds for i = 2:nx-1 # skip overlap
            step!(
                ctr[i, 1].w,
                ctr[i, 1].prim,
                ctr[i, 1].f,
                a1face[i, 1].fw,
                a1face[i, 1].ff,
                a1face[i+1, 1].fw,
                a1face[i+1, 1].ff,
                a2face[i, 1].fw,
                a2face[i, 1].ff,
                a2face[i, 2].fw,
                a2face[i, 2].ff,
                KS.vSpace.u,
                KS.vSpace.v,
                KS.vSpace.weights,
                KS.gas.γ,
                KS.gas.μᵣ,
                KS.gas.ω,
                KS.gas.Pr,
                dx[i, 1] * dy[i, 1],
                dt,
                resD,
                avgD,
                coll,
            )
        end
    end
    if bc[4] != :fix
        @inbounds for i = 2:nx-1 # skip overlap
            step!(
                ctr[i, ny].w,
                ctr[i, ny].prim,
                ctr[i, ny].f,
                a1face[i, ny].fw,
                a1face[i, ny].ff,
                a1face[i+1, ny].fw,
                a1face[i+1, ny].ff,
                a2face[i, ny].fw,
                a2face[i, ny].ff,
                a2face[i, ny+1].fw,
                a2face[i, ny+1].ff,
                KS.vSpace.u,
                KS.vSpace.v,
                KS.vSpace.weights,
                KS.gas.γ,
                KS.gas.μᵣ,
                KS.gas.ω,
                KS.gas.Pr,
                dx[i, ny] * dy[i, ny],
                dt,
                resU,
                avgU,
                coll,
            )
        end
    end

    for i in eachindex(residual)
        residual[i] +=
            sqrt((resL[i] + resR[i] + resU[i] + resD[i]) * 2) /
            (avgL[i] + avgR[i] + avgU[i] + avgD[i] + 1.e-7)
    end

    ngx = 1 - first(eachindex(KS.pSpace.x[:, 1]))
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

    ngy = 1 - first(eachindex(KS.pSpace.y[1, :]))
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
    KS::AbstractSolverSet,
    ctr::AM{TC},
    a1face::AM{TF},
    a2face::AM{TF},
    dt,
    residual;
    coll = symbolize(KS.set.collision)::Symbol,
    bc,
) where {TC<:Union{ControlVolume2F,ControlVolume2D2F},TF<:Union{Interface2F,Interface2D2F}}

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
            step!(
                ctr[1, j].w,
                ctr[1, j].prim,
                ctr[1, j].h,
                ctr[1, j].b,
                a1face[1, j].fw,
                a1face[1, j].fh,
                a1face[1, j].fb,
                a1face[2, j].fw,
                a1face[2, j].fh,
                a1face[2, j].fb,
                a2face[1, j].fw,
                a2face[1, j].fh,
                a2face[1, j].fb,
                a2face[1, j+1].fw,
                a2face[1, j+1].fh,
                a2face[1, j+1].fb,
                KS.vSpace.u,
                KS.vSpace.v,
                KS.vSpace.weights,
                KS.gas.K,
                KS.gas.γ,
                KS.gas.μᵣ,
                KS.gas.ω,
                KS.gas.Pr,
                dx[1, j] * dy[1, j],
                dt,
                resL,
                avgL,
                coll,
            )
        end
    end
    if bc[2] != :fix
        @inbounds for j = 1:ny
            step!(
                ctr[nx, j].w,
                ctr[nx, j].prim,
                ctr[nx, j].h,
                ctr[nx, j].b,
                a1face[nx, j].fw,
                a1face[nx, j].fh,
                a1face[nx, j].fb,
                a1face[nx+1, j].fw,
                a1face[nx+1, j].fh,
                a1face[nx+1, j].fb,
                a2face[nx, j].fw,
                a2face[nx, j].fh,
                a2face[nx, j].fb,
                a2face[nx, j+1].fw,
                a2face[nx, j+1].fh,
                a2face[nx, j+1].fb,
                KS.vSpace.u,
                KS.vSpace.v,
                KS.vSpace.weights,
                KS.gas.K,
                KS.gas.γ,
                KS.gas.μᵣ,
                KS.gas.ω,
                KS.gas.Pr,
                dx[nx, j] * dy[nx, j],
                dt,
                resR,
                avgR,
                coll,
            )
        end
    end
    if bc[3] != :fix
        @inbounds for i = 2:nx-1 # skip overlap
            step!(
                ctr[i, 1].w,
                ctr[i, 1].prim,
                ctr[i, 1].h,
                ctr[i, 1].b,
                a1face[i, 1].fw,
                a1face[i, 1].fh,
                a1face[i, 1].fb,
                a1face[i+1, 1].fw,
                a1face[i+1, 1].fh,
                a1face[i+1, 1].fb,
                a2face[i, 1].fw,
                a2face[i, 1].fh,
                a2face[i, 1].fb,
                a2face[i, 2].fw,
                a2face[i, 2].fh,
                a2face[i, 2].fb,
                KS.vSpace.u,
                KS.vSpace.v,
                KS.vSpace.weights,
                KS.gas.K,
                KS.gas.γ,
                KS.gas.μᵣ,
                KS.gas.ω,
                KS.gas.Pr,
                dx[i, 1] * dy[i, 1],
                dt,
                resD,
                avgD,
                coll,
            )
        end
    end
    if bc[4] != :fix
        @inbounds for i = 2:nx-1 # skip overlap
            step!(
                ctr[i, ny].w,
                ctr[i, ny].prim,
                ctr[i, ny].h,
                ctr[i, ny].b,
                a1face[i, ny].fw,
                a1face[i, ny].fh,
                a1face[i, ny].fb,
                a1face[i+1, ny].fw,
                a1face[i+1, ny].fh,
                a1face[i+1, ny].fb,
                a2face[i, ny].fw,
                a2face[i, ny].fh,
                a2face[i, ny].fb,
                a2face[i, ny+1].fw,
                a2face[i, ny+1].fh,
                a2face[i, ny+1].fb,
                KS.vSpace.u,
                KS.vSpace.v,
                KS.vSpace.weights,
                KS.gas.K,
                KS.gas.γ,
                KS.gas.μᵣ,
                KS.gas.ω,
                KS.gas.Pr,
                dx[i, ny] * dy[i, ny],
                dt,
                resU,
                avgU,
                coll,
            )
        end
    end

    for i in eachindex(residual)
        residual[i] +=
            sqrt((resL[i] + resR[i] + resU[i] + resD[i]) * 2) /
            (avgL[i] + avgR[i] + avgU[i] + avgD[i] + 1.e-7)
    end

    ngx = 1 - first(eachindex(KS.pSpace.x[:, 1]))
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

    ngy = 1 - first(eachindex(KS.pSpace.y[1, :]))
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
    KS::AbstractSolverSet,
    ctr::AV{TC},
    face::AV{TF},
    dt,
    residual;
    coll::Symbol,
    bc,
) where {TC<:ControlVolumeUS,TF<:Interface2D}

    for i in eachindex(KS.ps.cellType)
        if KS.ps.cellType[i] == 3
            ids = KS.ps.cellNeighbors[i, :]
            deleteat!(ids, findall(x -> x == -1, ids))
            id1, id2 = ids
            ctr[i].w .= 0.5 .* (ctr[id1].w .+ ctr[id2].w)
            ctr[i].prim .= KitBase.conserve_prim(ctr[i].w, KS.gas.γ)
        end
    end

    return nothing

end

function update_boundary!(
    KS::AbstractSolverSet,
    ctr::AV{TC},
    face::AV{TF},
    dt,
    residual;
    coll::Symbol,
    bc,
) where {TC<:ControlVolumeUS1F,TF<:Interface2D1F}

    for i in eachindex(KS.ps.cellType)
        if KS.ps.cellType[i] == 3
            ids = KS.ps.cellNeighbors[i, :]
            deleteat!(ids, findall(x -> x == -1, ids))
            id1, id2 = ids
            ctr[i].w .= 0.5 .* (ctr[id1].w .+ ctr[id2].w)
            ctr[i].f .= 0.5 .* (ctr[id1].f .+ ctr[id2].f)
            ctr[i].prim .= KitBase.conserve_prim(ctr[i].w, KS.gas.γ)
        end
    end

    return nothing

end

function update_boundary!(
    KS::AbstractSolverSet,
    ctr::AV{TC},
    face::AV{TF},
    dt,
    residual;
    coll::Symbol,
    bc,
) where {TC<:ControlVolumeUS2F,TF<:Interface2D2F}

    for i in eachindex(KS.ps.cellType)
        if KS.ps.cellType[i] == 3
            ids = KS.ps.cellNeighbors[i, :]
            deleteat!(ids, findall(x -> x == -1, ids))
            id1, id2 = ids
            ctr[i].w .= 0.5 .* (ctr[id1].w .+ ctr[id2].w)
            ctr[i].h .= 0.5 .* (ctr[id1].h .+ ctr[id2].h)
            ctr[i].b .= 0.5 .* (ctr[id1].b .+ ctr[id2].b)
            ctr[i].prim .= KitBase.conserve_prim(ctr[i].w, KS.gas.γ)
        end
    end

    return nothing

end
