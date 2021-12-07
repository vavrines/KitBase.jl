"""
    update!(
        KS::X,
        ctr::Y,
        face::Z,
        dt,
        residual; # 1D / 2D
        coll = :bgk::Symbol,
        bc = :fix::Symbol,
    ) where {
        X<:AbstractSolverSet,
        Y<:AA{ControlVolume1D1F,1},
        Z<:AA{Interface1D1F,1},
    }

    update!(
        KS::X,
        ctr::Y,
        face::Z,
        dt,
        residual; # 1D / 2D
        coll = :bgk::Symbol,
        bc = :extra::Symbol,
    ) where {
        X<:AbstractSolverSet,
        Y<:AA{ControlVolume1D2F,1},
        Z<:AA{Interface1D2F,1},
    }

    update!(
        KS::X,
        ctr::Y,
        face::Z,
        dt,
        residual; # 1D / 2D
        coll = :bgk::Symbol,
        bc = :extra::Symbol,
        isMHD = true::Bool,
    ) where {
        X<:AbstractSolverSet,
        Y<:AA{ControlVolume1D3F,1},
        Z<:AA{Interface1D3F,1},
    }

    update!(
        KS::X,
        ctr::Y,
        face::Z,
        dt,
        residual; # 1D / 2D
        coll = :bgk::Symbol,
        bc = :extra::Symbol,
        isMHD = true::Bool,
    ) where {
        X<:AbstractSolverSet,
        Y<:AA{ControlVolume1D4F,1},
        Z<:AA{Interface1D4F,1},
    }

    update!(
        KS::X,
        ctr::Y,
        a1face::Z,
        a2face::Z,
        dt,
        residual; # 1D / 2D
        coll = :bgk::Symbol,
        bc = :fix::Symbol,
    ) where {
        X<:AbstractSolverSet,
        Y<:AA{ControlVolume2D2F,2},
        Z<:AA{Interface2D2F,2},
    }

Update flow variables over control volumes

"""
function update!(
    KS::X,
    ctr::Y,
    face::Z,
    dt,
    residual::T; # 1D / 2D
    coll = symbolize(KS.set.collision),
    bc = symbolize(KS.set.boundary),
) where {
    X<:AbstractSolverSet,
    Y<:AA{ControlVolume1D,1},
    Z<:AA{Interface1D,1},
    T<:AbstractFloat,
}

    sumRes = zero(ctr[1].w)
    sumAvg = zero(ctr[1].w)

    @inbounds Threads.@threads for i = 1:KS.pSpace.nx
        ctr[i].w, sumRes, sumAvg = step!(
            face[i].fw,
            ctr[i].w,
            ctr[i].prim,
            face[i+1].fw,
            KS.gas.a,
            KS.ps.dx[i],
            sumRes,
            sumAvg,
        )
    end

    residual = sqrt(sumRes * KS.pSpace.nx) / (sumAvg + 1.e-7)

    if bc[1] == :period
        ctr[0].w = ctr[KS.pSpace.nx].w
        ctr[0].sw = ctr[KS.pSpace.nx].sw
        ctr[KS.pSpace.nx+1].w = ctr[1].w
        ctr[KS.pSpace.nx+1].sw = ctr[1].sw
    elseif bc[1] == :extra
        ctr[0].w = ctr[1].w
    end

    if bc[2] == :extra
        ctr[KS.pSpace.nx+1].w = ctr[KS.pSpace.nx].w
    end

    return residual

end

function update!(
    KS::X,
    ctr::Y,
    face::Z,
    dt,
    residual::AA{T}; # 1D / 2D
    coll = symbolize(KS.set.collision),
    bc = symbolize(KS.set.boundary),
) where {
    X<:AbstractSolverSet,
    Y<:AA{ControlVolume1D,1},
    Z<:AA{Interface1D,1},
    T<:AbstractFloat,
}

    sumRes = zero(ctr[1].w)
    sumAvg = zero(ctr[1].w)

    if ndims(sumRes) == 1
        @inbounds Threads.@threads for i = 2:KS.pSpace.nx-1
            step!(
                face[i].fw,
                ctr[i].w,
                ctr[i].prim,
                face[i+1].fw,
                KS.gas.γ,
                KS.ps.dx[i],
                sumRes,
                sumAvg,
            )
        end
    elseif ndims(sumRes) == 2
        @inbounds Threads.@threads for i = 2:KS.pSpace.nx-1
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
                sumRes,
                sumAvg,
            )
        end
    end

    for i in eachindex(residual)
        residual[i] = sqrt(sumRes[i] * KS.pSpace.nx) / (sumAvg[i] + 1.e-7)
    end

    if bc isa Symbol
        update_boundary!(KS, ctr, face, dt, residual; bc = [bc, bc])
    else
        update_boundary!(KS, ctr, face, dt, residual; bc = bc)
    end

    return nothing

end

function update!(
    KS::X,
    ctr::Y,
    face::Z,
    dt,
    residual; # 1D / 2D
    coll = symbolize(KS.set.collision),
    bc = symbolize(KS.set.boundary),
) where {
    X<:AbstractSolverSet,
    Y<:AA{ControlVolume1D1F,1},
    Z<:AA{Interface1D1F,1},
}

    sumRes = zero(ctr[1].w)
    sumAvg = zero(ctr[1].w)

    phase = KS.set.space[3:end]
    if phase == "1f1v"
        @inbounds Threads.@threads for i = 2:KS.pSpace.nx-1
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
                sumRes,
                sumAvg,
                coll,
            )
        end
    elseif phase == "1f3v"
        @inbounds Threads.@threads for i = 2:KS.pSpace.nx-1
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
                sumRes,
                sumAvg,
                coll,
            )
        end
    end

    for i in eachindex(residual)
        residual[i] = sqrt(sumRes[i] * KS.pSpace.nx) / (sumAvg[i] + 1.e-7)
    end

    if bc isa Symbol
        update_boundary!(KS, ctr, face, dt, residual; coll = coll, bc = [bc, bc])
    else
        update_boundary!(KS, ctr, face, dt, residual; coll = coll, bc = bc)
    end

    return nothing

end

function update!(
    KS::X,
    ctr::Y,
    face::Z,
    dt,
    residual; # 1D / 2D
    coll = symbolize(KS.set.collision),
    bc = symbolize(KS.set.boundary),
) where {
    X<:AbstractSolverSet,
    Y<:AA{ControlVolume1D2F,1},
    Z<:AA{Interface1D2F,1},
}

    sumRes = zero(ctr[1].w)
    sumAvg = zero(ctr[1].w)

    if ndims(sumRes) == 1
        @inbounds Threads.@threads for i = 2:KS.pSpace.nx-1
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
                sumRes,
                sumAvg,
                coll,
            )
        end
    elseif ndims(sumRes) == 2
        @inbounds Threads.@threads for i = 2:KS.pSpace.nx-1
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
                sumRes,
                sumAvg,
                coll,
            )
        end
    end

    for i in eachindex(residual)
        residual[i] = sqrt(sumRes[i] * KS.pSpace.nx) / (sumAvg[i] + 1.e-7)
    end

    if bc isa Symbol
        update_boundary!(KS, ctr, face, dt, residual; coll = coll, bc = [bc, bc])
    else
        update_boundary!(KS, ctr, face, dt, residual; coll = coll, bc = bc)
    end

    return nothing

end

function update!(
    KS::X,
    ctr::Y,
    face::Z,
    dt,
    residual; # 1D / 2D
    coll = symbolize(KS.set.collision),
    bc = symbolize(KS.set.boundary),
    isMHD = true::Bool,
) where {
    X<:AbstractSolverSet,
    Y<:AA{ControlVolume1D3F,1},
    Z<:AA{Interface1D3F,1},
}

    sumRes = zero(ctr[1].w)
    sumAvg = zero(ctr[1].w)

    @inbounds Threads.@threads for i = 2:KS.pSpace.nx-1
        step!(KS, face[i], ctr[i], face[i+1], KS.ps.dx[i], dt, sumRes, sumAvg, coll, isMHD)
    end

    for i in eachindex(residual)
        residual[i] = sqrt(sumRes[i] * KS.pSpace.nx) / (sumAvg[i] + 1.e-7)
    end

    if bc isa Symbol
        update_boundary!(
            KS,
            ctr,
            face,
            dt,
            residual;
            coll = coll,
            bc = [bc, bc],
            isMHD = isMHD,
        )
    else
        update_boundary!(KS, ctr, face, dt, residual; coll = coll, bc = bc, isMHD = isMHD)
    end

    #=
    ng = 1 - first(eachindex(KS.pSpace.x))
    if bc == :extra
        for i in 1:ng
            ctr[1-i].w .= ctr[1].w
            ctr[1-i].prim .= ctr[1].prim
            ctr[1-i].h0 .= ctr[1].h0
            ctr[1-i].h1 .= ctr[1].h1
            ctr[1-i].h2 .= ctr[1].h2
            ctr[1-i].E .= ctr[1].E
            ctr[1-i].B .= ctr[1].B
            ctr[1-i].ϕ = deepcopy(ctr[1].ϕ)
            ctr[1-i].ψ = deepcopy(ctr[1].ψ)
            ctr[1-i].lorenz .= ctr[1].lorenz

            ctr[KS.pSpace.nx+i].w .= ctr[KS.pSpace.nx].w
            ctr[KS.pSpace.nx+i].prim .= ctr[KS.pSpace.nx].prim
            ctr[KS.pSpace.nx+i].h0 .= ctr[KS.pSpace.nx].h0
            ctr[KS.pSpace.nx+i].h1 .= ctr[KS.pSpace.nx].h1
            ctr[KS.pSpace.nx+i].h2 .= ctr[KS.pSpace.nx].h2
            ctr[KS.pSpace.nx+i].E .= ctr[KS.pSpace.nx].E
            ctr[KS.pSpace.nx+i].B .= ctr[KS.pSpace.nx].B
            ctr[KS.pSpace.nx+i].ϕ = deepcopy(ctr[KS.pSpace.nx].ϕ)
            ctr[KS.pSpace.nx+i].ψ = deepcopy(ctr[KS.pSpace.nx].ψ)
            ctr[KS.pSpace.nx+i].lorenz .= ctr[KS.pSpace.nx].lorenz
        end
    elseif bc == :period
        for i in 1:ng
            ctr[1-i].w .= ctr[KS.pSpace.nx+1-i].w
            ctr[1-i].prim .= ctr[KS.pSpace.nx+1-i].prim
            ctr[1-i].h0 .= ctr[KS.pSpace.nx+1-i].h0
            ctr[1-i].h1 .= ctr[KS.pSpace.nx+1-i].h1
            ctr[1-i].h2 .= ctr[KS.pSpace.nx+1-i].h2
            ctr[1-i].E .= ctr[KS.pSpace.nx+1-i].E
            ctr[1-i].B .= ctr[KS.pSpace.nx+1-i].B
            ctr[1-i].ϕ = deepcopy(ctr[KS.pSpace.nx+1-i].ϕ)
            ctr[1-i].ψ = deepcopy(ctr[KS.pSpace.nx+1-i].ψ)
            ctr[1-i].lorenz .= ctr[KS.pSpace.nx+1-i].lorenz

            ctr[KS.pSpace.nx+i].w .= ctr[i].w
            ctr[KS.pSpace.nx+i].prim .= ctr[i].prim
            ctr[KS.pSpace.nx+i].h0 .= ctr[i].h0
            ctr[KS.pSpace.nx+i].h1 .= ctr[i].h1
            ctr[KS.pSpace.nx+i].h2 .= ctr[i].h2
            ctr[KS.pSpace.nx+i].E .= ctr[i].E
            ctr[KS.pSpace.nx+i].B .= ctr[i].B
            ctr[KS.pSpace.nx+i].ϕ = deepcopy(ctr[i].ϕ)
            ctr[KS.pSpace.nx+i].ψ = deepcopy(ctr[i].ψ)
            ctr[KS.pSpace.nx+i].lorenz .= ctr[i].lorenz
        end
    elseif bc == :balance
        @. ctr[0].w = 0.5 * (ctr[-1].w + ctr[1].w)
        @. ctr[0].prim = 0.5 * (ctr[-1].prim + ctr[1].prim)
        @. ctr[0].h0 = 0.5 * (ctr[-1].h0 + ctr[1].h0)
        @. ctr[0].h1 = 0.5 * (ctr[-1].h1 + ctr[1].h1)
        @. ctr[0].h2 = 0.5 * (ctr[-1].h2 + ctr[1].h2)
        @. ctr[0].E = 0.5 * (ctr[-1].E + ctr[1].E)
        @. ctr[0].B = 0.5 * (ctr[-1].B + ctr[1].B)
        ctr[0].ϕ = 0.5 * (ctr[-1].ϕ + ctr[1].ϕ)
        ctr[0].ψ = 0.5 * (ctr[-1].ψ + ctr[1].ψ)
        @. ctr[0].lorenz = 0.5 * (ctr[-1].lorenz + ctr[1].lorenz)

        @. ctr[KS.pSpace.nx+1].w = 0.5 * (ctr[KS.pSpace.nx].w + ctr[KS.pSpace.nx+2].w)
        @. ctr[KS.pSpace.nx+1].prim = 0.5 * (ctr[KS.pSpace.nx].prim + ctr[KS.pSpace.nx+2].prim)
        @. ctr[KS.pSpace.nx+1].h0 = 0.5 * (ctr[KS.pSpace.nx].h0 + ctr[KS.pSpace.nx+2].h0)
        @. ctr[KS.pSpace.nx+1].h1 = 0.5 * (ctr[KS.pSpace.nx].h1 + ctr[KS.pSpace.nx+2].h1)
        @. ctr[KS.pSpace.nx+1].h2 = 0.5 * (ctr[KS.pSpace.nx].h2 + ctr[KS.pSpace.nx+2].h2)
        @. ctr[KS.pSpace.nx+1].E = 0.5 * (ctr[KS.pSpace.nx].E + ctr[KS.pSpace.nx+2].E)
        @. ctr[KS.pSpace.nx+1].B = 0.5 * (ctr[KS.pSpace.nx].B + ctr[KS.pSpace.nx+2].B)
        ctr[KS.pSpace.nx+1].ϕ = 0.5 * (ctr[KS.pSpace.nx].ϕ + ctr[KS.pSpace.nx+2].ϕ)
        ctr[KS.pSpace.nx+1].ψ = 0.5 * (ctr[KS.pSpace.nx].ψ + ctr[KS.pSpace.nx+2].ψ)
        @. ctr[KS.pSpace.nx+1].lorenz = 0.5 * (ctr[KS.pSpace.nx].lorenz + ctr[KS.pSpace.nx+2].lorenz)
    else
    end
    =#
    return nothing

end

function update!(
    KS::X,
    ctr::Y,
    face::Z,
    dt,
    residual; # 1D / 2D
    coll = symbolize(KS.set.collision),
    bc = symbolize(KS.set.boundary),
    isMHD = true::Bool,
) where {
    X<:AbstractSolverSet,
    Y<:AA{ControlVolume1D4F,1},
    Z<:AA{Interface1D4F,1},
}

    sumRes = zero(ctr[1].w)
    sumAvg = zero(ctr[1].w)

    @inbounds Threads.@threads for i = 2:KS.pSpace.nx-1
        step!(KS, face[i], ctr[i], face[i+1], KS.ps.dx[i], dt, sumRes, sumAvg, coll, isMHD)
    end

    for i in eachindex(residual)
        residual[i] = sqrt(sumRes[i] * KS.pSpace.nx) / (sumAvg[i] + 1.e-7)
    end

    if bc isa Symbol
        update_boundary!(
            KS,
            ctr,
            face,
            dt,
            residual;
            coll = coll,
            bc = [bc, bc],
            isMHD = isMHD,
        )
    else
        update_boundary!(KS, ctr, face, dt, residual; coll = coll, bc = bc, isMHD = isMHD)
    end

    return nothing

end

function update!(
    KS::X,
    ctr::Y,
    a1face::Z,
    a2face::Z,
    dt,
    residual;
    coll = symbolize(KS.set.collision),
    bc = symbolize(KS.set.boundary),
) where {
    X<:AbstractSolverSet,
    Y<:AA{ControlVolume2D,2},
    Z<:AA{Interface2D,2},
}

    nx, ny, dx, dy = begin
        if KS.ps isa CSpace2D
            KS.ps.nr, KS.ps.nθ, KS.ps.dr, KS.ps.darc
        else
            KS.ps.nx, KS.ps.ny, KS.ps.dx, KS.ps.dy
        end
    end

    sumRes = zero(ctr[1].w)
    sumAvg = zero(ctr[1].w)

    if ndims(sumRes) == 1

        @inbounds for j = 2:ny-1
            for i = 2:nx-1
                step!(
                    ctr[i, j].w,
                    ctr[i, j].prim,
                    a1face[i, j].fw,
                    a1face[i+1, j].fw,
                    a2face[i, j].fw,
                    a2face[i, j+1].fw,
                    KS.gas.γ,
                    dx[i, j] * dy[i, j],
                    sumRes,
                    sumAvg,
                    coll,
                )
            end
        end

    end

    for i in eachindex(residual)
        residual[i] = sqrt(sumRes[i] * nx * ny) / (sumAvg[i] + 1.e-7)
    end

    if bc isa Symbol
        update_boundary!(
            KS,
            ctr,
            a1face,
            a2face,
            dt,
            residual;
            coll = coll,
            bc = [bc, bc, bc, bc],
        )
    else
        update_boundary!(KS, ctr, a1face, a2face, dt, residual; coll = coll, bc = bc)
    end

    return nothing

end

function update!(
    KS::X,
    ctr::Y,
    a1face::Z,
    a2face::Z,
    dt,
    residual;
    coll = symbolize(KS.set.collision),
    bc = symbolize(KS.set.boundary),
) where {
    X<:AbstractSolverSet,
    Y<:AA{ControlVolume2D1F,2},
    Z<:AA{Interface2D1F,2},
}

    nx, ny, dx, dy = begin
        if KS.ps isa CSpace2D
            KS.ps.nr, KS.ps.nθ, KS.ps.dr, KS.ps.darc
        else
            KS.ps.nx, KS.ps.ny, KS.ps.dx, KS.ps.dy
        end
    end

    sumRes = zero(ctr[1].w)
    sumAvg = zero(ctr[1].w)

    if ndims(sumRes) == 1

        @inbounds for j = 2:ny-1
            for i = 2:nx-1
                step!(
                    ctr[i, j].w,
                    ctr[i, j].prim,
                    ctr[i, j].f,
                    a1face[i, j].fw,
                    a1face[i, j].ff,
                    a1face[i+1, j].fw,
                    a1face[i+1, j].ff,
                    a2face[i, j].fw,
                    a2face[i, j].ff,
                    a2face[i, j+1].fw,
                    a2face[i, j+1].ff,
                    KS.vSpace.u,
                    KS.vSpace.v,
                    KS.vSpace.weights,
                    KS.gas.γ,
                    KS.gas.μᵣ,
                    KS.gas.ω,
                    KS.gas.Pr,
                    dx[i, j] * dy[i, j],
                    dt,
                    sumRes,
                    sumAvg,
                    coll,
                )
            end
        end

    end

    for i in eachindex(residual)
        residual[i] = sqrt(sumRes[i] * nx * ny) / (sumAvg[i] + 1.e-7)
    end

    if bc isa Symbol
        update_boundary!(
            KS,
            ctr,
            a1face,
            a2face,
            dt,
            residual;
            coll = coll,
            bc = [bc, bc, bc, bc],
        )
    else
        update_boundary!(KS, ctr, a1face, a2face, dt, residual; coll = coll, bc = bc)
    end

    return nothing

end

#--- 2d2f ---#
function update!(
    KS::X,
    ctr::Y,
    a1face::Z,
    a2face::Z,
    dt,
    residual;
    coll = symbolize(KS.set.collision),
    bc = symbolize(KS.set.boundary),
) where {
    X<:AbstractSolverSet,
    Y<:AA{ControlVolume2D2F,2},
    Z<:AA{Interface2D2F,2},
}

    nx, ny, dx, dy = begin
        if KS.ps isa CSpace2D
            KS.ps.nr, KS.ps.nθ, KS.ps.dr, KS.ps.darc
        else
            KS.ps.nx, KS.ps.ny, KS.ps.dx, KS.ps.dy
        end
    end

    sumRes = zero(ctr[1].w)
    sumAvg = zero(ctr[1].w)

    if ndims(sumRes) == 1

        @inbounds for j = 2:ny-1
            for i = 2:nx-1
                step!(
                    ctr[i, j].w,
                    ctr[i, j].prim,
                    ctr[i, j].h,
                    ctr[i, j].b,
                    a1face[i, j].fw,
                    a1face[i, j].fh,
                    a1face[i, j].fb,
                    a1face[i+1, j].fw,
                    a1face[i+1, j].fh,
                    a1face[i+1, j].fb,
                    a2face[i, j].fw,
                    a2face[i, j].fh,
                    a2face[i, j].fb,
                    a2face[i, j+1].fw,
                    a2face[i, j+1].fh,
                    a2face[i, j+1].fb,
                    KS.vSpace.u,
                    KS.vSpace.v,
                    KS.vSpace.weights,
                    KS.gas.K,
                    KS.gas.γ,
                    KS.gas.μᵣ,
                    KS.gas.ω,
                    KS.gas.Pr,
                    dx[i, j] * dy[i, j],
                    dt,
                    sumRes,
                    sumAvg,
                    coll,
                )
            end
        end

    end

    for i in eachindex(residual)
        residual[i] = sqrt(sumRes[i] * nx * ny) / (sumAvg[i] + 1.e-7)
    end

    if bc isa Symbol
        update_boundary!(
            KS,
            ctr,
            a1face,
            a2face,
            dt,
            residual;
            coll = coll,
            bc = [bc, bc, bc, bc],
        )
    else
        update_boundary!(KS, ctr, a1face, a2face, dt, residual; coll = coll, bc = bc)
    end

    return nothing

end

function update!(
    KS::X,
    ctr::Y,
    face::Z,
    dt,
    residual; # 1D / 2D
    coll = symbolize(KS.set.collision),
    bc = symbolize(KS.set.boundary),
) where {
    X<:AbstractSolverSet,
    Y<:AA{ControlVolumeUS,1},
    Z<:AA{Interface2D,1},
}

    sumRes = zero(ctr[1].w)
    sumAvg = zero(ctr[1].w)

    @inbounds Threads.@threads for i in eachindex(ctr)
        if KS.ps.cellType[i] in (0, 2)
            dirc = [sign(dot(ctr[i].n[j], face[KS.ps.cellFaces[i, j]].n)) for j = 1:3]

            step!(
                ctr[i].w,
                ctr[i].prim,
                face[KS.ps.cellFaces[i, 1]].fw,
                face[KS.ps.cellFaces[i, 2]].fw,
                face[KS.ps.cellFaces[i, 3]].fw,
                KS.gas.γ,
                KS.pSpace.cellArea[i],
                dirc,
                sumRes,
                sumAvg,
            )
        end
    end

    for i in eachindex(residual)
        residual[i] = sqrt(sumRes[i] * size(KS.ps.cellid, 1)) / (sumAvg[i] + 1.e-7)
    end

    update_boundary!(KS, ctr, face, dt, residual; coll = coll, bc = bc)

    return nothing

end

function update!(
    KS::X,
    ctr::Y,
    face::Z,
    dt,
    residual; # 1D / 2D
    coll = symbolize(KS.set.collision),
    bc = symbolize(KS.set.boundary),
) where {
    X<:AbstractSolverSet,
    Y<:AA{ControlVolumeUS1F,1},
    Z<:AA{Interface2D1F,1},
}

    sumRes = zero(ctr[1].w)
    sumAvg = zero(ctr[1].w)

    @inbounds Threads.@threads for i in eachindex(ctr)
        if KS.ps.cellType[i] in (0, 2)
            dirc = [sign(dot(ctr[i].n[j], face[KS.ps.cellFaces[i, j]].n)) for j = 1:3]

            step!(
                ctr[i].w,
                ctr[i].prim,
                ctr[i].f,
                face[KS.ps.cellFaces[i, 1]].fw,
                face[KS.ps.cellFaces[i, 1]].ff,
                face[KS.ps.cellFaces[i, 2]].fw,
                face[KS.ps.cellFaces[i, 2]].ff,
                face[KS.ps.cellFaces[i, 3]].fw,
                face[KS.ps.cellFaces[i, 3]].ff,
                KS.vSpace.u,
                KS.vSpace.v,
                KS.vSpace.weights,
                KS.gas.K,
                KS.gas.γ,
                KS.gas.μᵣ,
                KS.gas.ω,
                KS.gas.Pr,
                KS.pSpace.cellArea[i],
                dirc,
                dt,
                sumRes,
                sumAvg,
                coll,
            )
        end
    end

    for i in eachindex(residual)
        residual[i] = sqrt(sumRes[i] * size(KS.ps.cellid, 1)) / (sumAvg[i] + 1.e-7)
    end

    update_boundary!(KS, ctr, face, dt, residual; coll = coll, bc = bc)

    return nothing

end

function update!(
    KS::X,
    ctr::Y,
    face::Z,
    dt,
    residual; # 1D / 2D
    coll = symbolize(KS.set.collision),
    bc = symbolize(KS.set.boundary),
) where {
    X<:AbstractSolverSet,
    Y<:AA{ControlVolumeUS2F,1},
    Z<:AA{Interface2D2F,1},
}

    sumRes = zero(ctr[1].w)
    sumAvg = zero(ctr[1].w)

    @inbounds Threads.@threads for i in eachindex(ctr)
        if KS.ps.cellType[i] in (0, 2)
            dirc = [sign(dot(ctr[i].n[j], face[KS.ps.cellFaces[i, j]].n)) for j = 1:3]

            step!(
                ctr[i].w,
                ctr[i].prim,
                ctr[i].h,
                ctr[i].b,
                face[KS.ps.cellFaces[i, 1]].fw,
                face[KS.ps.cellFaces[i, 1]].fh,
                face[KS.ps.cellFaces[i, 1]].fb,
                face[KS.ps.cellFaces[i, 2]].fw,
                face[KS.ps.cellFaces[i, 2]].fh,
                face[KS.ps.cellFaces[i, 2]].fb,
                face[KS.ps.cellFaces[i, 3]].fw,
                face[KS.ps.cellFaces[i, 3]].fh,
                face[KS.ps.cellFaces[i, 3]].fb,
                KS.vSpace.u,
                KS.vSpace.v,
                KS.vSpace.weights,
                KS.gas.K,
                KS.gas.γ,
                KS.gas.μᵣ,
                KS.gas.ω,
                KS.gas.Pr,
                KS.pSpace.cellArea[i],
                dirc,
                dt,
                sumRes,
                sumAvg,
                coll,
            )
        end
    end

    for i in eachindex(residual)
        residual[i] = sqrt(sumRes[i] * size(KS.ps.cellid, 1)) / (sumAvg[i] + 1.e-7)
    end

    update_boundary!(KS, ctr, face, dt, residual; coll = coll, bc = bc)

    return nothing

end
