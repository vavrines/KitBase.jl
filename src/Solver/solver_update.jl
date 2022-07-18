"""
$(SIGNATURES)

Update flow variables

1D scalar
"""
function update!(
    KS,
    ctr::AV{TC},
    face::AV{TF},
    dt,
    residual::Number;
    coll = symbolize(KS.set.collision),
    bc = symbolize(KS.set.boundary),
    fn = step!,
) where {TC<:Union{ControlVolume,ControlVolume1D},TF<:Union{Interface,Interface1D}}

    sumRes = zero(ctr[1].w)
    sumAvg = zero(sumRes)

    @inbounds @threads for i = 1:KS.ps.nx
        ctr[i].w, sumRes, sumAvg = fn(
            KS,
            ctr[i],
            face[i],
            face[i+1],
            (dt, KS.ps.dx[i], sumRes, sumAvg),
            coll,
        )
    end

    residual = sqrt(sumRes * KS.ps.nx) / (sumAvg + 1.e-7)

    if bc[1] == :period
        ctr[0].w = ctr[KS.ps.nx].w
        ctr[0].sw = ctr[KS.ps.nx].sw
        ctr[KS.ps.nx+1].w = ctr[1].w
        ctr[KS.ps.nx+1].sw = ctr[1].sw
    elseif bc[1] == :extra
        ctr[0].w = ctr[1].w
    end

    if bc[2] == :extra
        ctr[KS.ps.nx+1].w = ctr[KS.ps.nx].w
    end

    return residual

end

"""
$(SIGNATURES)

- 1D0F pure & mixture
- 1D1F
- 1D2F pure & mixture
"""
function update!(
    KS::AbstractSolverSet,
    ctr::AV{TC},
    face,
    dt,
    residual::AA; # 1D / 2D
    coll = symbolize(KS.set.collision),
    bc = symbolize(KS.set.boundary),
    fn = step!,
) where {TC<:Union{ControlVolume,ControlVolume1D,ControlVolume1F,ControlVolume1D1F,ControlVolume2F,ControlVolume1D2F}}

    sumRes = zero(ctr[1].w)
    sumAvg = zero(ctr[1].w)

    @inbounds @threads for i = 2:KS.ps.nx-1
        fn(
            KS,
            ctr[i],
            face[i],
            face[i+1],
            (dt, KS.ps.dx[i], sumRes, sumAvg),
            coll,
        )
    end

    for i in eachindex(residual)
        residual[i] = sqrt(sumRes[i] * KS.ps.nx) / (sumAvg[i] + 1.e-7)
    end

    if bc isa Symbol
        update_boundary!(KS, ctr, face, dt, residual; bc = [bc, bc], fn = fn)
    else
        update_boundary!(KS, ctr, face, dt, residual; bc = bc, fn = fn)
    end

    return nothing

end

function update!(
    KS::AbstractSolverSet,
    ctr::AV{ControlVolume1D3F},
    face::AV{Interface1D3F},
    dt,
    residual; # 1D / 2D
    coll = symbolize(KS.set.collision),
    bc = symbolize(KS.set.boundary),
    isMHD = true::Bool,
    fn = step!,
)

    sumRes = zero(ctr[1].w)
    sumAvg = zero(ctr[1].w)

    @inbounds @threads for i = 2:KS.ps.nx-1
        fn(KS, ctr[i], face[i], face[i+1], KS.ps.dx[i], dt, sumRes, sumAvg, coll, isMHD)
    end

    for i in eachindex(residual)
        residual[i] = sqrt(sumRes[i] * KS.ps.nx) / (sumAvg[i] + 1.e-7)
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
            fn = fn,
        )
    else
        update_boundary!(KS, ctr, face, dt, residual; coll = coll, bc = bc, isMHD = isMHD, fn = fn)
    end

    #=
    ng = 1 - first(eachindex(KS.ps.x))
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

            ctr[KS.ps.nx+i].w .= ctr[KS.ps.nx].w
            ctr[KS.ps.nx+i].prim .= ctr[KS.ps.nx].prim
            ctr[KS.ps.nx+i].h0 .= ctr[KS.ps.nx].h0
            ctr[KS.ps.nx+i].h1 .= ctr[KS.ps.nx].h1
            ctr[KS.ps.nx+i].h2 .= ctr[KS.ps.nx].h2
            ctr[KS.ps.nx+i].E .= ctr[KS.ps.nx].E
            ctr[KS.ps.nx+i].B .= ctr[KS.ps.nx].B
            ctr[KS.ps.nx+i].ϕ = deepcopy(ctr[KS.ps.nx].ϕ)
            ctr[KS.ps.nx+i].ψ = deepcopy(ctr[KS.ps.nx].ψ)
            ctr[KS.ps.nx+i].lorenz .= ctr[KS.ps.nx].lorenz
        end
    elseif bc == :period
        for i in 1:ng
            ctr[1-i].w .= ctr[KS.ps.nx+1-i].w
            ctr[1-i].prim .= ctr[KS.ps.nx+1-i].prim
            ctr[1-i].h0 .= ctr[KS.ps.nx+1-i].h0
            ctr[1-i].h1 .= ctr[KS.ps.nx+1-i].h1
            ctr[1-i].h2 .= ctr[KS.ps.nx+1-i].h2
            ctr[1-i].E .= ctr[KS.ps.nx+1-i].E
            ctr[1-i].B .= ctr[KS.ps.nx+1-i].B
            ctr[1-i].ϕ = deepcopy(ctr[KS.ps.nx+1-i].ϕ)
            ctr[1-i].ψ = deepcopy(ctr[KS.ps.nx+1-i].ψ)
            ctr[1-i].lorenz .= ctr[KS.ps.nx+1-i].lorenz

            ctr[KS.ps.nx+i].w .= ctr[i].w
            ctr[KS.ps.nx+i].prim .= ctr[i].prim
            ctr[KS.ps.nx+i].h0 .= ctr[i].h0
            ctr[KS.ps.nx+i].h1 .= ctr[i].h1
            ctr[KS.ps.nx+i].h2 .= ctr[i].h2
            ctr[KS.ps.nx+i].E .= ctr[i].E
            ctr[KS.ps.nx+i].B .= ctr[i].B
            ctr[KS.ps.nx+i].ϕ = deepcopy(ctr[i].ϕ)
            ctr[KS.ps.nx+i].ψ = deepcopy(ctr[i].ψ)
            ctr[KS.ps.nx+i].lorenz .= ctr[i].lorenz
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

        @. ctr[KS.ps.nx+1].w = 0.5 * (ctr[KS.ps.nx].w + ctr[KS.ps.nx+2].w)
        @. ctr[KS.ps.nx+1].prim = 0.5 * (ctr[KS.ps.nx].prim + ctr[KS.ps.nx+2].prim)
        @. ctr[KS.ps.nx+1].h0 = 0.5 * (ctr[KS.ps.nx].h0 + ctr[KS.ps.nx+2].h0)
        @. ctr[KS.ps.nx+1].h1 = 0.5 * (ctr[KS.ps.nx].h1 + ctr[KS.ps.nx+2].h1)
        @. ctr[KS.ps.nx+1].h2 = 0.5 * (ctr[KS.ps.nx].h2 + ctr[KS.ps.nx+2].h2)
        @. ctr[KS.ps.nx+1].E = 0.5 * (ctr[KS.ps.nx].E + ctr[KS.ps.nx+2].E)
        @. ctr[KS.ps.nx+1].B = 0.5 * (ctr[KS.ps.nx].B + ctr[KS.ps.nx+2].B)
        ctr[KS.ps.nx+1].ϕ = 0.5 * (ctr[KS.ps.nx].ϕ + ctr[KS.ps.nx+2].ϕ)
        ctr[KS.ps.nx+1].ψ = 0.5 * (ctr[KS.ps.nx].ψ + ctr[KS.ps.nx+2].ψ)
        @. ctr[KS.ps.nx+1].lorenz = 0.5 * (ctr[KS.ps.nx].lorenz + ctr[KS.ps.nx+2].lorenz)
    else
    end
    =#
    return nothing

end

function update!(
    KS::AbstractSolverSet,
    ctr::AV{ControlVolume1D4F},
    face::AV{Interface1D4F},
    dt,
    residual; # 1D / 2D
    coll = symbolize(KS.set.collision),
    bc = symbolize(KS.set.boundary),
    isMHD = true::Bool,
    fn = step!,
)

    sumRes = zero(ctr[1].w)
    sumAvg = zero(ctr[1].w)

    @inbounds @threads for i = 2:KS.ps.nx-1
        fn(KS, ctr[i], face[i], face[i+1], KS.ps.dx[i], dt, sumRes, sumAvg, coll, isMHD)
    end

    for i in eachindex(residual)
        residual[i] = sqrt(sumRes[i] * KS.ps.nx) / (sumAvg[i] + 1.e-7)
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
            fn = fn,
        )
    else
        update_boundary!(KS, ctr, face, dt, residual; coll = coll, bc = bc, isMHD = isMHD, fn = fn)
    end

    return nothing

end

function update!(
    KS::AbstractSolverSet,
    ctr::AM{TC},
    a1face::AM{TF},
    a2face::AM{TF},
    dt,
    residual;
    coll = symbolize(KS.set.collision),
    bc = symbolize(KS.set.boundary),
    fn = step!,
) where {TC<:Union{ControlVolume,ControlVolume2D},TF<:Union{Interface,Interface2D}}

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

        @inbounds @threads for j = 2:ny-1
            for i = 2:nx-1
                fn(
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
            fn = fn,
        )
    else
        update_boundary!(KS, ctr, a1face, a2face, dt, residual; coll = coll, bc = bc, fn = fn)
    end

    return nothing

end

function update!(
    KS::AbstractSolverSet,
    ctr::AM{TC},
    a1face::AM{TF},
    a2face::AM{TF},
    dt,
    residual;
    coll = symbolize(KS.set.collision),
    bc = symbolize(KS.set.boundary),
    fn = step!,
) where {TC<:Union{ControlVolume1F,ControlVolume2D1F},TF<:Union{Interface1F,Interface2D1F}}

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

        @inbounds @threads for j = 2:ny-1
            for i = 2:nx-1
                fn(
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
                    KS.vs.u,
                    KS.vs.v,
                    KS.vs.weights,
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
            fn = fn,
        )
    else
        update_boundary!(KS, ctr, a1face, a2face, dt, residual; coll = coll, bc = bc, fn = fn)
    end

    return nothing

end

#--- 2d2f ---#
function update!(
    KS::AbstractSolverSet,
    ctr::AM{TC},
    a1face::AM{TF},
    a2face::AM{TF},
    dt,
    residual;
    coll = symbolize(KS.set.collision),
    bc = symbolize(KS.set.boundary),
    fn = step!,
) where {TC<:Union{ControlVolume2F,ControlVolume2D2F},TF<:Union{Interface2F,Interface2D2F}}

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

        @inbounds @threads for j = 2:ny-1
            for i = 2:nx-1
                fn(
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
                    KS.vs.u,
                    KS.vs.v,
                    KS.vs.weights,
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
            fn = fn,
        )
    else
        update_boundary!(KS, ctr, a1face, a2face, dt, residual; coll = coll, bc = bc, fn = fn)
    end

    return nothing

end

function update!(
    KS::AbstractSolverSet,
    ctr::AV{ControlVolumeUS},
    face::AV{T},
    dt,
    residual; # 1D / 2D
    coll = symbolize(KS.set.collision),
    bc = symbolize(KS.set.boundary),
    fn = step!,
) where {T<:Union{Interface,Interface2D}}

    sumRes = zero(ctr[1].w)
    sumAvg = zero(ctr[1].w)

    @inbounds @threads for i in eachindex(ctr)
        if KS.ps.cellType[i] in (0, 2)
            dirc = [sign(dot(ctr[i].n[j], face[KS.ps.cellFaces[i, j]].n)) for j = 1:3]

            fn(
                ctr[i].w,
                ctr[i].prim,
                face[KS.ps.cellFaces[i, 1]].fw,
                face[KS.ps.cellFaces[i, 2]].fw,
                face[KS.ps.cellFaces[i, 3]].fw,
                KS.gas.γ,
                KS.ps.cellArea[i],
                dirc,
                sumRes,
                sumAvg,
            )
        end
    end

    for i in eachindex(residual)
        residual[i] = sqrt(sumRes[i] * size(KS.ps.cellid, 1)) / (sumAvg[i] + 1.e-7)
    end

    update_boundary!(KS, ctr, face, dt, residual; coll = coll, bc = bc, fn = fn)

    return nothing

end

function update!(
    KS::AbstractSolverSet,
    ctr::AV{ControlVolumeUS1F},
    face::AV{T},
    dt,
    residual; # 1D / 2D
    coll = symbolize(KS.set.collision),
    bc = symbolize(KS.set.boundary),
    fn = step!,
) where {T<:Union{Interface1F,Interface2D1F}}

    sumRes = zero(ctr[1].w)
    sumAvg = zero(ctr[1].w)

    @inbounds @threads for i in eachindex(ctr)
        if KS.ps.cellType[i] in (0, 2)
            dirc = [sign(dot(ctr[i].n[j], face[KS.ps.cellFaces[i, j]].n)) for j = 1:3]

            fn(
                ctr[i].w,
                ctr[i].prim,
                ctr[i].f,
                face[KS.ps.cellFaces[i, 1]].fw,
                face[KS.ps.cellFaces[i, 1]].ff,
                face[KS.ps.cellFaces[i, 2]].fw,
                face[KS.ps.cellFaces[i, 2]].ff,
                face[KS.ps.cellFaces[i, 3]].fw,
                face[KS.ps.cellFaces[i, 3]].ff,
                KS.vs.u,
                KS.vs.v,
                KS.vs.weights,
                KS.gas.K,
                KS.gas.γ,
                KS.gas.μᵣ,
                KS.gas.ω,
                KS.gas.Pr,
                KS.ps.cellArea[i],
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

    update_boundary!(KS, ctr, face, dt, residual; coll = coll, bc = bc, fn = fn)

    return nothing

end

function update!(
    KS::AbstractSolverSet,
    ctr::AV{ControlVolumeUS2F},
    face::AV{T},
    dt,
    residual; # 1D / 2D
    coll = symbolize(KS.set.collision),
    bc = symbolize(KS.set.boundary),
    fn = step!,
) where {T<:Union{Interface2F,Interface2D2F}}

    sumRes = zero(ctr[1].w)
    sumAvg = zero(ctr[1].w)

    @inbounds @threads for i in eachindex(ctr)
        if KS.ps.cellType[i] in (0, 2)
            dirc = [sign(dot(ctr[i].n[j], face[KS.ps.cellFaces[i, j]].n)) for j = 1:3]

            fn(
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
                KS.vs.u,
                KS.vs.v,
                KS.vs.weights,
                KS.gas.K,
                KS.gas.γ,
                KS.gas.μᵣ,
                KS.gas.ω,
                KS.gas.Pr,
                KS.ps.cellArea[i],
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

    update_boundary!(KS, ctr, face, dt, residual; coll = coll, bc = bc, fn = fn)

    return nothing

end
