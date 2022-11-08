"""
$(SIGNATURES)

Update algorithm

## Arguments
``KS``: SolverSet
``ctr``: array of cell-centered solution
``face``: array of cell interface
``dt``: time step
``residual``: residual
``coll``: collision operator
``bc``: boundary condition
``fn``: update function
``st``: step function
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
    st = fn,
) where {TC<:Union{ControlVolume,ControlVolume1D},TF<:Union{Interface,Interface1D}}

    sumRes, sumAvg = 0.0, 0.0

    @inbounds @threads for i = 1:KS.ps.nx
        ctr[i].w, sumRes, sumAvg = fn(
            KS,
            ctr[i],
            face[i],
            face[i+1],
            (dt, KS.ps.dx[i], sumRes, sumAvg),
            coll;
            st = st,
        )
    end

    residual = sqrt(sumRes * KS.ps.nx) / (sumAvg + 1.e-7) |> Float64

    bcs = ifelse(bc isa Symbol, [bc, bc], bc)
    if bcs[1] == :period
        ctr[0].w = ctr[KS.ps.nx].w
        ctr[0].sw = ctr[KS.ps.nx].sw
        ctr[KS.ps.nx+1].w = ctr[1].w
        ctr[KS.ps.nx+1].sw = ctr[1].sw
    elseif bcs[1] == :extra
        ctr[0].w = ctr[1].w
    end

    if bcs[2] == :extra
        ctr[KS.ps.nx+1].w = ctr[KS.ps.nx].w
    end

    return residual

end

function update!(
    KS::AbstractSolverSet,
    ctr::AV{TC},
    face,
    dt,
    residual::AA; # 1D / 2D
    coll = symbolize(KS.set.collision),
    bc = symbolize(KS.set.boundary),
    fn = step!,
    st = fn,
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

    sumRes = zero(ctr[1].w)
    sumAvg = zero(ctr[1].w)

    @inbounds @threads for i = 2:KS.ps.nx-1
        fn(KS, ctr[i], face[i], face[i+1], (dt, KS.ps.dx[i], sumRes, sumAvg), coll; st = st)
    end

    for i in eachindex(residual)
        residual[i] = sqrt(sumRes[i] * KS.ps.nx) / (sumAvg[i] + 1.e-7)
    end

    update_boundary!(KS, ctr, face, dt, residual; bc = bc, fn = fn, st = st)

    return nothing

end

"""
$(SIGNATURES)
"""
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
    st = fn,
)

    sumRes = zero(ctr[1].w)
    sumAvg = zero(ctr[1].w)

    @inbounds @threads for i = 2:KS.ps.nx-1
        fn(KS, ctr[i], face[i], face[i+1], KS.ps.dx[i], dt, sumRes, sumAvg, coll, isMHD)
    end

    for i in eachindex(residual)
        residual[i] = sqrt(sumRes[i] * KS.ps.nx) / (sumAvg[i] + 1.e-7)
    end

    update_boundary!(
        KS,
        ctr,
        face,
        dt,
        residual;
        coll = coll,
        bc = bc,
        isMHD = isMHD,
        fn = fn,
        st = st,
    )

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
    st = fn,
)

    sumRes = zero(ctr[1].w)
    sumAvg = zero(ctr[1].w)

    @inbounds @threads for i = 2:KS.ps.nx-1
        fn(KS, ctr[i], face[i], face[i+1], KS.ps.dx[i], dt, sumRes, sumAvg, coll, isMHD)
    end

    for i in eachindex(residual)
        residual[i] = sqrt(sumRes[i] * KS.ps.nx) / (sumAvg[i] + 1.e-7)
    end

    update_boundary!(
        KS,
        ctr,
        face,
        dt,
        residual;
        coll = coll,
        bc = bc,
        isMHD = isMHD,
        fn = fn,
        st = st,
    )

    return nothing

end

"""
$(SIGNATURES)

Update algorithm

## Arguments
``KS``: SolverSet
``ctr``: array of cell-centered solution
``a1face``: array of cell interface perpendicular to `x` axis
``a2face``: array of cell interface perpendicular to `y` axis
``dt``: time step
``residual``: residual
``coll``: collision operator
``bc``: boundary condition
``fn``: update function
``st``: step function
"""
function update!(
    KS::AbstractSolverSet,
    ctr::AM{T},
    a1face,
    a2face,
    dt,
    residual;
    coll = symbolize(KS.set.collision),
    bc = symbolize(KS.set.boundary),
    fn = step!,
    st = fn,
) where {
    T<:Union{
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

    sumRes = zero(ctr[1].w)
    sumAvg = zero(ctr[1].w)

    @inbounds @threads for j = 2:ny-1
        for i = 2:nx-1
            fn(
                KS,
                ctr[i, j],
                a1face[i, j],
                a1face[i+1, j],
                a2face[i, j],
                a2face[i, j+1],
                (dt, dx[i, j] * dy[i, j], sumRes, sumAvg),
                coll;
                st = st,
            )
        end
    end

    for i in eachindex(residual)
        residual[i] = sqrt(sumRes[i] * nx * ny) / (sumAvg[i] + 1.e-7)
    end

    update_boundary!(
        KS,
        ctr,
        a1face,
        a2face,
        dt,
        residual;
        coll = coll,
        bc = bc,
        fn = fn,
        st = st,
    )

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
    st = fn,
) where {T<:Union{Interface,Interface2D}}

    sumRes = zero(ctr[1].w)
    sumAvg = zero(ctr[1].w)

    @inbounds @threads for i in eachindex(ctr)
        if KS.ps.cellType[i] in (0, 2)
            dirc = [sign(dot(ctr[i].n[j], face[KS.ps.cellFaces[i, j]].n)) for j = 1:3]

            st(
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

    update_boundary!(KS, ctr, face, dt, residual; coll = coll, bc = bc, fn = fn, st = st)

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
    st = fn,
) where {T<:Union{Interface1F,Interface2D1F}}

    sumRes = zero(ctr[1].w)
    sumAvg = zero(ctr[1].w)

    @inbounds @threads for i in eachindex(ctr)
        if KS.ps.cellType[i] in (0, 2)
            dirc = [sign(dot(ctr[i].n[j], face[KS.ps.cellFaces[i, j]].n)) for j = 1:3]

            st(
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

    update_boundary!(KS, ctr, face, dt, residual; coll = coll, bc = bc, fn = fn, st = st)

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
    st = fn,
) where {T<:Union{Interface2F,Interface2D2F}}

    sumRes = zero(ctr[1].w)
    sumAvg = zero(ctr[1].w)

    @inbounds @threads for i in eachindex(ctr)
        if KS.ps.cellType[i] in (0, 2)
            dirc = [sign(dot(ctr[i].n[j], face[KS.ps.cellFaces[i, j]].n)) for j = 1:3]

            st(
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

    update_boundary!(KS, ctr, face, dt, residual; coll = coll, bc = bc, fn = fn, st = st)

    return nothing

end
