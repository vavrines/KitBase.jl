
"""
Update flow variables

"""
function update!(
    KS::X,
    ctr::Y,
    face::Z,
    dt,
    residual; # 1D / 2D
    coll = :bgk::Symbol,
    bc = :fix::Symbol,
) where {
    X<:AbstractSolverSet,
    Y<:AbstractArray{ControlVolume1D1F,1},
    Z<:AbstractArray{Interface1D1F,1},
}

    sumRes = zeros(axes(KS.ib.wL))
    sumAvg = zeros(axes(KS.ib.wL))

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
                ctr[i].dx,
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
                ctr[i].dx,
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

    update_boundary!(KS, ctr, face, dt, residual; coll = coll, bc = bc, isMHD = false)

end

function update!(
    KS::X,
    ctr::Y,
    face::Z,
    dt,
    residual; # 1D / 2D
    coll = :bgk::Symbol,
    bc = :extra::Symbol,
) where {
    X<:AbstractSolverSet,
    Y<:AbstractArray{ControlVolume1D2F,1},
    Z<:AbstractArray{Interface1D2F,1},
}

    sumRes = zeros(axes(KS.ib.wL))
    sumAvg = zeros(axes(KS.ib.wL))

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
                ctr[i].dx,
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
                ctr[i].dx,
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

    update_boundary!(KS, ctr, face, dt, residual; coll = coll, bc = bc, isMHD = false)

end

function update!(
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
    Y<:AbstractArray{ControlVolume1D3F,1},
    Z<:AbstractArray{Interface1D3F,1},
}

    sumRes = zeros(axes(KS.ib.wL))
    sumAvg = zeros(axes(KS.ib.wL))

    @inbounds Threads.@threads for i = 2:KS.pSpace.nx-1
        step!(KS, face[i], ctr[i], face[i+1], dt, sumRes, sumAvg, coll, isMHD)
    end

    for i in eachindex(residual)
        residual[i] = sqrt(sumRes[i] * KS.pSpace.nx) / (sumAvg[i] + 1.e-7)
    end

    update_boundary!(KS, ctr, face, dt, residual; coll = coll, bc = bc, isMHD = isMHD)

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
end

function update!(
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
    Y<:AbstractArray{ControlVolume1D4F,1},
    Z<:AbstractArray{Interface1D4F,1},
}

    sumRes = zeros(axes(KS.ib.wL))
    sumAvg = zeros(axes(KS.ib.wL))

    @inbounds Threads.@threads for i = 2:KS.pSpace.nx-1
        step!(KS, face[i], ctr[i], face[i+1], dt, sumRes, sumAvg, coll, isMHD)
    end

    for i in eachindex(residual)
        residual[i] = sqrt(sumRes[i] * KS.pSpace.nx) / (sumAvg[i] + 1.e-7)
    end

    update_boundary!(KS, ctr, face, dt, residual; coll = coll, bc = bc, isMHD = isMHD)

end


function update_boundary!(
    KS::X,
    ctr::Y,
    face::Z,
    dt,
    residual;
    coll::Symbol,
    bc::Symbol,
    isMHD::Bool,
) where {
    X<:AbstractSolverSet,
    Y<:AbstractArray{<:AbstractControlVolume1D,1},
    Z<:AbstractArray{<:AbstractInterface1D,1},
}

    resL = zeros(axes(KS.ib.wL))
    avgL = zeros(axes(KS.ib.wL))
    resR = zeros(axes(KS.ib.wL))
    avgR = zeros(axes(KS.ib.wL))

    if bc != :fix
        i = 1
        j = KS.pSpace.nx
        if KS.set.space[3:4] == "0f"
            step!(
                face[i].fw,
                ctr[i].w,
                ctr[i].prim,
                face[i+1].fw,
                KS.gas.γ,
                ctr[i].dx,
                resL,
                avgL,
            )
            step!(
                face[j].fw,
                ctr[j].w,
                ctr[j].prim,
                face[j+1].fw,
                KS.gas.γ,
                ctr[j].dx,
                resR,
                avgR,
            )
        elseif KS.set.space[3:4] == "1f"
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
                    ctr[i].dx,
                    dt,
                    resL,
                    avgL,
                    Symbol(KS.set.collision),
                )
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
                    ctr[j].dx,
                    dt,
                    resR,
                    avgR,
                    Symbol(KS.set.collision),
                )
            elseif KS.set.space[5:6] == "3v"
                step!(
                    face[i].fw,
                    face[i].ff,
                    ctr[i].w,
                    ctr[i].prim,
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
                    ctr[i].dx,
                    dt,
                    resL,
                    avgL,
                    Symbol(KS.set.collision),
                )
                step!(
                    face[j].fw,
                    face[j].ff,
                    ctr[j].w,
                    ctr[j].prim,
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
                    ctr[j].dx,
                    dt,
                    resL,
                    avgL,
                    Symbol(KS.set.collision),
                )
            end
        elseif KS.set.space[3:4] == "2f"
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
                ctr[i].dx,
                dt,
                resL,
                avgL,
                Symbol(KS.set.collision),
            )
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
                ctr[j].dx,
                dt,
                resL,
                avgL,
                Symbol(KS.set.collision),
            )
        elseif KS.set.space[3:4] in ["3f", "4f"] # compact form
            step!(KS, face[i], ctr[i], face[i+1], dt, resL, avgL, coll, isMHD)
            step!(KS, face[j], ctr[j], face[j+1], dt, resR, avgR, coll, isMHD)
        else
            step!(KS, face[i], ctr[i], face[i+1], dt, resL, avgL, coll)
            step!(KS, face[j], ctr[j], face[j+1], dt, resR, avgR, coll)
        end
    end

    for i in eachindex(residual)
        residual[i] += sqrt((resL[i] + resR[i]) * 2) / (avgL[i] + avgR[i] + 1.e-7)
    end

    ng = 1 - first(eachindex(KS.pSpace.x))
    if bc == :extra
        for i = 1:ng
            ctr[1-i].w .= ctr[1].w
            ctr[1-i].prim .= ctr[1].prim
            ctr[KS.pSpace.nx+i].w .= ctr[KS.pSpace.nx].w
            ctr[KS.pSpace.nx+i].prim .= ctr[KS.pSpace.nx].prim

            if KS.set.space[3:4] == "1f"
                ctr[1-i].f .= ctr[1].f
                ctr[KS.pSpace.nx+i].f .= ctr[KS.pSpace.nx].f
            elseif KS.set.space[3:4] == "2f"
                ctr[1-i].h .= ctr[1].h
                ctr[1-i].b .= ctr[1].b
                ctr[KS.pSpace.nx+i].h .= ctr[KS.pSpace.nx].h
                ctr[KS.pSpace.nx+i].b .= ctr[KS.pSpace.nx].b
            elseif KS.set.space[3:4] == "3f"
                ctr[1-i].h0 .= ctr[1].h0
                ctr[1-i].h1 .= ctr[1].h1
                ctr[1-i].h2 .= ctr[1].h2
                ctr[1-i].E .= ctr[1].E
                ctr[1-i].B .= ctr[1].B
                ctr[1-i].ϕ = deepcopy(ctr[1].ϕ)
                ctr[1-i].ψ = deepcopy(ctr[1].ψ)
                ctr[1-i].lorenz .= ctr[1].lorenz

                ctr[KS.pSpace.nx+i].h0 .= ctr[KS.pSpace.nx].h0
                ctr[KS.pSpace.nx+i].h1 .= ctr[KS.pSpace.nx].h1
                ctr[KS.pSpace.nx+i].h2 .= ctr[KS.pSpace.nx].h2
                ctr[KS.pSpace.nx+i].E .= ctr[KS.pSpace.nx].E
                ctr[KS.pSpace.nx+i].B .= ctr[KS.pSpace.nx].B
                ctr[KS.pSpace.nx+i].ϕ = deepcopy(ctr[KS.pSpace.nx].ϕ)
                ctr[KS.pSpace.nx+i].ψ = deepcopy(ctr[KS.pSpace.nx].ψ)
                ctr[KS.pSpace.nx+i].lorenz .= ctr[KS.pSpace.nx].lorenz
            elseif KS.set.space[3:4] == "4f"
                ctr[1-i].h0 .= ctr[1].h0
                ctr[1-i].h1 .= ctr[1].h1
                ctr[1-i].h2 .= ctr[1].h2
                ctr[1-i].h3 .= ctr[1].h3
                ctr[1-i].E .= ctr[1].E
                ctr[1-i].B .= ctr[1].B
                ctr[1-i].ϕ = deepcopy(ctr[1].ϕ)
                ctr[1-i].ψ = deepcopy(ctr[1].ψ)
                ctr[1-i].lorenz .= ctr[1].lorenz

                ctr[KS.pSpace.nx+i].h0 .= ctr[KS.pSpace.nx].h0
                ctr[KS.pSpace.nx+i].h1 .= ctr[KS.pSpace.nx].h1
                ctr[KS.pSpace.nx+i].h2 .= ctr[KS.pSpace.nx].h2
                ctr[KS.pSpace.nx+i].h3 .= ctr[KS.pSpace.nx].h3
                ctr[KS.pSpace.nx+i].E .= ctr[KS.pSpace.nx].E
                ctr[KS.pSpace.nx+i].B .= ctr[KS.pSpace.nx].B
                ctr[KS.pSpace.nx+i].ϕ = deepcopy(ctr[KS.pSpace.nx].ϕ)
                ctr[KS.pSpace.nx+i].ψ = deepcopy(ctr[KS.pSpace.nx].ψ)
                ctr[KS.pSpace.nx+i].lorenz .= ctr[KS.pSpace.nx].lorenz
            else
                throw("incorrect amount of distribution functions")
            end
        end
    elseif bc == :period
        for i = 1:ng
            ctr[1-i].w .= ctr[KS.pSpace.nx+1-i].w
            ctr[1-i].prim .= ctr[KS.pSpace.nx+1-i].prim
            ctr[KS.pSpace.nx+i].w .= ctr[i].w
            ctr[KS.pSpace.nx+i].prim .= ctr[i].prim

            if KS.set.space[3:4] == "1f"
                ctr[1-i].f .= ctr[KS.pSpace.nx+1-i].f
                ctr[KS.pSpace.nx+i].f .= ctr[i].f
            elseif KS.set.space[3:4] == "2f"
                ctr[1-i].h .= ctr[KS.pSpace.nx+1-i].h
                ctr[1-i].b .= ctr[KS.pSpace.nx+1-i].b
                ctr[KS.pSpace.nx+i].h .= ctr[i].h
                ctr[KS.pSpace.nx+i].b .= ctr[i].b
            elseif KS.set.space[3:4] == "3f"
                ctr[1-i].h0 .= ctr[KS.pSpace.nx+1-i].h0
                ctr[1-i].h1 .= ctr[KS.pSpace.nx+1-i].h1
                ctr[1-i].h2 .= ctr[KS.pSpace.nx+1-i].h2
                ctr[1-i].E .= ctr[KS.pSpace.nx+1-i].E
                ctr[1-i].B .= ctr[KS.pSpace.nx+1-i].B
                ctr[1-i].ϕ = deepcopy(ctr[KS.pSpace.nx+1-i].ϕ)
                ctr[1-i].ψ = deepcopy(ctr[KS.pSpace.nx+1-i].ψ)
                ctr[1-i].lorenz .= ctr[KS.pSpace.nx+1-i].lorenz

                ctr[KS.pSpace.nx+i].h0 .= ctr[i].h0
                ctr[KS.pSpace.nx+i].h1 .= ctr[i].h1
                ctr[KS.pSpace.nx+i].h2 .= ctr[i].h2
                ctr[KS.pSpace.nx+i].E .= ctr[i].E
                ctr[KS.pSpace.nx+i].B .= ctr[i].B
                ctr[KS.pSpace.nx+i].ϕ = deepcopy(ctr[i].ϕ)
                ctr[KS.pSpace.nx+i].ψ = deepcopy(ctr[i].ψ)
                ctr[KS.pSpace.nx+i].lorenz .= ctr[i].lorenz
            elseif KS.set.space[3:4] == "4f"
                ctr[1-i].h0 .= ctr[KS.pSpace.nx+1-i].h0
                ctr[1-i].h1 .= ctr[KS.pSpace.nx+1-i].h1
                ctr[1-i].h2 .= ctr[KS.pSpace.nx+1-i].h2
                ctr[1-i].h3 .= ctr[KS.pSpace.nx+1-i].h3
                ctr[1-i].E .= ctr[KS.pSpace.nx+1-i].E
                ctr[1-i].B .= ctr[KS.pSpace.nx+1-i].B
                ctr[1-i].ϕ = deepcopy(ctr[KS.pSpace.nx+1-i].ϕ)
                ctr[1-i].ψ = deepcopy(ctr[KS.pSpace.nx+1-i].ψ)
                ctr[1-i].lorenz .= ctr[KS.pSpace.nx+1-i].lorenz

                ctr[KS.pSpace.nx+i].h0 .= ctr[i].h0
                ctr[KS.pSpace.nx+i].h1 .= ctr[i].h1
                ctr[KS.pSpace.nx+i].h2 .= ctr[i].h2
                ctr[KS.pSpace.nx+i].h3 .= ctr[i].h3
                ctr[KS.pSpace.nx+i].E .= ctr[i].E
                ctr[KS.pSpace.nx+i].B .= ctr[i].B
                ctr[KS.pSpace.nx+i].ϕ = deepcopy(ctr[i].ϕ)
                ctr[KS.pSpace.nx+i].ψ = deepcopy(ctr[i].ψ)
                ctr[KS.pSpace.nx+i].lorenz .= ctr[i].lorenz
            else
                throw("incorrect amount of distribution functions")
            end
        end
    elseif bc == :balance
        @. ctr[0].w = 0.5 * (ctr[-1].w + ctr[1].w)
        @. ctr[0].prim = 0.5 * (ctr[-1].prim + ctr[1].prim)
        @. ctr[KS.pSpace.nx+1].w = 0.5 * (ctr[KS.pSpace.nx].w + ctr[KS.pSpace.nx+2].w)
        @. ctr[KS.pSpace.nx+1].prim =
            0.5 * (ctr[KS.pSpace.nx].prim + ctr[KS.pSpace.nx+2].prim)

        if KS.set.space[3:4] == "1f"
            @. ctr[0].f = 0.5 * (ctr[-1].f + ctr[1].f)
            @. ctr[KS.pSpace.nx+1].f = 0.5 * (ctr[KS.pSpace.nx].f + ctr[KS.pSpace.nx+2].f)
        elseif KS.set.space[3:4] == "2f"
            @. ctr[0].h = 0.5 * (ctr[-1].h + ctr[1].h)
            @. ctr[0].b = 0.5 * (ctr[-1].b + ctr[1].b)
            @. ctr[KS.pSpace.nx+1].h = 0.5 * (ctr[KS.pSpace.nx].h + ctr[KS.pSpace.nx+2].h)
            @. ctr[KS.pSpace.nx+1].b = 0.5 * (ctr[KS.pSpace.nx].b + ctr[KS.pSpace.nx+2].b)
        elseif KS.set.space[3:4] == "3f"
            @. ctr[0].h0 = 0.5 * (ctr[-1].h0 + ctr[1].h0)
            @. ctr[0].h1 = 0.5 * (ctr[-1].h1 + ctr[1].h1)
            @. ctr[0].h2 = 0.5 * (ctr[-1].h2 + ctr[1].h2)
            @. ctr[0].E = 0.5 * (ctr[-1].E + ctr[1].E)
            @. ctr[0].B = 0.5 * (ctr[-1].B + ctr[1].B)
            ctr[0].ϕ = 0.5 * (ctr[-1].ϕ + ctr[1].ϕ)
            ctr[0].ψ = 0.5 * (ctr[-1].ψ + ctr[1].ψ)
            @. ctr[0].lorenz = 0.5 * (ctr[-1].lorenz + ctr[1].lorenz)

            @. ctr[KS.pSpace.nx+1].h0 =
                0.5 * (ctr[KS.pSpace.nx].h0 + ctr[KS.pSpace.nx+2].h0)
            @. ctr[KS.pSpace.nx+1].h1 =
                0.5 * (ctr[KS.pSpace.nx].h1 + ctr[KS.pSpace.nx+2].h1)
            @. ctr[KS.pSpace.nx+1].h2 =
                0.5 * (ctr[KS.pSpace.nx].h2 + ctr[KS.pSpace.nx+2].h2)
            @. ctr[KS.pSpace.nx+1].E = 0.5 * (ctr[KS.pSpace.nx].E + ctr[KS.pSpace.nx+2].E)
            @. ctr[KS.pSpace.nx+1].B = 0.5 * (ctr[KS.pSpace.nx].B + ctr[KS.pSpace.nx+2].B)
            ctr[KS.pSpace.nx+1].ϕ = 0.5 * (ctr[KS.pSpace.nx].ϕ + ctr[KS.pSpace.nx+2].ϕ)
            ctr[KS.pSpace.nx+1].ψ = 0.5 * (ctr[KS.pSpace.nx].ψ + ctr[KS.pSpace.nx+2].ψ)
            @. ctr[KS.pSpace.nx+1].lorenz =
                0.5 * (ctr[KS.pSpace.nx].lorenz + ctr[KS.pSpace.nx+2].lorenz)
        elseif KS.set.space[3:4] == "4f"
            @. ctr[0].h0 = 0.5 * (ctr[-1].h0 + ctr[1].h0)
            @. ctr[0].h1 = 0.5 * (ctr[-1].h1 + ctr[1].h1)
            @. ctr[0].h2 = 0.5 * (ctr[-1].h2 + ctr[1].h2)
            @. ctr[0].h3 = 0.5 * (ctr[-1].h3 + ctr[1].h3)
            @. ctr[0].E = 0.5 * (ctr[-1].E + ctr[1].E)
            @. ctr[0].B = 0.5 * (ctr[-1].B + ctr[1].B)
            ctr[0].ϕ = 0.5 * (ctr[-1].ϕ + ctr[1].ϕ)
            ctr[0].ψ = 0.5 * (ctr[-1].ψ + ctr[1].ψ)
            @. ctr[0].lorenz = 0.5 * (ctr[-1].lorenz + ctr[1].lorenz)

            @. ctr[KS.pSpace.nx+1].h0 =
                0.5 * (ctr[KS.pSpace.nx].h0 + ctr[KS.pSpace.nx+2].h0)
            @. ctr[KS.pSpace.nx+1].h1 =
                0.5 * (ctr[KS.pSpace.nx].h1 + ctr[KS.pSpace.nx+2].h1)
            @. ctr[KS.pSpace.nx+1].h2 =
                0.5 * (ctr[KS.pSpace.nx].h2 + ctr[KS.pSpace.nx+2].h2)
            @. ctr[KS.pSpace.nx+1].h3 =
                0.5 * (ctr[KS.pSpace.nx].h3 + ctr[KS.pSpace.nx+2].h3)
            @. ctr[KS.pSpace.nx+1].E = 0.5 * (ctr[KS.pSpace.nx].E + ctr[KS.pSpace.nx+2].E)
            @. ctr[KS.pSpace.nx+1].B = 0.5 * (ctr[KS.pSpace.nx].B + ctr[KS.pSpace.nx+2].B)
            ctr[KS.pSpace.nx+1].ϕ = 0.5 * (ctr[KS.pSpace.nx].ϕ + ctr[KS.pSpace.nx+2].ϕ)
            ctr[KS.pSpace.nx+1].ψ = 0.5 * (ctr[KS.pSpace.nx].ψ + ctr[KS.pSpace.nx+2].ψ)
            @. ctr[KS.pSpace.nx+1].lorenz =
                0.5 * (ctr[KS.pSpace.nx].lorenz + ctr[KS.pSpace.nx+2].lorenz)
        else
            throw("incorrect amount of distribution functions")
        end
    else
    end

end
