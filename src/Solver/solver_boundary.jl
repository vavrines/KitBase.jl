"""
Update solver for boundary cells

"""
function update_boundary!(
    KS::X,
    ctr::Y,
    face::Z,
    dt,
    residual;
    bc::Symbol,
) where {
    X<:AbstractSolverSet,
    Y<:AbstractArray{ControlVolume1D,1},
    Z<:AbstractArray{Interface1D,1},
}

    resL = zero(KS.ib.wL)
    avgL = zero(KS.ib.wL)
    resR = zero(KS.ib.wL)
    avgR = zero(KS.ib.wL)

    if bc != :fix
        i = 1
        j = KS.pSpace.nx
        if KS.set.nSpecies == 1
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
                ctr[i].dx,
                dt,
                resL,
                avgL,
            )
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
                ctr[j].dx,
                dt,
                resR,
                avgR,
            )
        end
    end

    @. residual += sqrt((resL + resR) * 2) / (avgL + avgR + 1.e-7)

    ng = 1 - first(eachindex(KS.pSpace.x))
    if bc == :extra
        for i = 1:ng
            ctr[1-i].w .= ctr[1].w
            ctr[1-i].prim .= ctr[1].prim
            ctr[KS.pSpace.nx+i].w .= ctr[KS.pSpace.nx].w
            ctr[KS.pSpace.nx+i].prim .= ctr[KS.pSpace.nx].prim
        end
    elseif bc == :period
        for i = 1:ng
            ctr[1-i].w .= ctr[KS.pSpace.nx+1-i].w
            ctr[1-i].prim .= ctr[KS.pSpace.nx+1-i].prim
            ctr[KS.pSpace.nx+i].w .= ctr[i].w
            ctr[KS.pSpace.nx+i].prim .= ctr[i].prim
        end
    elseif bc == :balance
        @. ctr[0].w = 0.5 * (ctr[-1].w + ctr[1].w)
        @. ctr[0].prim = 0.5 * (ctr[-1].prim + ctr[1].prim)
        @. ctr[KS.pSpace.nx+1].w = 0.5 * (ctr[KS.pSpace.nx].w + ctr[KS.pSpace.nx+2].w)
        @. ctr[KS.pSpace.nx+1].prim =
            0.5 * (ctr[KS.pSpace.nx].prim + ctr[KS.pSpace.nx+2].prim)
    end

end

function update_boundary!(
    KS::X,
    ctr::Y,
    face::Z,
    dt,
    residual;
    coll::Symbol,
    bc::Symbol,
    isMHD = false::Bool,
) where {
    X<:AbstractSolverSet,
    Y<:AbstractArray{ControlVolume1D1F,1},
    Z<:AbstractArray{Interface1D1F,1},
}

    resL = zero(KS.ib.wL)
    avgL = zero(KS.ib.wL)
    resR = zero(KS.ib.wL)
    avgR = zero(KS.ib.wL)

    if bc != :fix
        i = 1
        j = KS.pSpace.nx

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
            ctr[1-i].f .= ctr[1].f
            ctr[KS.pSpace.nx+i].f .= ctr[KS.pSpace.nx].f
        end
    elseif bc == :period
        for i = 1:ng
            ctr[1-i].w .= ctr[KS.pSpace.nx+1-i].w
            ctr[1-i].prim .= ctr[KS.pSpace.nx+1-i].prim
            ctr[KS.pSpace.nx+i].w .= ctr[i].w
            ctr[KS.pSpace.nx+i].prim .= ctr[i].prim
            ctr[1-i].f .= ctr[KS.pSpace.nx+1-i].f
            ctr[KS.pSpace.nx+i].f .= ctr[i].f
        end
    elseif bc == :balance
        @. ctr[0].w = 0.5 * (ctr[-1].w + ctr[1].w)
        @. ctr[0].prim = 0.5 * (ctr[-1].prim + ctr[1].prim)
        @. ctr[KS.pSpace.nx+1].w = 0.5 * (ctr[KS.pSpace.nx].w + ctr[KS.pSpace.nx+2].w)
        @. ctr[KS.pSpace.nx+1].prim =
            0.5 * (ctr[KS.pSpace.nx].prim + ctr[KS.pSpace.nx+2].prim)
        @. ctr[0].f = 0.5 * (ctr[-1].f + ctr[1].f)
        @. ctr[KS.pSpace.nx+1].f = 0.5 * (ctr[KS.pSpace.nx].f + ctr[KS.pSpace.nx+2].f)
    end

end

function update_boundary!(
    KS::X,
    ctr::Y,
    face::Z,
    dt,
    residual;
    coll::Symbol,
    bc::Symbol,
    isMHD = false::Bool,
) where {
    X<:AbstractSolverSet,
    Y<:AbstractArray{ControlVolume1D2F,1},
    Z<:AbstractArray{Interface1D2F,1},
}

    resL = zero(KS.ib.wL)
    avgL = zero(KS.ib.wL)
    resR = zero(KS.ib.wL)
    avgR = zero(KS.ib.wL)

    if bc != :fix
        i = 1
        j = KS.pSpace.nx
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
                KS.gas.mi,
                KS.gas.ni,
                KS.gas.me,
                KS.gas.ne,
                KS.gas.Kn[1],
                KS.gas.Pr,
                ctr[j].dx,
                dt,
                resL,
                avgL,
                Symbol(KS.set.collision),
            )
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

            ctr[1-i].h .= ctr[1].h
            ctr[1-i].b .= ctr[1].b
            ctr[KS.pSpace.nx+i].h .= ctr[KS.pSpace.nx].h
            ctr[KS.pSpace.nx+i].b .= ctr[KS.pSpace.nx].b
        end
    elseif bc == :period
        for i = 1:ng
            ctr[1-i].w .= ctr[KS.pSpace.nx+1-i].w
            ctr[1-i].prim .= ctr[KS.pSpace.nx+1-i].prim
            ctr[KS.pSpace.nx+i].w .= ctr[i].w
            ctr[KS.pSpace.nx+i].prim .= ctr[i].prim

            ctr[1-i].h .= ctr[KS.pSpace.nx+1-i].h
            ctr[1-i].b .= ctr[KS.pSpace.nx+1-i].b
            ctr[KS.pSpace.nx+i].h .= ctr[i].h
            ctr[KS.pSpace.nx+i].b .= ctr[i].b
        end
    elseif bc == :balance
        @. ctr[0].w = 0.5 * (ctr[-1].w + ctr[1].w)
        @. ctr[0].prim = 0.5 * (ctr[-1].prim + ctr[1].prim)
        @. ctr[KS.pSpace.nx+1].w = 0.5 * (ctr[KS.pSpace.nx].w + ctr[KS.pSpace.nx+2].w)
        @. ctr[KS.pSpace.nx+1].prim =
            0.5 * (ctr[KS.pSpace.nx].prim + ctr[KS.pSpace.nx+2].prim)

        @. ctr[0].h = 0.5 * (ctr[-1].h + ctr[1].h)
        @. ctr[0].b = 0.5 * (ctr[-1].b + ctr[1].b)
        @. ctr[KS.pSpace.nx+1].h = 0.5 * (ctr[KS.pSpace.nx].h + ctr[KS.pSpace.nx+2].h)
        @. ctr[KS.pSpace.nx+1].b = 0.5 * (ctr[KS.pSpace.nx].b + ctr[KS.pSpace.nx+2].b)
    end

end

function update_boundary!(
    KS::X,
    ctr::Y,
    face::Z,
    dt,
    residual;
    coll::Symbol,
    bc::Symbol,
    isMHD = false::Bool,
) where {
    X<:AbstractSolverSet,
    Y<:AbstractArray{ControlVolume1D3F,1},
    Z<:AbstractArray{Interface1D3F,1},
}

    resL = zero(KS.ib.wL)
    avgL = zero(KS.ib.wL)
    resR = zero(KS.ib.wL)
    avgR = zero(KS.ib.wL)

    if bc != :fix
        i = 1
        j = KS.pSpace.nx
        
        step!(KS, face[i], ctr[i], face[i+1], dt, resL, avgL, coll, isMHD)
        step!(KS, face[j], ctr[j], face[j+1], dt, resR, avgR, coll, isMHD)
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
        end
    elseif bc == :period
        for i = 1:ng
            ctr[1-i].w .= ctr[KS.pSpace.nx+1-i].w
            ctr[1-i].prim .= ctr[KS.pSpace.nx+1-i].prim
            ctr[KS.pSpace.nx+i].w .= ctr[i].w
            ctr[KS.pSpace.nx+i].prim .= ctr[i].prim

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
        end
    elseif bc == :balance
        @. ctr[0].w = 0.5 * (ctr[-1].w + ctr[1].w)
        @. ctr[0].prim = 0.5 * (ctr[-1].prim + ctr[1].prim)
        @. ctr[KS.pSpace.nx+1].w = 0.5 * (ctr[KS.pSpace.nx].w + ctr[KS.pSpace.nx+2].w)
        @. ctr[KS.pSpace.nx+1].prim =
            0.5 * (ctr[KS.pSpace.nx].prim + ctr[KS.pSpace.nx+2].prim)

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
    end

end

function update_boundary!(
    KS::X,
    ctr::Y,
    face::Z,
    dt,
    residual;
    coll::Symbol,
    bc::Symbol,
    isMHD = false::Bool,
) where {
    X<:AbstractSolverSet,
    Y<:AbstractArray{ControlVolume1D4F,1},
    Z<:AbstractArray{Interface1D4F,1},
}

    resL = zeros(axes(KS.ib.wL))
    avgL = zeros(axes(KS.ib.wL))
    resR = zeros(axes(KS.ib.wL))
    avgR = zeros(axes(KS.ib.wL))

    if bc != :fix
        i = 1
        j = KS.pSpace.nx
        
        step!(KS, face[i], ctr[i], face[i+1], dt, resL, avgL, coll, isMHD)
        step!(KS, face[j], ctr[j], face[j+1], dt, resR, avgR, coll, isMHD)
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
        end
    elseif bc == :period
        for i = 1:ng
            ctr[1-i].w .= ctr[KS.pSpace.nx+1-i].w
            ctr[1-i].prim .= ctr[KS.pSpace.nx+1-i].prim
            ctr[KS.pSpace.nx+i].w .= ctr[i].w
            ctr[KS.pSpace.nx+i].prim .= ctr[i].prim

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
        end
    elseif bc == :balance
        @. ctr[0].w = 0.5 * (ctr[-1].w + ctr[1].w)
        @. ctr[0].prim = 0.5 * (ctr[-1].prim + ctr[1].prim)
        @. ctr[KS.pSpace.nx+1].w = 0.5 * (ctr[KS.pSpace.nx].w + ctr[KS.pSpace.nx+2].w)
        @. ctr[KS.pSpace.nx+1].prim =
            0.5 * (ctr[KS.pSpace.nx].prim + ctr[KS.pSpace.nx+2].prim)

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
    end

end

function update_boundary!(
    KS::X,
    ctr::Y,
    a1face::Z,
    a2face::Z,
    dt,
    residual;
    coll::Symbol,
    bc::Symbol,
) where {
    X<:AbstractSolverSet,
    Y<:AbstractArray{ControlVolume2D1F,2},
    Z<:AbstractArray{Interface2D1F,2},
}

    resL = zero(KS.ib.wL)
    avgL = zero(KS.ib.wL)
    resR = zero(KS.ib.wL)
    avgR = zero(KS.ib.wL)
    resU = zero(KS.ib.wL)
    avgU = zero(KS.ib.wL)
    resD = zero(KS.ib.wL)
    avgD = zero(KS.ib.wL)

    if bc != :fix
        @inbounds for j = 1:KS.pSpace.ny
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
                ctr[1, j].dx * ctr[1, j].dy,
                dt,
                resL,
                avgL,
                coll,
            )

            step!(
                ctr[KS.pSpace.nx, j].w,
                ctr[KS.pSpace.nx, j].prim,
                ctr[KS.pSpace.nx, j].f,
                a1face[KS.pSpace.nx, j].fw,
                a1face[KS.pSpace.nx, j].ff,
                a1face[KS.pSpace.nx+1, j].fw,
                a1face[KS.pSpace.nx+1, j].ff,
                a2face[KS.pSpace.nx, j].fw,
                a2face[KS.pSpace.nx, j].ff,
                a2face[KS.pSpace.nx, j+1].fw,
                a2face[KS.pSpace.nx, j+1].ff,
                KS.vSpace.u,
                KS.vSpace.v,
                KS.vSpace.weights,
                KS.gas.γ,
                KS.gas.μᵣ,
                KS.gas.ω,
                KS.gas.Pr,
                ctr[KS.pSpace.nx, j].dx * ctr[KS.pSpace.nx, j].dy,
                dt,
                resR,
                avgR,
                coll,
            )
        end

        @inbounds for i = 2:KS.pSpace.nx-1 # skip overlap
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
                ctr[i, 1].dx * ctr[i, 1].dy,
                dt,
                resD,
                avgD,
                coll,
            )

            step!(
                ctr[i, KS.pSpace.ny].w,
                ctr[i, KS.pSpace.ny].prim,
                ctr[i, KS.pSpace.ny].f,
                a1face[i, KS.pSpace.ny].fw,
                a1face[i, KS.pSpace.ny].ff,
                a1face[i+1, KS.pSpace.ny].fw,
                a1face[i+1, KS.pSpace.ny].ff,
                a2face[i, KS.pSpace.ny].fw,
                a2face[i, KS.pSpace.ny].ff,
                a2face[i, KS.pSpace.ny+1].fw,
                a2face[i, KS.pSpace.ny+1].ff,
                KS.vSpace.u,
                KS.vSpace.v,
                KS.vSpace.weights,
                KS.gas.γ,
                KS.gas.μᵣ,
                KS.gas.ω,
                KS.gas.Pr,
                ctr[i, KS.pSpace.ny].dx * ctr[i, KS.pSpace.ny].dy,
                dt,
                resU,
                avgU,
                coll,
            )
        end
    end

    for i in eachindex(residual)
        residual[i] += sqrt((resL[i] + resR[i] + resU[i] + resD[i]) * 2) / (avgL[i] + avgR[i] + avgU[i] + avgD[i] + 1.e-7)
    end

    ngx = 1 - first(eachindex(KS.pSpace.x[:, 1]))
    ngy = 1 - first(eachindex(KS.pSpace.y[1, :]))
    if bc == :extra
        for i = 1:ngx, j = 1:KS.pSpace.ny
            ctr[1-i, j].w .= ctr[1, j].w
            ctr[1-i, j].prim .= ctr[1, j].prim
            ctr[KS.pSpace.nx+i, j].w .= ctr[KS.pSpace.nx, j].w
            ctr[KS.pSpace.nx+i, j].prim .= ctr[KS.pSpace.nx, j].prim

            ctr[1-i, j].f .= ctr[1, j].f
            ctr[KS.pSpace.nx+i, j].f .= ctr[KS.pSpace.nx, j].f           
        end

        for i = 1:KS.pSpace.nx, j = 1:ngy
            ctr[i, 1-j].w .= ctr[i, 1].w
            ctr[i, 1-j].prim .= ctr[i, 1].prim
            ctr[i, KS.pSpace.ny+j].w .= ctr[i, KS.pSpace.ny].w
            ctr[i, KS.pSpace.ny+j].prim .= ctr[i, KS.pSpace.ny].prim

            ctr[i, 1-j].f .= ctr[i, 1].f
            ctr[i, KS.pSpace.ny+j].f .= ctr[i, KS.pSpace.ny].f          
        end
    elseif bc == :period

    elseif bc == :balance
        
    end

end

function update_boundary!(
    KS::X,
    ctr::Y,
    a1face::Z,
    a2face::Z,
    dt,
    residual;
    coll::Symbol,
    bc::Symbol,
) where {
    X<:AbstractSolverSet,
    Y<:AbstractArray{ControlVolume2D2F,2},
    Z<:AbstractArray{Interface2D2F,2},
}

    resL = zero(KS.ib.wL)
    avgL = zero(KS.ib.wL)
    resR = zero(KS.ib.wL)
    avgR = zero(KS.ib.wL)
    resU = zero(KS.ib.wL)
    avgU = zero(KS.ib.wL)
    resD = zero(KS.ib.wL)
    avgD = zero(KS.ib.wL)

    if bc != :fix
        @inbounds for j = 1:KS.pSpace.ny
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
                ctr[1, j].dx * ctr[1, j].dy,
                dt,
                resL,
                avgL,
                coll,
            )

            step!(
                ctr[KS.pSpace.nx, j].w,
                ctr[KS.pSpace.nx, j].prim,
                ctr[KS.pSpace.nx, j].h,
                ctr[KS.pSpace.nx, j].b,
                a1face[KS.pSpace.nx, j].fw,
                a1face[KS.pSpace.nx, j].fh,
                a1face[KS.pSpace.nx, j].fb,
                a1face[KS.pSpace.nx+1, j].fw,
                a1face[KS.pSpace.nx+1, j].fh,
                a1face[KS.pSpace.nx+1, j].fb,
                a2face[KS.pSpace.nx, j].fw,
                a2face[KS.pSpace.nx, j].fh,
                a2face[KS.pSpace.nx, j].fb,
                a2face[KS.pSpace.nx, j+1].fw,
                a2face[KS.pSpace.nx, j+1].fh,
                a2face[KS.pSpace.nx, j+1].fb,
                KS.vSpace.u,
                KS.vSpace.v,
                KS.vSpace.weights,
                KS.gas.K,
                KS.gas.γ,
                KS.gas.μᵣ,
                KS.gas.ω,
                KS.gas.Pr,
                ctr[KS.pSpace.nx, j].dx * ctr[KS.pSpace.nx, j].dy,
                dt,
                resR,
                avgR,
                coll,
            )
        end

        @inbounds for i = 2:KS.pSpace.nx-1 # skip overlap
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
                ctr[i, 1].dx * ctr[i, 1].dy,
                dt,
                resD,
                avgD,
                coll,
            )

            step!(
                ctr[i, KS.pSpace.ny].w,
                ctr[i, KS.pSpace.ny].prim,
                ctr[i, KS.pSpace.ny].h,
                ctr[i, KS.pSpace.ny].b,
                a1face[i, KS.pSpace.ny].fw,
                a1face[i, KS.pSpace.ny].fh,
                a1face[i, KS.pSpace.ny].fb,
                a1face[i+1, KS.pSpace.ny].fw,
                a1face[i+1, KS.pSpace.ny].fh,
                a1face[i+1, KS.pSpace.ny].fb,
                a2face[i, KS.pSpace.ny].fw,
                a2face[i, KS.pSpace.ny].fh,
                a2face[i, KS.pSpace.ny].fb,
                a2face[i, KS.pSpace.ny+1].fw,
                a2face[i, KS.pSpace.ny+1].fh,
                a2face[i, KS.pSpace.ny+1].fb,
                KS.vSpace.u,
                KS.vSpace.v,
                KS.vSpace.weights,
                KS.gas.K,
                KS.gas.γ,
                KS.gas.μᵣ,
                KS.gas.ω,
                KS.gas.Pr,
                ctr[i, KS.pSpace.ny].dx * ctr[i, KS.pSpace.ny].dy,
                dt,
                resU,
                avgU,
                coll,
            )
        end
    end

    for i in eachindex(residual)
        residual[i] += sqrt((resL[i] + resR[i] + resU[i] + resD[i]) * 2) / (avgL[i] + avgR[i] + avgU[i] + avgD[i] + 1.e-7)
    end

    ngx = 1 - first(eachindex(KS.pSpace.x[:, 1]))
    ngy = 1 - first(eachindex(KS.pSpace.y[1, :]))
    if bc == :extra
        for i = 1:ngx, j = 1:KS.pSpace.ny
            ctr[1-i, j].w .= ctr[1, j].w
            ctr[1-i, j].prim .= ctr[1, j].prim
            ctr[KS.pSpace.nx+i, j].w .= ctr[KS.pSpace.nx, j].w
            ctr[KS.pSpace.nx+i, j].prim .= ctr[KS.pSpace.nx, j].prim

            ctr[1-i, j].h .= ctr[1, j].h
            ctr[1-i, j].b .= ctr[1, j].b
            ctr[KS.pSpace.nx+i, j].h .= ctr[KS.pSpace.nx, j].h
            ctr[KS.pSpace.nx+i, j].b .= ctr[KS.pSpace.nx, j].b            
        end

        for i = 1:KS.pSpace.nx, j = 1:ngy
            ctr[i, 1-j].w .= ctr[i, 1].w
            ctr[i, 1-j].prim .= ctr[i, 1].prim
            ctr[i, KS.pSpace.ny+j].w .= ctr[i, KS.pSpace.ny].w
            ctr[i, KS.pSpace.ny+j].prim .= ctr[i, KS.pSpace.ny].prim

            ctr[i, 1-j].h .= ctr[i, 1].h
            ctr[i, 1-j].b .= ctr[i, 1].b
            ctr[i, KS.pSpace.ny+j].h .= ctr[i, KS.pSpace.ny].h
            ctr[i, KS.pSpace.ny+j].b .= ctr[i, KS.pSpace.ny].b            
        end
    elseif bc == :period

    elseif bc == :balance
        
    end

end
