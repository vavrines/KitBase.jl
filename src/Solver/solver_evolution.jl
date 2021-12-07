"""
    evolve!(
        KS::SolverSet,
        ctr::T1,
        face::T2,
        dt;
        mode = Symbol(KS.set.flux)::Symbol,
        bc = :fix::Symbol,
    ) where {
        T1<:AA{<:AbstractControlVolume1D,1},
        T2<:AA{<:AbstractInterface1D,1},
    }
    
    evolve!(
        KS::SolverSet,
        ctr::T1,
        a1face::T2,
        a2face::T2,
        dt;
        mode = Symbol(KS.set.flux)::Symbol,
        bc = :fix::Symbol,
    ) where {
        T1<:AA{<:AbstractControlVolume2D,2},
        T2<:AA{<:AbstractInterface2D,2},
    }

Evolution of boundary fluxes

"""
function evolve!(
    KS::SolverSet,
    ctr::T1,
    face::T2,
    dt;
    mode = symbolize(KS.set.flux)::Symbol,
    bc = symbolize(KS.set.boundary),
) where {T1<:AA{ControlVolume1D,1},T2<:AA{Interface1D,1}}

    if firstindex(KS.pSpace.x) < 1
        idx0 = 1
        idx1 = KS.pSpace.nx + 1
    else
        idx0 = 2
        idx1 = KS.pSpace.nx
    end

    if mode == :gks

        if KS.set.nSpecies == 1

            if KS.set.matter == "scalar"
                @inbounds Threads.@threads for i = idx0:idx1
                    face[i].fw = KitBase.flux_gks(
                        ctr[i-1].w + 0.5 * KS.ps.dx[i-1] * ctr[i-1].sw,
                        ctr[i].w - 0.5 * KS.ps.dx[i] * ctr[i].sw,
                        KS.gas.μᵣ,
                        dt,
                        0.5 * KS.ps.dx[i-1],
                        0.5 * KS.ps.dx[i],
                        KS.gas.a,
                        ctr[i-1].sw,
                        ctr[i].sw,
                    )
                end
            else
                @inbounds Threads.@threads for i = idx0:idx1
                    flux_gks!(
                        face[i].fw,
                        ctr[i-1].w .+ 0.5 .* KS.ps.dx[i-1] .* ctr[i-1].sw,
                        ctr[i].w .- 0.5 .* KS.ps.dx[i] .* ctr[i].sw,
                        KS.gas.γ,
                        KS.gas.K,
                        KS.gas.μᵣ,
                        KS.gas.ω,
                        dt,
                        0.5 * KS.ps.dx[i-1],
                        0.5 * KS.ps.dx[i],
                        ctr[i-1].sw,
                        ctr[i].sw,
                    )
                end
            end

        elseif KS.set.nSpecies == 2

            @inbounds Threads.@threads for i = idx0:idx1
                flux_gks!(
                    face[i].fw,
                    ctr[i-1].w .+ 0.5 .* KS.ps.dx[i-1] .* ctr[i-1].sw,
                    ctr[i].w .- 0.5 .* KS.ps.dx[i] .* ctr[i].sw,
                    KS.gas.γ,
                    KS.gas.K,
                    KS.gas.mi,
                    KS.gas.ni,
                    KS.gas.me,
                    KS.gas.ne,
                    KS.gas.Kn[1],
                    dt,
                    0.5 * KS.ps.dx[i-1],
                    0.5 * KS.ps.dx[i],
                    ctr[i-1].sw,
                    ctr[i].sw,
                )
            end

        end

    elseif mode == :roe

        @inbounds Threads.@threads for i = idx0:idx1
            flux_roe!(
                face[i].fw,
                ctr[i-1].w .+ 0.5 .* KS.ps.dx[i-1] .* ctr[i-1].sw,
                ctr[i].w .- 0.5 .* KS.ps.dx[i] .* ctr[i].sw,
                KS.gas.γ,
                dt,
            )
        end

    end

end

function evolve!(
    KS::SolverSet,
    ctr::T1,
    face::T2,
    dt;
    mode = symbolize(KS.set.flux)::Symbol,
    bc = symbolize(KS.set.boundary),
) where {T1<:AA{ControlVolume1D1F,1},T2<:AA{Interface1D1F,1}}

    if firstindex(KS.pSpace.x) < 1
        idx0 = 1
        idx1 = KS.pSpace.nx + 1
    else
        idx0 = 2
        idx1 = KS.pSpace.nx
    end

    if KS.set.space[5:end] == "1v"

        if mode == :kfvs
            @inbounds Threads.@threads for i = idx0:idx1
                flux_kfvs!(
                    face[i].fw,
                    face[i].ff,
                    ctr[i-1].f .+ 0.5 .* KS.ps.dx[i-1] .* ctr[i-1].sf,
                    ctr[i].f .- 0.5 .* KS.ps.dx[i] .* ctr[i].sf,
                    KS.vSpace.u,
                    KS.vSpace.weights,
                    dt,
                    ctr[i-1].sf,
                    ctr[i].sf,
                )
            end
        elseif mode == :kcu
            @inbounds Threads.@threads for i = idx0:idx1
                flux_kcu!(
                    face[i].fw,
                    face[i].ff,
                    ctr[i-1].w .+ 0.5 .* KS.ps.dx[i-1] .* ctr[i-1].sw,
                    ctr[i-1].f .+ 0.5 .* KS.ps.dx[i-1] .* ctr[i-1].sf,
                    ctr[i].w .- 0.5 .* KS.ps.dx[i] .* ctr[i].sw,
                    ctr[i].f .- 0.5 .* KS.ps.dx[i] .* ctr[i].sf,
                    KS.vSpace.u,
                    KS.vSpace.weights,
                    KS.gas.K,
                    KS.gas.γ,
                    KS.gas.μᵣ,
                    KS.gas.ω,
                    KS.gas.Pr,
                    dt,
                )
            end
        end

    elseif KS.set.space[5:end] == "3v"

        if mode == :kfvs
            @inbounds Threads.@threads for i = idx0:idx1
                flux_kfvs!(
                    face[i].fw,
                    face[i].ff,
                    ctr[i-1].f .+ 0.5 .* KS.ps.dx[i-1] .* ctr[i-1].sf,
                    ctr[i].f .- 0.5 .* KS.ps.dx[i] .* ctr[i].sf,
                    KS.vSpace.u,
                    KS.vSpace.v,
                    KS.vSpace.w,
                    KS.vSpace.weights,
                    dt,
                    ctr[i-1].sf,
                    ctr[i].sf,
                )
            end
        elseif mode == :kcu
            @inbounds Threads.@threads for i = idx0:idx1
                flux_kcu!(
                    face[i].fw,
                    face[i].ff,
                    ctr[i-1].w .+ 0.5 .* KS.ps.dx[i-1] .* ctr[i-1].sw,
                    ctr[i-1].f .+ 0.5 .* KS.ps.dx[i-1] .* ctr[i-1].sf,
                    ctr[i].w .- 0.5 .* KS.ps.dx[i] .* ctr[i].sw,
                    ctr[i].f .- 0.5 .* KS.ps.dx[i] .* ctr[i].sf,
                    KS.vSpace.u,
                    KS.vSpace.v,
                    KS.vSpace.w,
                    KS.vSpace.weights,
                    KS.gas.K,
                    KS.gas.γ,
                    KS.gas.μᵣ,
                    KS.gas.ω,
                    KS.gas.Pr,
                    dt,
                )
            end
        end

    end

end

function evolve!(
    KS::SolverSet,
    ctr::T1,
    face::T2,
    dt;
    mode = symbolize(KS.set.flux)::Symbol,
    bc = symbolize(KS.set.boundary),
) where {T1<:AA{ControlVolume1D2F,1},T2<:AA{Interface1D2F,1}}

    if firstindex(KS.pSpace.x) < 1
        idx0 = 1
        idx1 = KS.pSpace.nx + 1
    else
        idx0 = 2
        idx1 = KS.pSpace.nx
    end

    if mode == :kfvs

        @inbounds Threads.@threads for i = idx0:idx1
            flux_kfvs!(
                face[i].fw,
                face[i].fh,
                face[i].fb,
                ctr[i-1].h .+ 0.5 .* KS.ps.dx[i-1] .* ctr[i-1].sh,
                ctr[i-1].b .+ 0.5 .* KS.ps.dx[i-1] .* ctr[i-1].sb,
                ctr[i].h .- 0.5 .* KS.ps.dx[i] .* ctr[i].sh,
                ctr[i].b .- 0.5 .* KS.ps.dx[i] .* ctr[i].sb,
                KS.vSpace.u,
                KS.vSpace.weights,
                dt,
                ctr[i-1].sh,
                ctr[i-1].sb,
                ctr[i].sh,
                ctr[i].sb,
            )
        end

    elseif mode == :kcu

        if KS.set.nSpecies == 1
            @inbounds Threads.@threads for i = idx0:idx1
                flux_kcu!(
                    face[i].fw,
                    face[i].fh,
                    face[i].fb,
                    ctr[i-1].w .+ 0.5 .* KS.ps.dx[i-1] .* ctr[i-1].sw,
                    ctr[i-1].h .+ 0.5 .* KS.ps.dx[i-1] .* ctr[i-1].sh,
                    ctr[i-1].b .+ 0.5 .* KS.ps.dx[i-1] .* ctr[i-1].sb,
                    ctr[i].w .- 0.5 .* KS.ps.dx[i] .* ctr[i].sw,
                    ctr[i].h .- 0.5 .* KS.ps.dx[i] .* ctr[i].sh,
                    ctr[i].b .- 0.5 .* KS.ps.dx[i] .* ctr[i].sb,
                    KS.vSpace.u,
                    KS.vSpace.weights,
                    KS.gas.K,
                    KS.gas.γ,
                    KS.gas.μᵣ,
                    KS.gas.ω,
                    KS.gas.Pr,
                    dt,
                )
            end
        elseif KS.set.nSpecies == 2
            @inbounds Threads.@threads for i = idx0:idx1
                flux_kcu!(
                    face[i].fw,
                    face[i].fh,
                    face[i].fb,
                    ctr[i-1].w .+ 0.5 .* KS.ps.dx[i-1] .* ctr[i-1].sw,
                    ctr[i-1].h .+ 0.5 .* KS.ps.dx[i-1] .* ctr[i-1].sh,
                    ctr[i-1].b .+ 0.5 .* KS.ps.dx[i-1] .* ctr[i-1].sb,
                    ctr[i].w .- 0.5 .* KS.ps.dx[i] .* ctr[i].sw,
                    ctr[i].h .- 0.5 .* KS.ps.dx[i] .* ctr[i].sh,
                    ctr[i].b .- 0.5 .* KS.ps.dx[i] .* ctr[i].sb,
                    KS.vSpace.u,
                    KS.vSpace.weights,
                    KS.gas.K,
                    KS.gas.γ,
                    KS.gas.mi,
                    KS.gas.ni,
                    KS.gas.me,
                    KS.gas.ne,
                    KS.gas.Kn[1],
                    dt,
                )
            end
        end

    end

    bcs = ifelse(bc isa Symbol, [bc, bc], bc)
    if bcs[1] == :maxwell
        flux_boundary_maxwell!(
            face[1].fw,
            face[1].fh,
            face[1].fb,
            KS.ib.bc(KS.ps.x0),
            ctr[1].h,
            ctr[1].b,
            KS.vSpace.u,
            KS.vSpace.weights,
            KS.gas.inK,
            dt,
            1,
        )
    end
    if bcs[2] == :maxwell
        flux_boundary_maxwell!(
            face[KS.pSpace.nx+1].fw,
            face[KS.pSpace.nx+1].fh,
            face[KS.pSpace.nx+1].fb,
            KS.ib.bc(KS.ps.x1),
            ctr[KS.pSpace.nx].h,
            ctr[KS.pSpace.nx].b,
            KS.vSpace.u,
            KS.vSpace.weights,
            KS.gas.inK,
            dt,
            -1,
        )
    end

end

function evolve!(
    KS::SolverSet,
    ctr::T1,
    face::T2,
    dt;
    mode = symbolize(KS.set.flux)::Symbol,
    bc = symbolize(KS.set.boundary),
    isPlasma = false::Bool,
    isMHD = false::Bool,
) where {T1<:AA{ControlVolume1D4F,1},T2<:AA{Interface1D4F,1}}

    if firstindex(KS.pSpace.x) < 1
        idx0 = 1
        idx1 = KS.pSpace.nx + 1
    else
        idx0 = 2
        idx1 = KS.pSpace.nx
    end

    if mode == :kcu
        @inbounds Threads.@threads for i = idx0:idx1
            flux_kcu!(
                face[i].fw,
                face[i].fh0,
                face[i].fh1,
                face[i].fh2,
                face[i].fh3,
                ctr[i-1].w .+ 0.5 .* KS.ps.dx[i-1] .* ctr[i-1].sw,
                ctr[i-1].h0 .+ 0.5 .* KS.ps.dx[i-1] .* ctr[i-1].sh0,
                ctr[i-1].h1 .+ 0.5 .* KS.ps.dx[i-1] .* ctr[i-1].sh1,
                ctr[i-1].h2 .+ 0.5 .* KS.ps.dx[i-1] .* ctr[i-1].sh2,
                ctr[i-1].h3 .+ 0.5 .* KS.ps.dx[i-1] .* ctr[i-1].sh3,
                ctr[i].w .- 0.5 .* KS.ps.dx[i] .* ctr[i].sw,
                ctr[i].h0 .- 0.5 .* KS.ps.dx[i] .* ctr[i].sh0,
                ctr[i].h1 .- 0.5 .* KS.ps.dx[i] .* ctr[i].sh1,
                ctr[i].h2 .- 0.5 .* KS.ps.dx[i] .* ctr[i].sh2,
                ctr[i].h3 .- 0.5 .* KS.ps.dx[i] .* ctr[i].sh3,
                KS.vSpace.u,
                KS.vSpace.weights,
                KS.gas.K,
                KS.gas.γ,
                KS.gas.mi,
                KS.gas.ni,
                KS.gas.me,
                KS.gas.ne,
                KS.gas.Kn[1],
                dt,
                isMHD,
            )
        end
    elseif mode == :kfvs
        @inbounds Threads.@threads for i = idx0:idx1
            flux_kfvs!(
                face[i].fw,
                face[i].fh0,
                face[i].fh1,
                face[i].fh2,
                face[i].fh3,
                ctr[i-1].h0 .+ 0.5 .* KS.ps.dx[i-1] .* ctr[i-1].sh0,
                ctr[i-1].h1 .+ 0.5 .* KS.ps.dx[i-1] .* ctr[i-1].sh1,
                ctr[i-1].h2 .+ 0.5 .* KS.ps.dx[i-1] .* ctr[i-1].sh2,
                ctr[i-1].h3 .+ 0.5 .* KS.ps.dx[i-1] .* ctr[i-1].sh3,
                ctr[i].h0 .- 0.5 .* KS.ps.dx[i] .* ctr[i].sh0,
                ctr[i].h1 .- 0.5 .* KS.ps.dx[i] .* ctr[i].sh1,
                ctr[i].h2 .- 0.5 .* KS.ps.dx[i] .* ctr[i].sh2,
                ctr[i].h3 .- 0.5 .* KS.ps.dx[i] .* ctr[i].sh3,
                KS.vSpace.u,
                KS.vSpace.weights,
                dt,
                ctr[i-1].sh0,
                ctr[i-1].sh1,
                ctr[i-1].sh2,
                ctr[i-1].sh3,
                ctr[i].sh0,
                ctr[i].sh1,
                ctr[i].sh2,
                ctr[i].sh3,
            )
        end
    end

    if isPlasma
        @inbounds Threads.@threads for i = idx0:idx1
            flux_em!(
                face[i].femL,
                face[i].femR,
                ctr[i-2].E,
                ctr[i-2].B,
                ctr[i-1].E,
                ctr[i-1].B,
                ctr[i].E,
                ctr[i].B,
                ctr[i+1].E,
                ctr[i+1].B,
                ctr[i-1].ϕ,
                ctr[i].ϕ,
                ctr[i-1].ψ,
                ctr[i].ψ,
                KS.ps.dx[i-1],
                KS.ps.dx[i],
                KS.gas.Ap,
                KS.gas.An,
                KS.gas.D,
                KS.gas.sol,
                KS.gas.χ,
                KS.gas.ν,
                dt,
            )
        end
    end

end

function evolve!(
    KS::SolverSet,
    ctr::T1,
    face::T2,
    dt;
    mode = symbolize(KS.set.flux)::Symbol,
    bc = symbolize(KS.set.boundary),
    isPlasma = false::Bool,
    isMHD = false::Bool,
) where {T1<:AA{ControlVolume1D3F,1},T2<:AA{Interface1D3F,1}}

    if firstindex(KS.pSpace.x) < 1
        idx0 = 1
        idx1 = KS.pSpace.nx + 1
    else
        idx0 = 2
        idx1 = KS.pSpace.nx
    end

    if mode == :kcu
        @inbounds Threads.@threads for i = idx0:idx1
            flux_kcu!(
                face[i].fw,
                face[i].fh0,
                face[i].fh1,
                face[i].fh2,
                ctr[i-1].w .+ 0.5 .* KS.ps.dx[i-1] .* ctr[i-1].sw,
                ctr[i-1].h0 .+ 0.5 .* KS.ps.dx[i-1] .* ctr[i-1].sh0,
                ctr[i-1].h1 .+ 0.5 .* KS.ps.dx[i-1] .* ctr[i-1].sh1,
                ctr[i-1].h2 .+ 0.5 .* KS.ps.dx[i-1] .* ctr[i-1].sh2,
                ctr[i].w .- 0.5 .* KS.ps.dx[i] .* ctr[i].sw,
                ctr[i].h0 .- 0.5 .* KS.ps.dx[i] .* ctr[i].sh0,
                ctr[i].h1 .- 0.5 .* KS.ps.dx[i] .* ctr[i].sh1,
                ctr[i].h2 .- 0.5 .* KS.ps.dx[i] .* ctr[i].sh2,
                KS.vSpace.u,
                KS.vSpace.v,
                KS.vSpace.weights,
                KS.gas.K,
                KS.gas.γ,
                KS.gas.mi,
                KS.gas.ni,
                KS.gas.me,
                KS.gas.ne,
                KS.gas.Kn[1],
                dt,
                1.0,
                isMHD,
            )
        end
    elseif mode == :kfvs
        @inbounds Threads.@threads for i = idx0:idx1
            flux_kfvs!(
                face[i].fw,
                face[i].fh0,
                face[i].fh1,
                face[i].fh2,
                ctr[i-1].h0 .+ 0.5 .* KS.ps.dx[i-1] .* ctr[i-1].sh0,
                ctr[i-1].h1 .+ 0.5 .* KS.ps.dx[i-1] .* ctr[i-1].sh1,
                ctr[i-1].h2 .+ 0.5 .* KS.ps.dx[i-1] .* ctr[i-1].sh2,
                ctr[i].h0 .- 0.5 .* KS.ps.dx[i] .* ctr[i].sh0,
                ctr[i].h1 .- 0.5 .* KS.ps.dx[i] .* ctr[i].sh1,
                ctr[i].h2 .- 0.5 .* KS.ps.dx[i] .* ctr[i].sh2,
                KS.vSpace.u,
                KS.vSpace.v,
                KS.vSpace.weights,
                dt,
                1.0,
                ctr[i-1].sh0,
                ctr[i-1].sh1,
                ctr[i-1].sh2,
                ctr[i].sh0,
                ctr[i].sh1,
                ctr[i].sh2,
            )
        end
    elseif mode == :ugks
        @inbounds Threads.@threads for i = idx0:idx1
            flux_ugks!(
                face[i].fw,
                face[i].fh0,
                face[i].fh1,
                face[i].fh2,
                ctr[i-1].w .+ 0.5 .* KS.ps.dx[i-1] .* ctr[i-1].sw,
                ctr[i-1].h0 .+ 0.5 .* KS.ps.dx[i-1] .* ctr[i-1].sh0,
                ctr[i-1].h1 .+ 0.5 .* KS.ps.dx[i-1] .* ctr[i-1].sh1,
                ctr[i-1].h2 .+ 0.5 .* KS.ps.dx[i-1] .* ctr[i-1].sh2,
                ctr[i].w .- 0.5 .* KS.ps.dx[i] .* ctr[i].sw,
                ctr[i].h0 .- 0.5 .* KS.ps.dx[i] .* ctr[i].sh0,
                ctr[i].h1 .- 0.5 .* KS.ps.dx[i] .* ctr[i].sh1,
                ctr[i].h2 .- 0.5 .* KS.ps.dx[i] .* ctr[i].sh2,
                KS.vSpace.u,
                KS.vSpace.v,
                KS.vSpace.weights,
                KS.gas.K,
                KS.gas.γ,
                KS.gas.mi,
                KS.gas.ni,
                KS.gas.me,
                KS.gas.ne,
                KS.gas.Kn[1],
                dt,
                0.5 * KS.ps.dx[i-1],
                0.5 * KS.ps.dx[i],
                1.0,
                ctr[i-1].sh0,
                ctr[i-1].sh1,
                ctr[i-1].sh2,
                ctr[i].sh0,
                ctr[i].sh1,
                ctr[i].sh2,
            )
        end
    end

    if isPlasma
        @inbounds Threads.@threads for i = idx0:idx1
            flux_em!(
                face[i].femL,
                face[i].femR,
                ctr[i-2].E,
                ctr[i-2].B,
                ctr[i-1].E,
                ctr[i-1].B,
                ctr[i].E,
                ctr[i].B,
                ctr[i+1].E,
                ctr[i+1].B,
                ctr[i-1].ϕ,
                ctr[i].ϕ,
                ctr[i-1].ψ,
                ctr[i].ψ,
                KS.ps.dx[i-1],
                KS.ps.dx[i],
                KS.gas.Ap,
                KS.gas.An,
                KS.gas.D,
                KS.gas.sol,
                KS.gas.χ,
                KS.gas.ν,
                dt,
            )
        end
    end

end

function evolve!(
    KS::SolverSet,
    ctr::T1,
    a1face::T2,
    a2face::T2,
    dt;
    mode = symbolize(KS.set.flux)::Symbol,
    bc = symbolize(KS.set.boundary),
) where {T1<:AA{ControlVolume2D,2},T2<:AA{Interface2D,2}}

    nx, ny, dx, dy = begin
        if KS.ps isa CSpace2D
            KS.ps.nr, KS.ps.nθ, KS.ps.dr, KS.ps.darc
        else
            KS.ps.nx, KS.ps.ny, KS.ps.dx, KS.ps.dy
        end
    end

    if firstindex(KS.pSpace.x[:, 1]) < 1
        idx0 = 1
        idx1 = nx + 1
    else
        idx0 = 2
        idx1 = nx
    end
    if firstindex(KS.pSpace.y[1, :]) < 1
        idy0 = 1
        idy1 = ny + 1
    else
        idy0 = 2
        idy1 = ny
    end

    if mode == :hll

        # x direction
        @inbounds Threads.@threads for j = 1:ny
            for i = idx0:idx1
                flux_hll!(
                    a1face[i, j].fw,
                    local_frame(
                        ctr[i-1, j].w .+ 0.5 .* dx[i-1, j] .* ctr[i-1, j].sw[:, 1],
                        a1face[i, j].n[1],
                        a1face[i, j].n[2],
                    ),
                    local_frame(
                        ctr[i, j].w .- 0.5 .* dx[i, j] .* ctr[i, j].sw[:, 1],
                        a1face[i, j].n[1],
                        a1face[i, j].n[2],
                    ),
                    KS.gas.γ,
                    dt,
                )
                a1face[i, j].fw .=
                    global_frame(a1face[i, j].fw, a1face[i, j].n[1], a1face[i, j].n[2]) .*
                    a1face[i, j].len
            end
        end

        # y direction
        @inbounds Threads.@threads for j = idy0:idy1
            for i = 1:nx
                flux_hll!(
                    a2face[i, j].fw,
                    local_frame(
                        ctr[i, j-1].w .+ 0.5 .* dy[i, j-1] .* ctr[i, j-1].sw[:, 2],
                        a2face[i, j].n[1],
                        a2face[i, j].n[2],
                    ),
                    local_frame(
                        ctr[i, j].w .- 0.5 .* dy[i, j] .* ctr[i, j].sw[:, 2],
                        a2face[i, j].n[1],
                        a2face[i, j].n[2],
                    ),
                    KS.gas.γ,
                    dt,
                )
                a2face[i, j].fw .=
                    global_frame(a2face[i, j].fw, a2face[i, j].n[1], a2face[i, j].n[2]) .*
                    a2face[i, j].len
            end
        end

    end

end

function evolve!(
    KS::SolverSet,
    ctr::T1,
    a1face::T2,
    a2face::T2,
    dt;
    mode = symbolize(KS.set.flux)::Symbol,
    bc = symbolize(KS.set.boundary),
) where {T1<:AA{ControlVolume2D1F,2},T2<:AA{Interface2D1F,2}}

    nx, ny, dx, dy = begin
        if KS.ps isa CSpace2D
            KS.ps.nr, KS.ps.nθ, KS.ps.dr, KS.ps.darc
        else
            KS.ps.nx, KS.ps.ny, KS.ps.dx, KS.ps.dy
        end
    end

    if firstindex(KS.pSpace.x[:, 1]) < 1
        idx0 = 1
        idx1 = nx + 1
    else
        idx0 = 2
        idx1 = nx
    end
    if firstindex(KS.pSpace.y[1, :]) < 1
        idy0 = 1
        idy1 = ny + 1
    else
        idy0 = 2
        idy1 = ny
    end

    if mode == :kfvs

        # x direction
        @inbounds Threads.@threads for j = 1:ny
            for i = idx0:idx1
                vn = KS.vSpace.u .* a1face[i, j].n[1] .+ KS.vSpace.v .* a1face[i, j].n[2]
                vt = KS.vSpace.v .* a1face[i, j].n[1] .- KS.vSpace.u .* a1face[i, j].n[2]

                flux_kfvs!(
                    a1face[i, j].fw,
                    a1face[i, j].ff,
                    ctr[i-1, j].f .+ 0.5 .* dx[i-1, j] .* ctr[i-1, j].sf[:, :, 1],
                    ctr[i, j].f .- 0.5 .* dx[i, j] .* ctr[i, j].sf[:, :, 1],
                    vn,
                    vt,
                    KS.vSpace.weights,
                    dt,
                    a1face[i, j].len,
                    ctr[i-1, j].sf[:, :, 1],
                    ctr[i, j].sf[:, :, 1],
                )
                a1face[i, j].fw .=
                    global_frame(a1face[i, j].fw, a1face[i, j].n[1], a1face[i, j].n[2])
            end
        end

        # y direction
        @inbounds Threads.@threads for j = idy0:idy1
            for i = 1:nx
                vn = KS.vSpace.u .* a2face[i, j].n[1] .+ KS.vSpace.v .* a2face[i, j].n[2]
                vt = KS.vSpace.v .* a2face[i, j].n[1] .- KS.vSpace.u .* a2face[i, j].n[2]

                flux_kfvs!(
                    a2face[i, j].fw,
                    a2face[i, j].ff,
                    ctr[i, j-1].f .+ 0.5 .* dy[i, j-1] .* ctr[i, j-1].sf[:, :, 2],
                    ctr[i, j].f .- 0.5 .* dy[i, j] .* ctr[i, j].sf[:, :, 2],
                    vn,
                    vt,
                    KS.vSpace.weights,
                    dt,
                    a2face[i, j].len,
                    ctr[i, j-1].sf[:, :, 2],
                    ctr[i, j].sf[:, :, 2],
                )
                a2face[i, j].fw .=
                    global_frame(a2face[i, j].fw, a2face[i, j].n[1], a2face[i, j].n[2])
            end
        end

    elseif mode == :kcu

        # x direction
        @inbounds Threads.@threads for j = 1:ny
            for i = idx0:idx1
                vn = KS.vSpace.u .* a1face[i, j].n[1] .+ KS.vSpace.v .* a1face[i, j].n[2]
                vt = KS.vSpace.v .* a1face[i, j].n[1] .- KS.vSpace.u .* a1face[i, j].n[2]

                flux_kcu!(
                    a1face[i, j].fw,
                    a1face[i, j].ff,
                    local_frame(
                        ctr[i-1, j].w .+ 0.5 .* dx[i-1, j] .* ctr[i-1, j].sw[:, 1],
                        a1face[i, j].n[1],
                        a1face[i, j].n[2],
                    ),
                    ctr[i-1, j].f .+ 0.5 .* dx[i-1, j] .* ctr[i-1, j].sf[:, :, 1],
                    local_frame(
                        ctr[i, j].w .- 0.5 .* dx[i, j] .* ctr[i, j].sw[:, 1],
                        a1face[i, j].n[1],
                        a1face[i, j].n[2],
                    ),
                    ctr[i, j].f .- 0.5 .* dx[i, j] .* ctr[i, j].sf[:, :, 1],
                    vn,
                    vt,
                    KS.vSpace.weights,
                    KS.gas.K,
                    KS.gas.γ,
                    KS.gas.μᵣ,
                    KS.gas.ω,
                    KS.gas.Pr,
                    dt,
                    a1face[i, j].len,
                )
                a1face[i, j].fw .=
                    global_frame(a1face[i, j].fw, a1face[i, j].n[1], a1face[i, j].n[2])
            end
        end

        # y direction
        @inbounds Threads.@threads for j = idy0:idy1
            for i = 1:nx
                vn = KS.vSpace.u .* a2face[i, j].n[1] .+ KS.vSpace.v .* a2face[i, j].n[2]
                vt = KS.vSpace.v .* a2face[i, j].n[1] .- KS.vSpace.u .* a2face[i, j].n[2]

                flux_kcu!(
                    a2face[i, j].fw,
                    a2face[i, j].ff,
                    local_frame(
                        ctr[i, j-1].w .+ 0.5 .* dy[i, j-1] .* ctr[i, j-1].sw[:, 2],
                        a2face[i, j].n[1],
                        a2face[i, j].n[2],
                    ),
                    ctr[i, j-1].f .+ 0.5 .* dy[i, j-1] .* ctr[i, j-1].sf[:, :, 2],
                    local_frame(
                        ctr[i, j].w .- 0.5 .* dy[i, j] .* ctr[i, j].sw[:, 2],
                        a2face[i, j].n[1],
                        a2face[i, j].n[2],
                    ),
                    ctr[i, j].f .- 0.5 .* dy[i, j] .* ctr[i, j].sf[:, :, 2],
                    vn,
                    vt,
                    KS.vSpace.weights,
                    KS.gas.K,
                    KS.gas.γ,
                    KS.gas.μᵣ,
                    KS.gas.ω,
                    KS.gas.Pr,
                    dt,
                    a2face[i, j].len,
                )
                a2face[i, j].fw .=
                    global_frame(a2face[i, j].fw, a2face[i, j].n[1], a2face[i, j].n[2])
            end
        end

    end

    bcs = ifelse(bc isa Symbol, [bc, bc, bc, bc], bc)
    if bcs[1] == :maxwell
        @inbounds Threads.@threads for j = 1:ny
            vn = KS.vSpace.u .* a1face[1, j].n[1] .+ KS.vSpace.v .* a1face[1, j].n[2]
            vt = KS.vSpace.v .* a1face[1, j].n[1] .- KS.vSpace.u .* a1face[1, j].n[2]
            xc = (KS.ps.vertices[1, j, 1, 1] + KS.ps.vertices[1, j, 4, 1]) / 2
            yc = (KS.ps.vertices[1, j, 1, 2] + KS.ps.vertices[1, j, 4, 2]) / 2
            bcL = local_frame(KS.ib.bc(xc, yc), a1face[1, j].n[1], a1face[1, j].n[2])

            flux_boundary_maxwell!(
                a1face[1, j].fw,
                a1face[1, j].ff,
                bcL, # left
                ctr[1, j].f,
                vn,
                vt,
                KS.vSpace.weights,
                dt,
                KS.ps.dy[1, j],
                1,
            )
            a1face[1, j].fw .=
                global_frame(a1face[1, j].fw, a1face[1, j].n[1], a1face[1, j].n[2])
        end
    end
    if bcs[2] == :maxwell
        @inbounds Threads.@threads for j = 1:ny
            vn = KS.vSpace.u .* a1face[nx+1, j].n[1] .+ KS.vSpace.v .* a1face[nx+1, j].n[2]
            vt = KS.vSpace.v .* a1face[nx+1, j].n[1] .- KS.vSpace.u .* a1face[nx+1, j].n[2]
            xc = (KS.ps.vertices[nx, j, 2, 1] + KS.ps.vertices[nx, j, 3, 1]) / 2
            yc = (KS.ps.vertices[nx, j, 2, 2] + KS.ps.vertices[nx, j, 3, 2]) / 2
            bcR = local_frame(KS.ib.bc(xc, yc), a1face[nx+1, j].n[1], a1face[nx+1, j].n[2])

            flux_boundary_maxwell!(
                a1face[nx+1, j].fw,
                a1face[nx+1, j].ff,
                bcR, # right
                ctr[nx, j].f,
                vn,
                vt,
                KS.vSpace.weights,
                dt,
                KS.ps.dy[nx, j],
                -1,
            )
            a1face[nx+1, j].fw .=
                global_frame(a1face[nx+1, j].fw, a1face[nx+1, j].n[1], a1face[nx+1, j].n[2])
        end
    end
    if bcs[3] == :maxwell
        @inbounds Threads.@threads for i = 1:nx
            vn = KS.vSpace.u .* a2face[i, 1].n[1] .+ KS.vSpace.v .* a2face[i, 1].n[2]
            vt = KS.vSpace.v .* a2face[i, 1].n[1] .- KS.vSpace.u .* a2face[i, 1].n[2]
            xc = (KS.ps.vertices[i, 1, 1, 1] + KS.ps.vertices[i, 1, 2, 1]) / 2
            yc = (KS.ps.vertices[i, 1, 1, 2] + KS.ps.vertices[i, 1, 2, 2]) / 2
            bcD = local_frame(KS.ib.bc(xc, yc), a2face[i, 1].n[1], a2face[i, 1].n[2])

            flux_boundary_maxwell!(
                a2face[i, 1].fw,
                a2face[i, 1].ff,
                bcD,
                ctr[i, 1].f,
                vn,
                vt,
                KS.vSpace.weights,
                dt,
                dx[i, 1],
                1,
            )
            a2face[i, 1].fw .=
                global_frame(a2face[i, 1].fw, a2face[i, 1].n[1], a2face[i, 1].n[2])
        end
    end
    if bcs[4] == :maxwell
        @inbounds Threads.@threads for i = 1:nx
            vn = KS.vSpace.u .* a2face[i, ny+1].n[1] .+ KS.vSpace.v .* a2face[i, ny+1].n[2]
            vt = KS.vSpace.v .* a2face[i, ny+1].n[1] .- KS.vSpace.u .* a2face[i, ny+1].n[2]
            xc = (KS.ps.vertices[i, ny, 3, 1] + KS.ps.vertices[i, ny, 4, 1]) / 2
            yc = (KS.ps.vertices[i, ny, 3, 2] + KS.ps.vertices[i, ny, 4, 2]) / 2
            bcU = local_frame(KS.ib.bc(xc, yc), a2face[i, ny+1].n[1], a2face[i, ny+1].n[2])

            flux_boundary_maxwell!(
                a2face[i, ny+1].fw,
                a2face[i, ny+1].ff,
                bcU,
                ctr[i, ny].f,
                vn,
                vt,
                KS.vSpace.weights,
                dt,
                KS.ps.dx[i, ny],
                -1,
            )
            a2face[i, ny+1].fw .=
                global_frame(a2face[i, ny+1].fw, a2face[i, ny+1].n[1], a2face[i, ny+1].n[2])
        end
    end

end

function evolve!(
    KS::SolverSet,
    ctr::T1,
    a1face::T2,
    a2face::T2,
    dt;
    mode = symbolize(KS.set.flux)::Symbol,
    bc = symbolize(KS.set.boundary),
) where {T1<:AA{ControlVolume2D2F,2},T2<:AA{Interface2D2F,2}}

    nx, ny, dx, dy = begin
        if KS.ps isa CSpace2D
            KS.ps.nr, KS.ps.nθ, KS.ps.dr, KS.ps.darc
        else
            KS.ps.nx, KS.ps.ny, KS.ps.dx, KS.ps.dy
        end
    end

    if firstindex(KS.pSpace.x[:, 1]) < 1
        idx0 = 1
        idx1 = nx + 1
    else
        idx0 = 2
        idx1 = nx
    end
    if firstindex(KS.pSpace.y[1, :]) < 1
        idy0 = 1
        idy1 = ny + 1
    else
        idy0 = 2
        idy1 = ny
    end

    if mode == :kfvs

        # x direction
        @inbounds Threads.@threads for j = 1:ny
            for i = idx0:idx1
                vn = KS.vSpace.u .* a1face[i, j].n[1] .+ KS.vSpace.v .* a1face[i, j].n[2]
                vt = KS.vSpace.v .* a1face[i, j].n[1] .- KS.vSpace.u .* a1face[i, j].n[2]

                flux_kfvs!(
                    a1face[i, j].fw,
                    a1face[i, j].fh,
                    a1face[i, j].fb,
                    ctr[i-1, j].h .+ 0.5 .* dx[i-1, j] .* ctr[i-1, j].sh[:, :, 1],
                    ctr[i-1, j].b .+ 0.5 .* dx[i-1, j] .* ctr[i-1, j].sb[:, :, 1],
                    ctr[i, j].h .- 0.5 .* dx[i, j] .* ctr[i, j].sh[:, :, 1],
                    ctr[i, j].b .- 0.5 .* dx[i, j] .* ctr[i, j].sb[:, :, 1],
                    vn,
                    vt,
                    KS.vSpace.weights,
                    dt,
                    a1face[i, j].len,
                    ctr[i-1, j].sh[:, :, 1],
                    ctr[i-1, j].sb[:, :, 1],
                    ctr[i, j].sh[:, :, 1],
                    ctr[i, j].sb[:, :, 1],
                )
                a1face[i, j].fw .=
                    global_frame(a1face[i, j].fw, a1face[i, j].n[1], a1face[i, j].n[2])
            end
        end

        # y direction
        @inbounds Threads.@threads for j = idy0:idy1
            for i = 1:nx
                vn = KS.vSpace.u .* a2face[i, j].n[1] .+ KS.vSpace.v .* a2face[i, j].n[2]
                vt = KS.vSpace.v .* a2face[i, j].n[1] .- KS.vSpace.u .* a2face[i, j].n[2]

                flux_kfvs!(
                    a2face[i, j].fw,
                    a2face[i, j].fh,
                    a2face[i, j].fb,
                    ctr[i, j-1].h .+ 0.5 .* dy[i, j-1] .* ctr[i, j-1].sh[:, :, 2],
                    ctr[i, j-1].b .+ 0.5 .* dy[i, j-1] .* ctr[i, j-1].sb[:, :, 2],
                    ctr[i, j].h .- 0.5 .* dy[i, j] .* ctr[i, j].sh[:, :, 2],
                    ctr[i, j].b .- 0.5 .* dy[i, j] .* ctr[i, j].sb[:, :, 2],
                    vn,
                    vt,
                    KS.vSpace.weights,
                    dt,
                    a2face[i, j].len,
                    ctr[i, j-1].sh[:, :, 2],
                    ctr[i, j-1].sb[:, :, 2],
                    ctr[i, j].sh[:, :, 2],
                    ctr[i, j].sb[:, :, 2],
                )
                a2face[i, j].fw .=
                    global_frame(a2face[i, j].fw, a2face[i, j].n[1], a2face[i, j].n[2])
            end
        end

    elseif mode == :kcu

        # x direction
        @inbounds Threads.@threads for j = 1:ny
            for i = idx0:idx1
                vn = KS.vSpace.u .* a1face[i, j].n[1] .+ KS.vSpace.v .* a1face[i, j].n[2]
                vt = KS.vSpace.v .* a1face[i, j].n[1] .- KS.vSpace.u .* a1face[i, j].n[2]

                flux_kcu!(
                    a1face[i, j].fw,
                    a1face[i, j].fh,
                    a1face[i, j].fb,
                    local_frame(
                        ctr[i-1, j].w .+ 0.5 .* dx[i-1, j] .* ctr[i-1, j].sw[:, 1],
                        a1face[i, j].n[1],
                        a1face[i, j].n[2],
                    ),
                    ctr[i-1, j].h .+ 0.5 .* dx[i-1, j] .* ctr[i-1, j].sh[:, :, 1],
                    ctr[i-1, j].b .+ 0.5 .* dx[i-1, j] .* ctr[i-1, j].sb[:, :, 1],
                    local_frame(
                        ctr[i, j].w .- 0.5 .* dx[i, j] .* ctr[i, j].sw[:, 1],
                        a1face[i, j].n[1],
                        a1face[i, j].n[2],
                    ),
                    ctr[i, j].h .- 0.5 .* dx[i, j] .* ctr[i, j].sh[:, :, 1],
                    ctr[i, j].b .- 0.5 .* dx[i, j] .* ctr[i, j].sb[:, :, 1],
                    vn,
                    vt,
                    KS.vSpace.weights,
                    KS.gas.K,
                    KS.gas.γ,
                    KS.gas.μᵣ,
                    KS.gas.ω,
                    KS.gas.Pr,
                    dt,
                    a1face[i, j].len,
                )
                a1face[i, j].fw .=
                    global_frame(a1face[i, j].fw, a1face[i, j].n[1], a1face[i, j].n[2])
            end
        end

        # y direction
        @inbounds Threads.@threads for j = idy0:idy1
            for i = 1:nx
                vn = KS.vSpace.u .* a2face[i, j].n[1] .+ KS.vSpace.v .* a2face[i, j].n[2]
                vt = KS.vSpace.v .* a2face[i, j].n[1] .- KS.vSpace.u .* a2face[i, j].n[2]

                flux_kcu!(
                    a2face[i, j].fw,
                    a2face[i, j].fh,
                    a2face[i, j].fb,
                    local_frame(
                        ctr[i, j-1].w .+ 0.5 .* dy[i, j-1] .* ctr[i, j-1].sw[:, 2],
                        a2face[i, j].n[1],
                        a2face[i, j].n[2],
                    ),
                    ctr[i, j-1].h .+ 0.5 .* dy[i, j-1] .* ctr[i, j-1].sh[:, :, 2],
                    ctr[i, j-1].b .+ 0.5 .* dy[i, j-1] .* ctr[i, j-1].sb[:, :, 2],
                    local_frame(
                        ctr[i, j].w .- 0.5 .* dy[i, j] .* ctr[i, j].sw[:, 2],
                        a2face[i, j].n[1],
                        a2face[i, j].n[2],
                    ),
                    ctr[i, j].h .- 0.5 .* dy[i, j] .* ctr[i, j].sh[:, :, 2],
                    ctr[i, j].b .- 0.5 .* dy[i, j] .* ctr[i, j].sb[:, :, 2],
                    vn,
                    vt,
                    KS.vSpace.weights,
                    KS.gas.K,
                    KS.gas.γ,
                    KS.gas.μᵣ,
                    KS.gas.ω,
                    KS.gas.Pr,
                    dt,
                    a2face[i, j].len,
                )
                a2face[i, j].fw .=
                    global_frame(a2face[i, j].fw, a2face[i, j].n[1], a2face[i, j].n[2])
            end
        end

    end

    evolve_boundary!(KS, ctr, a1face, a2face, dt, mode, bc)

    return nothing

end

function evolve!(
    KS::SolverSet,
    ctr::T1,
    face::T2,
    dt;
    mode = symbolize(KS.set.flux)::Symbol,
    bc = symbolize(KS.set.boundary),
) where {T1<:AbstractVector{ControlVolumeUS},T2<:AbstractVector{Interface2D}}

    if mode == :hll

        @inbounds Threads.@threads for i in eachindex(face)
            vn = KS.vSpace.u .* face[i].n[1] .+ KS.vSpace.v .* face[i].n[2]
            vt = KS.vSpace.v .* face[i].n[1] .- KS.vSpace.u .* face[i].n[2]

            if !(-1 in KS.ps.faceCells[i, :]) # exclude boundary face
                flux_hll!(
                    face[i].fw,
                    ctr[KS.ps.faceCells[i, 1]].w .+
                    ctr[KS.ps.faceCells[i, 1]].sw[:, 1] .*
                    (KS.ps.faceCenter[i, 1] - KS.ps.cellCenter[KS.ps.faceCells[i, 1], 1]) .+
                    ctr[KS.ps.faceCells[i, 1]].sw[:, 2] .*
                    (KS.ps.faceCenter[i, 2] - KS.ps.cellCenter[KS.ps.faceCells[i, 1], 2]),
                    ctr[KS.ps.faceCells[i, 2]].w .+
                    ctr[KS.ps.faceCells[i, 2]].sw[:, 1] .*
                    (KS.ps.faceCenter[i, 1] - KS.ps.cellCenter[KS.ps.faceCells[i, 2], 1]) .+
                    ctr[KS.ps.faceCells[i, 2]].sw[:, 2] .*
                    (KS.ps.faceCenter[i, 2] - KS.ps.cellCenter[KS.ps.faceCells[i, 2], 2]),
                    KS.gas.γ,
                    dt,
                )
                face[i].fw .=
                    KitBase.global_frame(face[i].fw, face[i].n[1], face[i].n[2]) .*
                    face[i].len
            else
                idx = ifelse(KS.ps.faceCells[i, 1] != -1, 1, 2)

                if KS.ps.cellType[KS.ps.faceCells[i, idx]] == 2
                    prim = KitBase.local_frame(
                        ctr[KS.ps.faceCells[i, idx]].prim,
                        face[i].n[1],
                        face[i].n[2],
                    )
                    #bc = copy(KS.ib.bcR)
                    bc0 = KS.ib.bc(KS.ps.faceCenter[i, 1], KS.ps.faceCenter[i, 2])
                    bc = zero(bc)
                    bc[2] = -prim[2]
                    bc[3] = -prim[3]
                    bc[4] = 2.0 * bc0[4] - prim[4]
                    bc[1] = (2.0 * bc0[4] - prim[4]) / prim[4] * prim[1]
                    bc_w = KitBase.prim_conserve(bc, KS.gas.γ)

                    if idx == 1
                        wL = bc_w
                        wR = KitBase.local_frame(
                            ctr[KS.ps.faceCells[i, idx]].w,
                            face[i].n[1],
                            face[i].n[2],
                        )
                    else
                        wL = KitBase.local_frame(
                            ctr[KS.ps.faceCells[i, idx]].w,
                            face[i].n[1],
                            face[i].n[2],
                        )
                        wR = bc_w
                    end

                    KitBase.flux_hll!(face[i].fw, wL, wR, KS.gas.γ, dt)

                    face[i].fw .=
                        KitBase.global_frame(face[i].fw, face[i].n[1], face[i].n[2]) .*
                        face[i].len
                end
            end
        end

    end

end

function evolve!(
    KS::SolverSet,
    ctr::T1,
    face::T2,
    dt;
    mode = symbolize(KS.set.flux)::Symbol,
    bc = symbolize(KS.set.boundary),
) where {T1<:AbstractVector{ControlVolumeUS1F},T2<:AbstractVector{Interface2D1F}}

    if mode == :kfvs

        @inbounds Threads.@threads for i in eachindex(face)
            vn = KS.vSpace.u .* face[i].n[1] .+ KS.vSpace.v .* face[i].n[2]
            vt = KS.vSpace.v .* face[i].n[1] .- KS.vSpace.u .* face[i].n[2]

            if !(-1 in KS.ps.faceCells[i, :]) # exclude boundary face
                flux_kfvs!(
                    face[i].fw,
                    face[i].ff,
                    ctr[KS.ps.faceCells[i, 1]].f .+
                    ctr[KS.ps.faceCells[i, 1]].sf[:, :, 1] .*
                    (KS.ps.faceCenter[i, 1] - KS.ps.cellCenter[KS.ps.faceCells[i, 1], 1]) .+
                    ctr[KS.ps.faceCells[i, 1]].sf[:, :, 2] .*
                    (KS.ps.faceCenter[i, 2] - KS.ps.cellCenter[KS.ps.faceCells[i, 1], 2]),
                    ctr[KS.ps.faceCells[i, 2]].f .+
                    ctr[KS.ps.faceCells[i, 2]].sf[:, :, 1] .*
                    (KS.ps.faceCenter[i, 1] - KS.ps.cellCenter[KS.ps.faceCells[i, 2], 1]) .+
                    ctr[KS.ps.faceCells[i, 2]].sf[:, :, 2] .*
                    (KS.ps.faceCenter[i, 2] - KS.ps.cellCenter[KS.ps.faceCells[i, 2], 2]),
                    vn,
                    vt,
                    KS.vSpace.weights,
                    dt,
                    face[i].len,
                )
                face[i].fw .= KitBase.global_frame(face[i].fw, face[i].n[1], face[i].n[2])
            else
                idx = ifelse(KS.ps.faceCells[i, 1] != -1, 1, 2)

                if KS.ps.cellType[KS.ps.faceCells[i, idx]] == 2
                    bc = KitBase.local_frame(
                        KS.ib.bc(KS.ps.faceCenter[i, 1], KS.ps.faceCenter[i, 2]),
                        face[i].n[1],
                        face[i].n[2],
                    )

                    KitBase.flux_boundary_maxwell!(
                        face[i].fw,
                        face[i].ff,
                        bc,
                        ctr[KS.ps.faceCells[i, idx]].f .+
                        ctr[KS.ps.faceCells[i, idx]].sf[:, :, 1] .* (
                            KS.ps.faceCenter[i, 1] -
                            KS.ps.cellCenter[KS.ps.faceCells[i, idx], 1]
                        ) .+
                        ctr[KS.ps.faceCells[i, idx]].sf[:, :, 2] .* (
                            KS.ps.faceCenter[i, 2] -
                            KS.ps.cellCenter[KS.ps.faceCells[i, idx], 2]
                        ),
                        vn,
                        vt,
                        KS.vSpace.weights,
                        KS.gas.K,
                        dt,
                        face[i].len,
                    )

                    face[i].fw .=
                        KitBase.global_frame(face[i].fw, face[i].n[1], face[i].n[2])
                end
            end
        end

    end

end

function evolve!(
    KS::SolverSet,
    ctr::T1,
    face::T2,
    dt;
    mode = symbolize(KS.set.flux)::Symbol,
    bc = symbolize(KS.set.boundary),
) where {T1<:AbstractVector{ControlVolumeUS2F},T2<:AbstractVector{Interface2D2F}}

    if mode == :kfvs

        @inbounds Threads.@threads for i in eachindex(face)
            vn = KS.vSpace.u .* face[i].n[1] .+ KS.vSpace.v .* face[i].n[2]
            vt = KS.vSpace.v .* face[i].n[1] .- KS.vSpace.u .* face[i].n[2]

            if !(-1 in KS.ps.faceCells[i, :]) # exclude boundary face
                flux_kfvs!(
                    face[i].fw,
                    face[i].fh,
                    face[i].fb,
                    ctr[KS.ps.faceCells[i, 1]].h .+
                    ctr[KS.ps.faceCells[i, 1]].sh[:, :, 1] .*
                    (KS.ps.faceCenter[i, 1] - KS.ps.cellCenter[KS.ps.faceCells[i, 1], 1]) .+
                    ctr[KS.ps.faceCells[i, 1]].sh[:, :, 2] .*
                    (KS.ps.faceCenter[i, 2] - KS.ps.cellCenter[KS.ps.faceCells[i, 1], 2]),
                    ctr[KS.ps.faceCells[i, 1]].b .+
                    ctr[KS.ps.faceCells[i, 1]].sb[:, :, 1] .*
                    (KS.ps.faceCenter[i, 1] - KS.ps.cellCenter[KS.ps.faceCells[i, 1], 1]) .+
                    ctr[KS.ps.faceCells[i, 1]].sb[:, :, 2] .*
                    (KS.ps.faceCenter[i, 2] - KS.ps.cellCenter[KS.ps.faceCells[i, 1], 2]),
                    ctr[KS.ps.faceCells[i, 2]].h .+
                    ctr[KS.ps.faceCells[i, 2]].sh[:, :, 1] .*
                    (KS.ps.faceCenter[i, 1] - KS.ps.cellCenter[KS.ps.faceCells[i, 2], 1]) .+
                    ctr[KS.ps.faceCells[i, 2]].sh[:, :, 2] .*
                    (KS.ps.faceCenter[i, 2] - KS.ps.cellCenter[KS.ps.faceCells[i, 2], 2]),
                    ctr[KS.ps.faceCells[i, 2]].b .+
                    ctr[KS.ps.faceCells[i, 2]].sb[:, :, 1] .*
                    (KS.ps.faceCenter[i, 1] - KS.ps.cellCenter[KS.ps.faceCells[i, 2], 1]) .+
                    ctr[KS.ps.faceCells[i, 2]].sb[:, :, 2] .*
                    (KS.ps.faceCenter[i, 2] - KS.ps.cellCenter[KS.ps.faceCells[i, 2], 2]),
                    vn,
                    vt,
                    KS.vSpace.weights,
                    dt,
                    face[i].len,
                )
                face[i].fw .= KitBase.global_frame(face[i].fw, face[i].n[1], face[i].n[2])
            else
                idx = ifelse(KS.ps.faceCells[i, 1] != -1, 1, 2)

                if KS.ps.cellType[KS.ps.faceCells[i, idx]] == 2
                    bc = KitBase.local_frame(
                        KS.ib.bc(KS.ps.faceCenter[i, 1], KS.ps.faceCenter[i, 2]),
                        face[i].n[1],
                        face[i].n[2],
                    )

                    KitBase.flux_boundary_maxwell!(
                        face[i].fw,
                        face[i].fh,
                        face[i].fb,
                        bc,
                        ctr[KS.ps.faceCells[i, idx]].h .+
                        ctr[KS.ps.faceCells[i, idx]].sh[:, :, 1] .* (
                            KS.ps.faceCenter[i, 1] -
                            KS.ps.cellCenter[KS.ps.faceCells[i, idx], 1]
                        ) .+
                        ctr[KS.ps.faceCells[i, idx]].sh[:, :, 2] .* (
                            KS.ps.faceCenter[i, 2] -
                            KS.ps.cellCenter[KS.ps.faceCells[i, idx], 2]
                        ),
                        ctr[KS.ps.faceCells[i, idx]].b .+
                        ctr[KS.ps.faceCells[i, idx]].sb[:, :, 1] .* (
                            KS.ps.faceCenter[i, 1] -
                            KS.ps.cellCenter[KS.ps.faceCells[i, idx], 1]
                        ) .+
                        ctr[KS.ps.faceCells[i, idx]].sb[:, :, 2] .* (
                            KS.ps.faceCenter[i, 2] -
                            KS.ps.cellCenter[KS.ps.faceCells[i, idx], 2]
                        ),
                        vn,
                        vt,
                        KS.vSpace.weights,
                        KS.gas.K,
                        dt,
                        face[i].len,
                    )

                    face[i].fw .=
                        KitBase.global_frame(face[i].fw, face[i].n[1], face[i].n[2])
                end
            end
        end

    end

end
