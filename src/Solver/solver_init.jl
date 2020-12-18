# ============================================================
# Initializer of Simulation
# ============================================================


"""
Initialize solver from input file

"""
function initialize(configfilename::T) where {T<:AbstractString}

    println("==============================================================")
    println("Kinetic.jl")
    println("A Lightweight Toolbox for Kinetic Modeling and Simulation")
    println("==============================================================")
    println("")
    println("reading configurations from $(configfilename)")
    println("")
    println("initializing solver: ")

    if configfilename[end-2:end] == "txt"

        ks = SolverSet(configfilename)
        ctr, face = init_fvm(ks)

        return ks, ctr, face, 0.0

    elseif configfilename[end-3:end] == "jld2"

        _1, _2, _3 = @load configfilename KS ctr t
        ks, ctr, simTime = eval(_1), eval(_2), eval(_3)

        face = init_fvm(ks)[2]

        return ks, ctr, face, simTime

    end

end


"""
Initialize finite volume method

"""
function init_fvm(KS::T) where {T<:AbstractSolverSet}

    if KS.set.space[1:2] == "1d"

        if KS.set.space[3:4] == "1f"

            ctr = OffsetArray{ControlVolume1D1F}(undef, eachindex(KS.pSpace.x))
            face = Array{Interface1D1F}(undef, KS.pSpace.nx + 1)

            for i in eachindex(ctr)
                if i <= KS.pSpace.nx ÷ 2
                    ctr[i] = ControlVolume1D1F(
                        KS.pSpace.x[i],
                        KS.pSpace.dx[i],
                        KS.ib.wL,
                        KS.ib.primL,
                        KS.ib.fL,
                    )
                else
                    ctr[i] = ControlVolume1D1F(
                        KS.pSpace.x[i],
                        KS.pSpace.dx[i],
                        KS.ib.wR,
                        KS.ib.primR,
                        KS.ib.fR,
                    )
                end
            end

            for i = 1:KS.pSpace.nx+1
                face[i] = Interface1D1F(KS.ib.wL, KS.ib.fL)
            end

        elseif KS.set.space[3:4] == "2f"

            ctr = OffsetArray{ControlVolume1D2F}(undef, eachindex(KS.pSpace.x))
            face = Array{Interface1D2F}(undef, KS.pSpace.nx + 1)

            for i in eachindex(ctr)
                if i <= KS.pSpace.nx ÷ 2
                    ctr[i] = ControlVolume1D2F(
                        KS.pSpace.x[i],
                        KS.pSpace.dx[i],
                        KS.ib.wL,
                        KS.ib.primL,
                        KS.ib.hL,
                        KS.ib.bL,
                    )
                else
                    ctr[i] = ControlVolume1D2F(
                        KS.pSpace.x[i],
                        KS.pSpace.dx[i],
                        KS.ib.wR,
                        KS.ib.primR,
                        KS.ib.hR,
                        KS.ib.bR,
                    )
                end
            end

            for i = 1:KS.pSpace.nx+1
                face[i] = Interface1D2F(KS.ib.wL, KS.ib.hL)
            end

        elseif KS.set.space[3:4] == "3f"

            ctr = OffsetArray{ControlVolume1D3F}(undef, eachindex(KS.pSpace.x))
            face = Array{Interface1D3F}(undef, KS.pSpace.nx + 1)

            for i in eachindex(ctr)
                if i <= KS.pSpace.nx ÷ 2
                    ctr[i] = ControlVolume1D3F(
                        KS.pSpace.x[i],
                        KS.pSpace.dx[i],
                        KS.ib.wL,
                        KS.ib.primL,
                        KS.ib.h0L,
                        KS.ib.h1L,
                        KS.ib.h2L,
                        KS.ib.EL,
                        KS.ib.BL,
                        KS.ib.lorenzL,
                    )
                else
                    ctr[i] = ControlVolume1D3F(
                        KS.pSpace.x[i],
                        KS.pSpace.dx[i],
                        KS.ib.wR,
                        KS.ib.primR,
                        KS.ib.h0R,
                        KS.ib.h1R,
                        KS.ib.h2R,
                        KS.ib.ER,
                        KS.ib.BR,
                        KS.ib.lorenzR,
                    )
                end
            end

            for i = 1:KS.pSpace.nx+1
                face[i] = Interface1D3F(KS.ib.wL, KS.ib.h0L, KS.ib.EL)
            end

        elseif KS.set.space[3:4] == "4f"

            ctr = OffsetArray{ControlVolume1D4F}(undef, eachindex(KS.pSpace.x))
            face = Array{Interface1D4F}(undef, KS.pSpace.nx + 1)

            for i in eachindex(ctr)
                if i <= KS.pSpace.nx ÷ 2
                    ctr[i] = ControlVolume1D4F(
                        KS.pSpace.x[i],
                        KS.pSpace.dx[i],
                        KS.ib.wL,
                        KS.ib.primL,
                        KS.ib.h0L,
                        KS.ib.h1L,
                        KS.ib.h2L,
                        KS.ib.h3L,
                        KS.ib.EL,
                        KS.ib.BL,
                        KS.ib.lorenzL,
                    )
                else
                    ctr[i] = ControlVolume1D4F(
                        KS.pSpace.x[i],
                        KS.pSpace.dx[i],
                        KS.ib.wR,
                        KS.ib.primR,
                        KS.ib.h0R,
                        KS.ib.h1R,
                        KS.ib.h2R,
                        KS.ib.h3R,
                        KS.ib.ER,
                        KS.ib.BR,
                        KS.ib.lorenzR,
                    )
                end
            end

            for i = 1:KS.pSpace.nx+1
                face[i] = Interface1D4F(KS.ib.wL, KS.ib.h0L, KS.ib.EL)
            end

        end

        return ctr, face

    elseif KS.set.space[1:2] == "2d"

        if KS.set.space[3:4] == "1f"

            ctr = OffsetArray{ControlVolume2D1F}(
                undef,
                axes(KS.pSpace.x, 1),
                axes(KS.pSpace.y, 2),
            )
            a1face = Array{Interface2D1F}(undef, KS.pSpace.nx + 1, KS.pSpace.ny)
            a2face = Array{Interface2D1F}(undef, KS.pSpace.nx, KS.pSpace.ny + 1)

            for j in axes(ctr, 2), i in axes(ctr, 1)
                if i <= KS.pSpace.nx ÷ 2
                    ctr[i, j] = ControlVolume2D1F(
                        KS.pSpace.x[i, j],
                        KS.pSpace.y[i, j],
                        KS.pSpace.dx[i, j],
                        KS.pSpace.dy[i, j],
                        KS.ib.wL,
                        KS.ib.primL,
                        KS.ib.fL,
                    )
                else
                    ctr[i] = ControlVolume2D1F(
                        KS.pSpace.x[i, j],
                        KS.pSpace.y[i, j],
                        KS.pSpace.dx[i, j],
                        KS.pSpace.dy[i, j],
                        KS.ib.wR,
                        KS.ib.primR,
                        KS.ib.fR,
                    )
                end
            end

            for j = 1:KS.pSpace.ny, i = 1:KS.pSpace.nx+1
                a1face[i] = Interface2D1F(KS.pSpace.dy[i, j], 0.0, 1.0, KS.ib.wL, KS.ib.hL)
            end
            for j = 1:KS.pSpace.ny+1, i = 1:KS.pSpace.nx
                a2face[i] = Interface2D1F(KS.pSpace.dx[i, j], 1.0, 0.0, KS.ib.wL, KS.ib.hL)
            end

        elseif KS.set.space[3:4] == "2f"

            ctr = OffsetArray{ControlVolume2D2F}(
                undef,
                axes(KS.pSpace.x, 1),
                axes(KS.pSpace.y, 2),
            )
            a1face = Array{Interface2D2F}(undef, KS.pSpace.nx + 1, KS.pSpace.ny)
            a2face = Array{Interface2D2F}(undef, KS.pSpace.nx, KS.pSpace.ny + 1)

            for j in axes(ctr, 2), i in axes(ctr, 1)
                if i <= KS.pSpace.nx ÷ 2
                    ctr[i, j] = ControlVolume2D2F(
                        KS.pSpace.x[i, j],
                        KS.pSpace.y[i, j],
                        KS.pSpace.dx[i, j],
                        KS.pSpace.dy[i, j],
                        KS.ib.wL,
                        KS.ib.primL,
                        KS.ib.hL,
                        KS.ib.bL,
                    )
                else
                    ctr[i, j] = ControlVolume2D2F(
                        KS.pSpace.x[i, j],
                        KS.pSpace.y[i, j],
                        KS.pSpace.dx[i, j],
                        KS.pSpace.dy[i, j],
                        KS.ib.wR,
                        KS.ib.primR,
                        KS.ib.hR,
                        KS.ib.bR,
                    )
                end
            end

            for j = 1:KS.pSpace.ny, i = 1:KS.pSpace.nx+1
                a1face[i] = Interface2D2F(KS.pSpace.dy[i, j], 0.0, 1.0, KS.ib.wL, KS.ib.hL)
            end
            for j = 1:KS.pSpace.ny+1, i = 1:KS.pSpace.nx
                a2face[i] = Interface2D2F(KS.pSpace.dx[i, j], 1.0, 0.0, KS.ib.wL, KS.ib.hL)
            end

        end

        return ctr, a1face, a2face

    end

end


"""
Initialize particles

"""
function init_ptc!(KS::SolverSet, ctr::T; mode = :soa::Symbol) where {T<:AbstractArray{<:AbstractControlVolume1D,1}}
    if mode == :soa
        init_ptc_soa!(KS, ctr)
    elseif mode == :aos
        init_ptc_aos!(KS, ctr)
    end
end


"""
Array of structs

"""
function init_ptc_aos!(KS::SolverSet, ctr::T) where {T<:AbstractArray{<:AbstractControlVolume1D,1}}

    np = 0
    for i in eachindex(ctr)
        np += Int(round(ctr[i].w[1] * ctr[i].dx / KS.gas.m))
    end

    ptc = Array{Particle1D}(undef, 2 * np)
    for i in eachindex(ptc)
        m = KS.gas.m
        x = 0.0
        v = zeros(3)
        e = 0.0
        idx = -7
        flag = 0
        tc = 0.0

        ptc[i] = Particle1D(m, x, v, e, idx, flag, tc)
    end

    np_tmp = 0
    for i in eachindex(ctr)
        npl = Int(round(ctr[i].w[1] * ctr[i].dx / KS.gas.m))
        for j = 1:npl
            np_tmp += 1
            sample_particle!(ptc[np_tmp], KS, ctr[i], i)
        end
    end

    return ptc

end


"""
Struct of arrays

"""
function init_ptc_soa!(KS::SolverSet, ctr::T) where {T<:AbstractArray{<:AbstractControlVolume1D,1}}

    np = 0
    for i in 1:KS.pSpace.nx
        np += round(ctr[i].prim[1] * ctr[i].dx / KS.gas.m) |> Int
    end
    KS.gas.np = np

    m = zeros(np)
    x = zeros(np)
    v = zeros(np, 3)
    e = zeros(np)
    idx = zeros(Int, np)
    ref = zeros(Int, np) # default
    flag = zeros(Int, np) # default
    tc = zeros(np)

    for i in 1:np
        m[i] = KS.gas.m
        x[i] = KS.pSpace.x0 + rand() * (KS.pSpace.x1 - KS.pSpace.x0)
        idx[i] = find_idx(KS.pSpace.x[1:end], x[i], mode=:uniform)
        v[i, :] .= sample_maxwell(ctr[idx[i]].prim)
        e[i] = 0.5 / ctr[idx[i]].prim[end]
        τ = vhs_collision_time(ctr[idx[i]].prim, KS.gas.μᵣ, KS.gas.ω)
        tc[i] = next_collision_time(τ)
    end

    ptc = Particle(m, x, v, e, idx, ref, flag, tc)

    return ptc

end