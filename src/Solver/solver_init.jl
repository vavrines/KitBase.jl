# ============================================================
# Solver Initializer
# ============================================================

"""
    initialize(configfilename::T) where {T<:AbstractString}
    initialize(config::T) where {T<:AbstractDict}

Initialize solver from input file or dictionary.
This can also be done from a Julia script directly.

"""
function initialize(configfilename::T) where {T<:AbstractString}
    println("==============================================================")
    println("Kinetic.jl")
    println("Portable Kinetic Simulation and Scientific Machine Learning")
    println("==============================================================")
    println("")
    @info "initializing solver"
    println("")

    if configfilename[end-3:end] == "jld2"
        D = load(configfilename)
        ks = D["set"]
        ctr = D["ctr"]
        t = D["t"]

        face = init_fvm(ks, ks.pSpace)[2]

        return ks, ctr, face, t
    else
        ks = SolverSet(configfilename)

        if ks.set.space[1:2] == "1d"
            ctr, face = init_fvm(ks, ks.pSpace)
            return ks, ctr, face, 0.0
        elseif ks.set.space[1:2] == "2d"
            ctr, a1face, a2face = init_fvm(ks, ks.pSpace)
            return ks, ctr, a1face, a2face, 0.0
        end
    end
end

function initialize(config::T) where {T<:AbstractDict}
    println("==============================================================")
    println("Kinetic.jl")
    println("A Lightweight Toolbox for Kinetic Modeling and Simulation")
    println("==============================================================")
    println("")
    @info "initializing solver"
    println("")

    ks = SolverSet(config)

    if ks.set.space[1:2] == "1d"
        ctr, face = init_fvm(ks, ks.pSpace)
        return ks, ctr, face, 0.0
    elseif ks.set.space[1:2] == "2d"
        ctr, a1face, a2face = init_fvm(ks, ks.pSpace)
        return ks, ctr, a1face, a2face, 0.0
    end
end


"""
    init_fvm(KS::T, ps::T1) where {T<:AbstractSolverSet,T1<:AbstractPhysicalSpace}

Initialize finite volume method

"""
function init_fvm(KS::T, ps::T1, array=:static_array) where {T<:AbstractSolverSet,T1<:AbstractPhysicalSpace1D}
    funcar = eval(array)
    
    if KS.set.space[3:4] == "0f"

        ctr = OffsetArray{ControlVolume1D}(undef, eachindex(KS.pSpace.x))
        face = Array{Interface1D}(undef, KS.pSpace.nx + 1)

        for i in eachindex(ctr)
            if i <= KS.pSpace.nx ÷ 2
                ctr[i] =
                    ControlVolume1D(KS.pSpace.x[i], KS.pSpace.dx[i], funcar(KS.ib.wL), funcar(KS.ib.primL))
            else
                ctr[i] =
                    ControlVolume1D(KS.pSpace.x[i], KS.pSpace.dx[i], funcar(KS.ib.wR), funcar(KS.ib.primR))
            end
        end

        for i = 1:KS.pSpace.nx+1
            face[i] = Interface1D(KS.ib.wL)
        end

    elseif KS.set.space[3:4] == "1f"

        ctr = OffsetArray{ControlVolume1D1F}(undef, eachindex(KS.pSpace.x))
        face = Array{Interface1D1F}(undef, KS.pSpace.nx + 1)

        for i in eachindex(ctr)
            if i <= KS.pSpace.nx ÷ 2
                ctr[i] = ControlVolume1D1F(
                    KS.pSpace.x[i],
                    KS.pSpace.dx[i],
                    funcar(KS.ib.wL),
                    funcar(KS.ib.primL),
                    funcar(KS.ib.fL),
                )
            else
                ctr[i] = ControlVolume1D1F(
                    KS.pSpace.x[i],
                    KS.pSpace.dx[i],
                    funcar(KS.ib.wR),
                    funcar(KS.ib.primR),
                    funcar(KS.ib.fR),
                )
            end
        end

        for i = 1:KS.pSpace.nx+1
            face[i] = Interface1D1F(funcar(KS.ib.wL), funcar(KS.ib.fL))
        end

    elseif KS.set.space[3:4] == "2f"

        ctr = OffsetArray{ControlVolume1D2F}(undef, eachindex(KS.pSpace.x))
        face = Array{Interface1D2F}(undef, KS.pSpace.nx + 1)

        for i in eachindex(ctr)
            if i <= KS.pSpace.nx ÷ 2
                ctr[i] = ControlVolume1D2F(
                    KS.pSpace.x[i],
                    KS.pSpace.dx[i],
                    funcar(KS.ib.wL),
                    funcar(KS.ib.primL),
                    funcar(KS.ib.hL),
                    funcar(KS.ib.bL),
                )
            else
                ctr[i] = ControlVolume1D2F(
                    KS.pSpace.x[i],
                    KS.pSpace.dx[i],
                    funcar(KS.ib.wR),
                    funcar(KS.ib.primR),
                    funcar(KS.ib.hR),
                    funcar(KS.ib.bR),
                )
            end
        end

        for i = 1:KS.pSpace.nx+1
            face[i] = Interface1D2F(funcar(KS.ib.wL), funcar(KS.ib.hL))
        end

    elseif KS.set.space[3:4] == "3f"

        ctr = OffsetArray{ControlVolume1D3F}(undef, eachindex(KS.pSpace.x))
        face = Array{Interface1D3F}(undef, KS.pSpace.nx + 1)

        for i in eachindex(ctr)
            if i <= KS.pSpace.nx ÷ 2
                ctr[i] = ControlVolume1D3F(
                    KS.pSpace.x[i],
                    KS.pSpace.dx[i],
                    funcar(KS.ib.wL),
                    funcar(KS.ib.primL),
                    funcar(KS.ib.h0L),
                    funcar(KS.ib.h1L),
                    funcar(KS.ib.h2L),
                    funcar(KS.ib.EL),
                    funcar(KS.ib.BL),
                    funcar(KS.ib.lorenzL),
                )
            else
                ctr[i] = ControlVolume1D3F(
                    KS.pSpace.x[i],
                    KS.pSpace.dx[i],
                    funcar(KS.ib.wR),
                    funcar(KS.ib.primR),
                    funcar(KS.ib.h0R),
                    funcar(KS.ib.h1R),
                    funcar(KS.ib.h2R),
                    funcar(KS.ib.ER),
                    funcar(KS.ib.BR),
                    funcar(KS.ib.lorenzR),
                )
            end
        end

        for i = 1:KS.pSpace.nx+1
            face[i] = Interface1D3F(funcar(KS.ib.wL), funcar(KS.ib.h0L), funcar(KS.ib.EL))
        end

    elseif KS.set.space[3:4] == "4f"

        ctr = OffsetArray{ControlVolume1D4F}(undef, eachindex(KS.pSpace.x))
        face = Array{Interface1D4F}(undef, KS.pSpace.nx + 1)

        for i in eachindex(ctr)
            if i <= KS.pSpace.nx ÷ 2
                ctr[i] = ControlVolume1D4F(
                    KS.pSpace.x[i],
                    KS.pSpace.dx[i],
                    funcar(KS.ib.wL),
                    funcar(KS.ib.primL),
                    funcar(KS.ib.h0L),
                    funcar(KS.ib.h1L),
                    funcar(KS.ib.h2L),
                    funcar(KS.ib.h3L),
                    funcar(KS.ib.EL),
                    funcar(KS.ib.BL),
                    funcar(KS.ib.lorenzL),
                )
            else
                ctr[i] = ControlVolume1D4F(
                    KS.pSpace.x[i],
                    KS.pSpace.dx[i],
                    funcar(KS.ib.wR),
                    funcar(KS.ib.primR),
                    funcar(KS.ib.h0R),
                    funcar(KS.ib.h1R),
                    funcar(KS.ib.h2R),
                    funcar(KS.ib.h3R),
                    funcar(KS.ib.ER),
                    funcar(KS.ib.BR),
                    funcar(KS.ib.lorenzR),
                )
            end
        end

        for i = 1:KS.pSpace.nx+1
            face[i] = Interface1D4F(funcar(KS.ib.wL), funcar(KS.ib.h0L), funcar(KS.ib.EL))
        end

    end

    return ctr, face
end

function init_fvm(KS::T, ps::T1, array=:static_array) where {T<:AbstractSolverSet,T1<:AbstractPhysicalSpace2D}
    funcar = eval(array)

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
                    funcar(KS.ib.wL),
                    funcar(KS.ib.primL),
                    funcar(KS.ib.fL),
                )
            else
                ctr[i, j] = ControlVolume2D1F(
                    KS.pSpace.x[i, j],
                    KS.pSpace.y[i, j],
                    KS.pSpace.dx[i, j],
                    KS.pSpace.dy[i, j],
                    funcar(KS.ib.wR),
                    funcar(KS.ib.primR),
                    funcar(KS.ib.fR),
                )
            end
        end

        for j = 1:KS.pSpace.ny
            for i = 1:KS.pSpace.nx
                a1face[i, j] =
                    Interface2D1F(KS.pSpace.dy[i, j], 1.0, 0.0, funcar(KS.ib.wL), funcar(KS.ib.fL))
            end
            a1face[KS.pSpace.nx+1, j] =
                Interface2D1F(KS.pSpace.dy[KS.pSpace.nx, j], 1.0, 0.0, funcar(KS.ib.wL), funcar(KS.ib.fL))
        end
        for i = 1:KS.pSpace.nx
            for j = 1:KS.pSpace.ny
                a2face[i, j] =
                    Interface2D1F(KS.pSpace.dx[i, j], 0.0, 1.0, funcar(KS.ib.wL), funcar(KS.ib.fL))
            end
            a2face[i, KS.pSpace.ny+1] =
                Interface2D1F(KS.pSpace.dx[i, KS.pSpace.ny], 0.0, 1.0, funcar(KS.ib.wL), funcar(KS.ib.fL))
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
                    funcar(KS.ib.wL),
                    funcar(KS.ib.primL),
                    funcar(KS.ib.hL),
                    funcar(KS.ib.bL),
                )
            else
                ctr[i, j] = ControlVolume2D2F(
                    KS.pSpace.x[i, j],
                    KS.pSpace.y[i, j],
                    KS.pSpace.dx[i, j],
                    KS.pSpace.dy[i, j],
                    funcar(KS.ib.wR),
                    funcar(KS.ib.primR),
                    funcar(KS.ib.hR),
                    funcar(KS.ib.bR),
                )
            end
        end

        for j = 1:KS.pSpace.ny
            for i = 1:KS.pSpace.nx
                a1face[i, j] =
                    Interface2D2F(KS.pSpace.dy[i, j], 1.0, 0.0, funcar(KS.ib.wL), funcar(KS.ib.hL))
            end
            a1face[KS.pSpace.nx+1, j] =
                Interface2D2F(KS.pSpace.dy[KS.pSpace.nx, j], 1.0, 0.0, funcar(KS.ib.wL), funcar(KS.ib.hL))
        end
        for i = 1:KS.pSpace.nx
            for j = 1:KS.pSpace.ny
                a2face[i, j] =
                    Interface2D2F(KS.pSpace.dx[i, j], 0.0, 1.0, funcar(KS.ib.wL), funcar(KS.ib.hL))
            end
            a2face[i, KS.pSpace.ny+1] =
                Interface2D2F(KS.pSpace.dx[i, KS.pSpace.ny], 0.0, 1.0, funcar(KS.ib.wL), funcar(KS.ib.hL))
        end

    end

    return ctr, a1face, a2face
end

function init_fvm(KS::T, ps::UnstructPSpace) where {T<:AbstractSolverSet}
    if KS.set.space[3:4] == "0f"
        
        ctr = Array{KitBase.ControlVolumeUS}(undef, size(ps.cellid, 1))
        for i in eachindex(ctr)
            n = Vector{Float64}[]
            for j = 1:3
                push!(
                    n,
                    KitBase.unit_normal(
                        ps.points[ps.facePoints[ps.cellFaces[i, j], 1], :],
                        ps.points[ps.facePoints[ps.cellFaces[i, j], 2], :],
                    ) |> funcar,
                )

                if dot(
                    ps.faceCenter[ps.cellFaces[i, j], 1:2] .- ps.cellCenter[i, 1:2],
                    n[j],
                ) < 0
                    n[j] .= -n[j]
                end
            end

            dx = [
                KitBase.point_distance(
                    ps.cellCenter[i, :],
                    ps.points[ps.cellid[i, 1], :],
                    ps.points[ps.cellid[i, 2], :],
                ),
                KitBase.point_distance(
                    ps.cellCenter[i, :],
                    ps.points[ps.cellid[i, 2], :],
                    ps.points[ps.cellid[i, 3], :],
                ),
                KitBase.point_distance(
                    ps.cellCenter[i, :],
                    ps.points[ps.cellid[i, 3], :],
                    ps.points[ps.cellid[i, 1], :],
                ),
            ] |> funcar

            if ps.cellCenter[i, 1] <=
               minimum(ps.cellCenter[:, 1]) +
               (maximum(ps.cellCenter[:, 1]) - minimum(ps.cellCenter[:, 1])) / 2
                ctr[i] = KitBase.ControlVolumeUS(
                    n,
                    funcar(ps.cellCenter[i, :]),
                    dx,
                    funcar(KS.ib.wL),
                    funcar(KS.ib.primL),
                )
            else
                ctr[i] = KitBase.ControlVolumeUS(
                    n,
                    funcar(ps.cellCenter[i, :]),
                    dx,
                    funcar(KS.ib.wR),
                    funcar(KS.ib.primR),
                )
            end
        end

        face = Array{KitBase.Interface2D}(undef, size(ps.facePoints, 1))
        for i in eachindex(face)
            len =
                norm(ps.points[ps.facePoints[i, 1], :] .- ps.points[ps.facePoints[i, 2], :])
            n = KitBase.unit_normal(
                ps.points[ps.facePoints[i, 1], :],
                ps.points[ps.facePoints[i, 2], :],
            )

            if !(-1 in ps.faceCells[i, :])
                n0 =
                    ps.cellCenter[ps.faceCells[i, 2], :] .-
                    ps.cellCenter[ps.faceCells[i, 1], :]
            else
                idx =
                    ifelse(ps.faceCells[i, 1] != -1, ps.faceCells[i, 1], ps.faceCells[i, 2])
                n0 = ps.cellCenter[idx, :] .- ps.faceCenter[i, :]
            end
            if dot(n, n0[1:2]) < 0
                n .= -n
            end

            fw = zero(KS.ib.wL) |> funcar

            face[i] = KitBase.Interface2D(len, n[1], n[2], fw)
        end
    
    elseif KS.set.space[3:4] == "1f"
    
        ctr = Array{KitBase.ControlVolumeUS1F}(undef, size(ps.cellid, 1))
        for i in eachindex(ctr)
            n = Vector{Float64}[]
            for j = 1:3
                push!(
                    n,
                    KitBase.unit_normal(
                        ps.points[ps.facePoints[ps.cellFaces[i, j], 1], :],
                        ps.points[ps.facePoints[ps.cellFaces[i, j], 2], :],
                    ) |> funcar,
                )

                if dot(
                    ps.faceCenter[ps.cellFaces[i, j], 1:2] .- ps.cellCenter[i, 1:2],
                    n[j],
                ) < 0
                    n[j] .= -n[j]
                end
            end

            dx = [
                KitBase.point_distance(
                    ps.cellCenter[i, :],
                    ps.points[ps.cellid[i, 1], :],
                    ps.points[ps.cellid[i, 2], :],
                ),
                KitBase.point_distance(
                    ps.cellCenter[i, :],
                    ps.points[ps.cellid[i, 2], :],
                    ps.points[ps.cellid[i, 3], :],
                ),
                KitBase.point_distance(
                    ps.cellCenter[i, :],
                    ps.points[ps.cellid[i, 3], :],
                    ps.points[ps.cellid[i, 1], :],
                ),
            ] |> funcar

            if ps.cellCenter[i, 1] <=
               minimum(ps.cellCenter[:, 1]) +
               (maximum(ps.cellCenter[:, 1]) - minimum(ps.cellCenter[:, 1])) / 2
                ctr[i] = KitBase.ControlVolumeUS1F(
                    n,
                    ps.cellCenter[i, :],
                    dx,
                    KS.ib.wL,
                    KS.ib.primL,
                    KS.ib.fL,
                )
            else
                ctr[i] = KitBase.ControlVolumeUS1F(
                    n,
                    funcar(ps.cellCenter[i, :]),
                    dx,
                    KS.ib.wR,
                    funcar(KS.ib.primR),
                    funcar(KS.ib.fR),
                )
            end
        end

        face = Array{KitBase.Interface2D1F}(undef, size(ps.facePoints, 1))
        for i in eachindex(face)
            len =
                norm(ps.points[ps.facePoints[i, 1], :] .- ps.points[ps.facePoints[i, 2], :])
            n = KitBase.unit_normal(
                ps.points[ps.facePoints[i, 1], :],
                ps.points[ps.facePoints[i, 2], :],
            )

            if !(-1 in ps.faceCells[i, :])
                n0 =
                    ps.cellCenter[ps.faceCells[i, 2], :] .-
                    ps.cellCenter[ps.faceCells[i, 1], :]
            else
                idx =
                    ifelse(ps.faceCells[i, 1] != -1, ps.faceCells[i, 1], ps.faceCells[i, 2])
                n0 = ps.cellCenter[idx, :] .- ps.faceCenter[i, :]
            end
            if dot(n, n0[1:2]) < 0
                n .= -n
            end

            fw = zero(KS.ib.wL)
            ff = zero(KS.ib.fL)

            face[i] = KitBase.Interface2D1F(len, n[1], n[2], funcar(fw), funcar(ff))
        end
    
    elseif KS.set.space[3:4] == "2f"

        ctr = Array{KitBase.ControlVolumeUS2F}(undef, size(ps.cellid, 1))
        for i in eachindex(ctr)
            n = Vector{Float64}[]
            for j = 1:3
                push!(
                    n,
                    KitBase.unit_normal(
                        ps.points[ps.facePoints[ps.cellFaces[i, j], 1], :],
                        ps.points[ps.facePoints[ps.cellFaces[i, j], 2], :],
                    ) |> funcar,
                )

                if dot(
                    ps.faceCenter[ps.cellFaces[i, j], 1:2] .- ps.cellCenter[i, 1:2],
                    n[j],
                ) < 0
                    n[j] .= -n[j]
                end
            end

            dx = [
                KitBase.point_distance(
                    ps.cellCenter[i, :],
                    ps.points[ps.cellid[i, 1], :],
                    ps.points[ps.cellid[i, 2], :],
                ),
                KitBase.point_distance(
                    ps.cellCenter[i, :],
                    ps.points[ps.cellid[i, 2], :],
                    ps.points[ps.cellid[i, 3], :],
                ),
                KitBase.point_distance(
                    ps.cellCenter[i, :],
                    ps.points[ps.cellid[i, 3], :],
                    ps.points[ps.cellid[i, 1], :],
                ),
            ] |> funcar

            if ps.cellCenter[i, 1] <=
               minimum(ps.cellCenter[:, 1]) +
               (maximum(ps.cellCenter[:, 1]) - minimum(ps.cellCenter[:, 1])) / 2
                ctr[i] = KitBase.ControlVolumeUS2F(
                    n,
                    funcar(ps.cellCenter[i, :]),
                    dx,
                    funcar(KS.ib.wL),
                    funcar(KS.ib.primL),
                    funcar(KS.ib.hL),
                    funcar(KS.ib.bL),
                )
            else
                ctr[i] = KitBase.ControlVolumeUS2F(
                    n,
                    funcar(ps.cellCenter[i, :]),
                    dx,
                    funcar(KS.ib.wR),
                    funcar(KS.ib.primR),
                    funcar(KS.ib.hR),
                    funcar(KS.ib.bR),
                )
            end
        end

        face = Array{KitBase.Interface2D2F}(undef, size(ps.facePoints, 1))
        for i in eachindex(face)
            len =
                norm(ps.points[ps.facePoints[i, 1], :] .- ps.points[ps.facePoints[i, 2], :])
            n = KitBase.unit_normal(
                ps.points[ps.facePoints[i, 1], :],
                ps.points[ps.facePoints[i, 2], :],
            )

            if !(-1 in ps.faceCells[i, :])
                n0 =
                    ps.cellCenter[ps.faceCells[i, 2], :] .-
                    ps.cellCenter[ps.faceCells[i, 1], :]
            else
                idx =
                    ifelse(ps.faceCells[i, 1] != -1, ps.faceCells[i, 1], ps.faceCells[i, 2])
                n0 = ps.cellCenter[idx, :] .- ps.faceCenter[i, :]
            end
            if dot(n, n0[1:2]) < 0
                n .= -n
            end

            fw = zero(KS.ib.wL)
            fh = zero(KS.ib.hL)

            face[i] = KitBase.Interface2D2F(len, n[1], n[2], funcar(fw), funcar(fh))
        end

    end

    return ctr, face
end


"""
    init_ptc!(
        KS::SolverSet,
        ctr::T;
        mode = :soa::Symbol,
        factor = 1::Real,
    ) where {T<:AbstractArray{<:AbstractControlVolume1D,1}}

Initialize particles based on flow conditions

"""
function init_ptc!(
    KS::SolverSet,
    ctr::T;
    mode = :soa::Symbol,
    factor = 1::Real,
) where {T<:AbstractArray{<:AbstractControlVolume1D,1}}
    if mode == :soa
        init_ptc_soa!(KS, ctr, factor)
    elseif mode == :aos
        init_ptc_aos!(KS, ctr, factor)
    end
end


"""
    init_ptc_aos!(
        KS::SolverSet,
        ctr::T,
        factor = 1,
    ) where {T<:AbstractArray{<:AbstractControlVolume1D,1}}

Initialize particles with array of structs

"""
function init_ptc_aos!(
    KS::SolverSet,
    ctr::T,
    factor = 1,
) where {T<:AbstractArray{<:AbstractControlVolume1D,1}}

    np = 0
    for i in eachindex(ctr)
        np += Int(round(ctr[i].w[1] * ctr[i].dx / KS.gas.m))
    end
    KS.gas.np = np
    np *= round(factor) |> Int

    ptc = Array{Particle1D}(undef, np)
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
    init_ptc_soa!(
        KS::SolverSet,
        ctr::T,
        factor = 1,
    ) where {T<:AbstractArray{<:AbstractControlVolume1D,1}}

Initialize particles with struct of arrays

"""
function init_ptc_soa!(
    KS::SolverSet,
    ctr::T,
    factor = 1,
) where {T<:AbstractArray{<:AbstractControlVolume1D,1}}

    np = 0
    for i in eachindex(ctr)
        np += round(ctr[i].prim[1] * ctr[i].dx / KS.gas.m) |> Int
    end
    KS.gas.np = np
    np *= round(factor) |> Int
    np_tmp = 0

    m = zeros(np)
    x = zeros(np)
    v = zeros(np, 3)
    e = zeros(np)
    idx = zeros(Int, np)
    ref = zeros(Int, np) # default
    flag = zeros(Int, np) # default
    tc = zeros(np)

    # in-cell particles
    for i in eachindex(ctr)
        npl = Int(round(ctr[i].w[1] * ctr[i].dx / KS.gas.m))
        for j = 1:npl
            np_tmp += 1

            m[np_tmp] = KS.gas.m
            x[np_tmp] = ctr[i].x + (rand() - 0.5) * ctr[i].dx
            v[np_tmp, :] .= sample_maxwell(ctr[i].prim)
            e[np_tmp] = 0.5 / ctr[i].prim[end]
            idx[np_tmp] = i
            flag[np_tmp] = 0
            if i < 1
                flag[np_tmp] = 1
            elseif i > KS.pSpace.nx
                flag[np_tmp] = 2
            end
            τ = vhs_collision_time(ctr[i].prim, KS.gas.μᵣ, KS.gas.ω)
            tc[np_tmp] = next_collision_time(τ)
        end
    end

    # placeholder particles
    for i = KS.gas.np+1:np
        m[i] = KS.gas.m
        x[i] = 0.0
        idx[i] = -7
        v[i, :] .= 0.0
        e[i] = 0.0
        τ = -1e8
        tc[i] = 1e8
    end

    ptc = Particle(m, x, v, e, idx, ref, flag, tc)

    return ptc

end
