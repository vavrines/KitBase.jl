# ============================================================
# Solver Initializer
# ============================================================

"""
$(TYPEDSIGNATURES)

Initialize solver...
"""
function initialize(configfilename::T) where {T<:AbstractString}
    println("==============================================================")
    println("Kinetic.jl")
    println("A Portable Toolbox for Scientific and Neural Computing")
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
        cftuple = init_fvm(ks)

        return (ks, cftuple..., 0.0)
    end
end

"""
$(TYPEDSIGNATURES)
"""
function initialize(config::T) where {T<:AbstractDict}
    println("==============================================================")
    println("Kinetic.jl")
    println("A Portable Toolbox for Scientific and Neural Computing")
    println("==============================================================")
    println("")
    @info "initializing solver"
    println("")

    ks = SolverSet(config)
    cftuple = init_fvm(ks)

    return (ks, cftuple..., 0.0)
end


"""
$(TYPEDSIGNATURES)

Initialize finite volume method...
"""
function init_fvm(
    KS::T,
    array = :dynamic_array;
    structarray = false,
) where {T<:AbstractSolverSet}
    return init_fvm(KS, KS.ps, array; structarray = structarray)
end

function init_fvm(
    KS::T,
    ps::T1,
    array = :dynamic_array;
    structarray = false,
) where {T<:AbstractSolverSet,T1<:AbstractPhysicalSpace1D}

    funcar = eval(array)
    funcst = ifelse(structarray, StructArray, dynamic_array)

    funcprim = ifelse(KS.set.nSpecies == 1, conserve_prim, mixture_conserve_prim)
    γ = begin
        if KS.gas isa Scalar
            KS.gas.a
        else
            KS.gas.γ
        end
    end

    if KS.set.space[3:4] == "0f"

        ctr = OffsetArray{ControlVolume}(undef, eachindex(KS.pSpace.x))
        face = Array{Interface}(undef, KS.pSpace.nx + 1)

        for i in eachindex(ctr)
            w = KS.ib.fw(KS.ps.x[i], KS.ib.p)
            prim = begin
                if γ == 0
                    funcprim(w)
                else
                    funcprim(w, γ)
                end
            end
            ctr[i] = ControlVolume(funcar(w), funcar(prim), 1)
        end

        for i = 1:KS.pSpace.nx+1
            fw = deepcopy(KS.ib.fw(KS.ps.x[1], KS.ib.p)) |> funcar
            face[i] = Interface(fw, 1)
        end

    elseif KS.set.space[3:4] == "1f"

        ctr = OffsetArray{ControlVolume1F}(undef, eachindex(KS.pSpace.x))
        face = Array{Interface1F}(undef, KS.pSpace.nx + 1)

        for i in eachindex(ctr)
            w = KS.ib.fw(KS.ps.x[i], KS.ib.p)
            prim = funcprim(w, γ)
            f = KS.ib.ff(KS.ps.x[i], KS.ib.p)
            ctr[i] = ControlVolume(funcar(w), funcar(prim), funcar(f), 1)
        end

        for i = 1:KS.pSpace.nx+1
            fw = deepcopy(KS.ib.fw(KS.ps.x[1], KS.ib.p)) |> funcar
            ff = deepcopy(KS.ib.ff(KS.ps.x[1], KS.ib.p)) |> funcar
            face[i] = Interface(fw, ff, 1)
        end

    elseif KS.set.space[3:4] == "2f"

        ctr = OffsetArray{ControlVolume2F}(undef, eachindex(KS.pSpace.x))
        face = Array{Interface2F}(undef, KS.pSpace.nx + 1)

        for i in eachindex(ctr)
            w = KS.ib.fw(KS.ps.x[i], KS.ib.p)
            prim = funcprim(w, γ)
            h, b = KS.ib.ff(KS.ps.x[i], KS.ib.p)
            ctr[i] = ControlVolume(funcar(w), funcar(prim), funcar(h), funcar(b), 1)
        end

        for i = 1:KS.pSpace.nx+1
            fw = deepcopy(KS.ib.fw(KS.ps.x[1], KS.ib.p)) |> funcar
            ff = deepcopy(KS.ib.ff(KS.ps.x[1], KS.ib.p)[1]) |> funcar
            face[i] = Interface(fw, ff, ff, 1)
        end

    elseif KS.set.space[3:4] == "3f"

        ctr = OffsetArray{ControlVolume1D3F}(undef, eachindex(KS.pSpace.x))
        face = Array{Interface1D3F}(undef, KS.pSpace.nx + 1)

        for i in eachindex(ctr)
            w = KS.ib.fw(KS.ps.x[i], KS.ib.p)
            prim = funcprim(w, γ)
            h0, h1, h2 = KS.ib.ff(KS.ps.x[i], KS.ib.p)
            E = KS.ib.fE(KS.ps.x[i], KS.ib.p)
            B = KS.ib.fB(KS.ps.x[i], KS.ib.p)
            L = KS.ib.fL(KS.ps.x[i], KS.ib.p)

            ctr[i] = ControlVolume1D3F(
                funcar(w),
                funcar(prim),
                funcar(h0),
                funcar(h1),
                funcar(h2),
                funcar(E),
                funcar(B),
                funcar(L),
            )
        end

        for i = 1:KS.pSpace.nx+1
            fw = deepcopy(KS.ib.fw(KS.ps.x[1], KS.ib.p)) |> funcar
            ff = deepcopy(KS.ib.ff(KS.ps.x[1], KS.ib.p)[1]) |> funcar
            fe = deepcopy(KS.ib.fE(KS.ps.x[1], KS.ib.p)) |> funcar
            face[i] = Interface1D3F(fw, ff, fe)
        end

    elseif KS.set.space[3:4] == "4f"

        ctr = OffsetArray{ControlVolume1D4F}(undef, eachindex(KS.pSpace.x))
        face = Array{Interface1D4F}(undef, KS.pSpace.nx + 1)

        for i in eachindex(ctr)
            w = KS.ib.fw(KS.ps.x[i], KS.ib.p)
            prim = funcprim(w, γ)
            h0, h1, h2, h3 = KS.ib.ff(KS.ps.x[i], KS.ib.p)
            E = KS.ib.fE(KS.ps.x[i], KS.ib.p)
            B = KS.ib.fB(KS.ps.x[i], KS.ib.p)
            L = KS.ib.fL(KS.ps.x[i], KS.ib.p)

            ctr[i] = ControlVolume1D4F(
                funcar(w),
                funcar(prim),
                funcar(h0),
                funcar(h1),
                funcar(h2),
                funcar(h3),
                funcar(E),
                funcar(B),
                funcar(L),
            )
        end

        for i = 1:KS.pSpace.nx+1
            fw = deepcopy(KS.ib.fw(KS.ps.x[1], KS.ib.p)) |> funcar
            ff = deepcopy(KS.ib.ff(KS.ps.x[1], KS.ib.p)[1]) |> funcar
            fe = deepcopy(KS.ib.fE(KS.ps.x[1], KS.ib.p)) |> funcar
            face[i] = Interface1D4F(fw, ff, fe)
        end

    end

    return ctr |> funcst, face |> funcst
end

function init_fvm(
    KS::T,
    ps::T1,
    array = :dynamic_array;
    structarray = false,
) where {T<:AbstractSolverSet,T1<:AbstractPhysicalSpace2D}

    funcar = eval(array)
    funcst = ifelse(structarray, StructArray, dynamic_array)
    funcprim = ifelse(KS.set.nSpecies == 1, conserve_prim, mixture_conserve_prim)

    nx, ny, dx, dy = begin
        if ps isa CSpace2D
            ps.nr, ps.nθ, ps.dr, ps.darc
        else
            ps.nx, ps.ny, ps.dx, ps.dy
        end
    end

    if KS.set.space[3:4] == "0f"

        ctr = OffsetArray{ControlVolume}(undef, axes(KS.pSpace.x, 1), axes(KS.pSpace.y, 2))
        a1face = Array{Interface}(undef, nx + 1, ny)
        a2face = Array{Interface}(undef, nx, ny + 1)

        for j in axes(ctr, 2), i in axes(ctr, 1)
            w = KS.ib.fw(KS.ps.x[i, j], KS.ps.y[i, j], KS.ib.p)
            prim = funcprim(w, KS.gas.γ)

            ctr[i, j] = ControlVolume(funcar(w), funcar(prim), 2)
        end

        for j = 1:ny
            for i = 1:nx+1
                a1face[i, j] =
                    Interface(funcar(KS.ib.fw(KS.ps.x[1], KS.ps.y[1], KS.ib.p)), 2)
            end
        end
        for i = 1:nx
            for j = 1:ny+1
                a2face[i, j] =
                    Interface(funcar(KS.ib.fw(KS.ps.x[1], KS.ps.y[1], KS.ib.p)), 2)
            end
        end

    elseif KS.set.space[3:4] == "1f"

        ctr =
            OffsetArray{ControlVolume1F}(undef, axes(KS.pSpace.x, 1), axes(KS.pSpace.y, 2))
        a1face = Array{Interface1F}(undef, nx + 1, ny)
        a2face = Array{Interface1F}(undef, nx, ny + 1)

        for j in axes(ctr, 2), i in axes(ctr, 1)
            w = KS.ib.fw(KS.ps.x[i, j], KS.ps.y[i, j], KS.ib.p)
            prim = funcprim(w, KS.gas.γ)
            h = KS.ib.ff(KS.ps.x[i, j], KS.ps.y[i, j], KS.ib.p)
            ctr[i, j] = ControlVolume(funcar(w), funcar(prim), funcar(h), 2)
        end

        for j = 1:ny
            for i = 1:nx+1
                a1face[i, j] = Interface(
                    funcar(KS.ib.fw(KS.ps.x[1], KS.ps.y[1], KS.ib.p)),
                    funcar(KS.ib.ff(KS.ps.x[1], KS.ps.y[1], KS.ib.p)),
                    2,
                )
            end
        end
        for i = 1:nx
            for j = 1:ny+1
                a2face[i, j] = Interface(
                    funcar(KS.ib.fw(KS.ps.x[1], KS.ps.y[1], KS.ib.p)),
                    funcar(KS.ib.ff(KS.ps.x[1], KS.ps.y[1], KS.ib.p)),
                    2,
                )
            end
        end

    elseif KS.set.space[3:4] == "2f"

        ctr =
            OffsetArray{ControlVolume2F}(undef, axes(KS.pSpace.x, 1), axes(KS.pSpace.y, 2))
        a1face = Array{Interface2F}(undef, nx + 1, ny)
        a2face = Array{Interface2F}(undef, nx, ny + 1)

        for j in axes(ctr, 2), i in axes(ctr, 1)
            w = KS.ib.fw(KS.ps.x[i, j], KS.ps.y[i, j], KS.ib.p)
            prim = funcprim(w, KS.gas.γ)
            h, b = KS.ib.ff(KS.ps.x[i, j], KS.ps.y[i, j], KS.ib.p)

            ctr[i, j] = ControlVolume(funcar(w), funcar(prim), funcar(h), funcar(b), 2)
        end

        for j = 1:ny
            for i = 1:nx+1
                a1face[i, j] = Interface(
                    funcar(KS.ib.fw(KS.ps.x[1], KS.ps.y[1], KS.ib.p)),
                    funcar(KS.ib.ff(KS.ps.x[1], KS.ps.y[1], KS.ib.p)[1]),
                    funcar(KS.ib.ff(KS.ps.x[1], KS.ps.y[1], KS.ib.p)[2]),
                    2,
                )
            end
        end
        for i = 1:nx
            for j = 1:ny+1
                a2face[i, j] = Interface(
                    funcar(KS.ib.fw(KS.ps.x[1], KS.ps.y[1], KS.ib.p)),
                    funcar(KS.ib.ff(KS.ps.x[1], KS.ps.y[1], KS.ib.p)[1]),
                    funcar(KS.ib.ff(KS.ps.x[1], KS.ps.y[1], KS.ib.p)[2]),
                    2,
                )
            end
        end

    end

    return ctr |> funcst, a1face |> funcst, a2face |> funcst

end

function init_fvm(
    KS::T,
    ps::UnstructPSpace,
    array = :dynamic_array;
    structarray = false,
) where {T<:AbstractSolverSet}

    funcar = eval(array)
    funcst = ifelse(structarray, StructArray, dynamic_array)
    funcprim = ifelse(KS.set.nSpecies == 1, conserve_prim, mixture_conserve_prim)

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

            dx =
                [
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

            w = KS.ib.fw(ps.cellCenter[i, 1], ps.cellCenter[i, 2], KS.ib.p)
            prim = funcprim(w, KS.gas.γ)

            ctr[i] = KitBase.ControlVolumeUS(
                n,
                funcar(ps.cellCenter[i, :]),
                dx,
                funcar(w),
                funcar(prim),
            )
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

            fw = zero(KS.ib.fw(ctr[1].x[1], ctr[1].x[2], KS.ib.p)) |> funcar

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

            dx =
                [
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

            w = KS.ib.fw(ps.cellCenter[i, 1], ps.cellCenter[i, 2], KS.ib.p)
            prim = funcprim(w, KS.gas.γ)
            h = KS.ib.ff(ps.cellCenter[i, 1], ps.cellCenter[i, 2], KS.ib.p)

            ctr[i] = KitBase.ControlVolumeUS1F(
                n,
                ps.cellCenter[i, :],
                dx,
                funcar(w),
                funcar(prim),
                funcar(h),
            )

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

            fw = zero(KS.ib.fw(ctr[1].x[1], ctr[1].x[2], KS.ib.p))
            ff = zero(KS.ib.ff(ctr[1].x[1], ctr[1].x[2], KS.ib.p))

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

            dx =
                [
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

            w = KS.ib.fw(ps.cellCenter[i, 1], ps.cellCenter[i, 2], KS.ib.p)
            prim = funcprim(w, KS.gas.γ)
            h, b = KS.ib.ff(ps.cellCenter[i, 1], ps.cellCenter[i, 2], KS.ib.p)

            ctr[i] = KitBase.ControlVolumeUS2F(
                n,
                funcar(ps.cellCenter[i, :]),
                dx,
                funcar(w),
                funcar(prim),
                funcar(h),
                funcar(b),
            )
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

            fw = zero(KS.ib.fw(ctr[1].x[1], ctr[1].x[2], KS.ib.p))
            fh = zero(KS.ib.ff(ctr[1].x[1], ctr[1].x[2], KS.ib.p)[1])

            face[i] = KitBase.Interface2D2F(len, n[1], n[2], funcar(fw), funcar(fh))
        end

    end

    return ctr |> funcst, face |> funcst

end


"""
$(SIGNATURES)

Initialize particles based on flow conditions
"""
function init_ptc!(
    KS::SolverSet,
    ctr::AV{T};
    mode = :soa::Symbol,
    factor = 1::Real,
) where {T<:AbstractControlVolume}
    if mode == :soa
        init_ptc_soa!(KS, ctr, factor)
    elseif mode == :aos
        init_ptc_aos!(KS, ctr, factor)
    end
end


"""
$(SIGNATURES)

Initialize particles with array of structs
"""
function init_ptc_aos!(
    KS::SolverSet,
    ctr::AV{T},
    factor = 1,
) where {T<:AbstractControlVolume}

    np = 0
    for i in eachindex(ctr)
        np += Int(round(ctr[i].w[1] * KS.ps.dx[i] / KS.gas.m))
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
        npl = Int(round(ctr[i].w[1] * KS.ps.dx[i] / KS.gas.m))
        for j = 1:npl
            np_tmp += 1
            sample_particle!(ptc[np_tmp], KS, ctr[i], i, KS.ps.x[i], KS.ps.dx[i])
        end
    end

    return ptc |> StructArray

end


"""
$(SIGNATURES)

Initialize particles with struct of arrays
"""
function init_ptc_soa!(
    KS::SolverSet,
    ctr::AV{T},
    factor = 1,
) where {T<:AbstractControlVolume}

    np = 0
    for i in eachindex(ctr)
        np += round(ctr[i].prim[1] * KS.ps.dx[i] / KS.gas.m) |> Int
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
        npl = Int(round(ctr[i].w[1] * KS.ps.dx[i] / KS.gas.m))
        for j = 1:npl
            np_tmp += 1

            m[np_tmp] = KS.gas.m
            x[np_tmp] = KS.ps.x[i] + (rand() - 0.5) * KS.ps.dx[i]
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
