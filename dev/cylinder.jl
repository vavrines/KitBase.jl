using LinearAlgebra, WriteVTK, ProgressMeter, JLD2, PyCall
import KitBase

cd(@__DIR__)
meshio = pyimport("meshio")
m = meshio.read("../assets/mesh/cylinder_circle.msh")

D = KitBase.read_dict("cylinder.txt")
set = KitBase.set_setup(D)

ps = KitBase.set_geometry(D)
for i in eachindex(ps.edgeType)
    i1 = ps.edgePoints[i, 1]
    i2 = ps.edgePoints[i, 2]

    if i1 in [ps.cells.index[1]; ps.cells.index[2]] && i2 in [ps.cells.index[1]; ps.cells.index[2]]
        ps.edgeType[i] = 2
        c1 = ps.edgeCells[i, 1]
        c2 = ps.edgeCells[i, 2]

        if c1 != -1
            ps.cellType[c1] = 2
        elseif c2 != -1
            ps.cellType[c2] = 2
        else
            throw("index error")
        end
    end
end

vs = KitBase.set_velocity(D)
gas = KitBase.set_property(D)

begin
    primL = [1.0, KitBase.sound_speed(1.0, gas.γ) * gas.Ma, 0.0, 1.0]
    wL = KitBase.prim_conserve(primL, gas.γ)
    hL = KitBase.maxwellian(vs.u, vs.v, primL)
    bL = @. hL * gas.K / 2 / primL[end]
    primR = [1.0, 0.0, 0.0, 1.0]
    wR = KitBase.prim_conserve(primR, gas.γ)
    hR = KitBase.maxwellian(vs.u, vs.v, primR)
    bR = @. hR * gas.K / 2 / primR[end]
    ib = KitBase.IB2F(
        wL,
        primL,
        hL,
        bL,
        primL,
        wR,
        primR,
        hR,
        bR,
        primR,
    )
end

ks = KitBase.SolverSet(set, ps, vs, gas, ib, @__DIR__)

ctr = Array{KitBase.ControlVolumeUS2F}(undef, size(ps.cellid, 1))
for i in eachindex(ctr)
    n = Vector{Float64}[]
    for j = 1:3
        push!(n, KitBase.unit_normal(ps.points[ps.edgePoints[ps.cellEdges[i, j], 1], :], ps.points[ps.edgePoints[ps.cellEdges[i, j], 2], :]))

        if dot(ps.edgeCenter[ps.cellEdges[i, j], 1:2] .- ps.cellCenter[i, 1:2], n[j]) < 0
            n[j] .= -n[j]
        end
    end

    dx = [
        KitBase.point_distance(ps.cellCenter[i, :], ps.points[ps.cellid[i, 1], :], ps.points[ps.cellid[i, 2], :]),
        KitBase.point_distance(ps.cellCenter[i, :], ps.points[ps.cellid[i, 2], :], ps.points[ps.cellid[i, 3], :]),
        KitBase.point_distance(ps.cellCenter[i, :], ps.points[ps.cellid[i, 3], :], ps.points[ps.cellid[i, 1], :]),
    ]

    ctr[i] = KitBase.ControlVolumeUS2F(
        n,
        ps.cellCenter[i, :],
        dx,
        wL,
        primL,
        hL,
        bL,
    )
end

face = Array{KitBase.Interface2D2F}(undef, size(ps.edgePoints, 1))
for i in eachindex(face)
    len = norm(ps.points[ps.edgePoints[i, 1], :] .- ps.points[ps.edgePoints[i, 2], :])
    n = KitBase.unit_normal(ps.points[ps.edgePoints[i, 1], :], ps.points[ps.edgePoints[i, 2], :])
    
    if !(-1 in ps.edgeCells[i, :])
        n0 = ps.cellCenter[ps.edgeCells[i, 2], :] .- ps.cellCenter[ps.edgeCells[i, 1], :]
    else
        idx = ifelse(ps.edgeCells[i, 1] != -1, ps.edgeCells[i, 1], ps.edgeCells[i, 2])
        n0 = ps.cellCenter[idx, :] .- ps.edgeCenter[i, :]
    end
    if dot(n, n0[1:2]) < 0
        n .= -n
    end
    
    fw = zero(ks.ib.wL)
    fh = zero(ks.ib.hL)

    face[i] = KitBase.Interface2D2F(
        len,
        n[1],
        n[2],
        fw,
        fh,
    )
end