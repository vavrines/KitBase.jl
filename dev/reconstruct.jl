using LinearAlgebra, ProgressMeter, JLD2
import KitBase

cd(@__DIR__)
D = KitBase.read_dict("naca_coarse.txt")
set = KitBase.set_setup(D)

ps = KitBase.set_geometry(D)
for i in eachindex(ps.cellType)
    if ps.cellType[i] == 1 &&
       -0.5 < ps.cellCenter[i, 1] < 1.5 &&
       -0.1 < ps.cellCenter[i, 2] < 0.1
        ps.cellType[i] = 2

        for j = 1:3
            if ps.edgeType[ps.cellEdges[i, j]] == 1
                ps.edgeType[ps.cellEdges[i, j]] = 2
            end
        end
    elseif ps.cellType[i] == 1 && ps.cellCenter[i, 1] > 0.5
        ps.cellType[i] = 3

        for j = 1:3
            if ps.edgeType[ps.cellEdges[i, j]] == 1
                ps.edgeType[ps.cellEdges[i, j]] = 3
            end
        end
    end
end

vs = KitBase.set_velocity(D)
gas = KitBase.set_property(D)
c0 = KitBase.sound_speed(1.0, gas.γ) * gas.Ma
α = π / 180 * 7.0

begin
    primL = [1.0, c0 * cos(α), c0 * sin(α), 1.0]
    #primL = [1.0, KitBase.sound_speed(1.0, gas.γ) * gas.Ma, 0.0, 1.0]
    wL = KitBase.prim_conserve(primL, gas.γ)
    hL = KitBase.maxwellian(vs.u, vs.v, primL)
    bL = @. hL * gas.K / 2 / primL[end]
    primR = [1.0, 0.0, 0.0, 1.0]
    wR = KitBase.prim_conserve(primR, gas.γ)
    hR = KitBase.maxwellian(vs.u, vs.v, primR)
    bR = @. hR * gas.K / 2 / primR[end]
    ib = KitBase.IB2F(wL, primL, hL, bL, primL, wL, primL, hL, bL, primR)
end

ks = KitBase.SolverSet(set, ps, vs, gas, ib, @__DIR__)
ctr, face = KitBase.init_fvm(ks, ks.pSpace)

KitBase.extract_last(randn(3, 2), 1, mode = :view)


KitBase.reconstruct!(ks, ctr)



#@inbounds Threads.@threads for i in eachindex(ctr)
i = 1
    ids = ps.cellNeighbors[i, :]
    deleteat!(ids, findall(x -> x == -1, ids))
    id1, id2 = ids[1:2]

    sw = KitBase.extract_last(ctr[id1].sw, 1, mode = :view)

    wL = ctr[id1].w
    wR = ctr[id2].w

    dxL = ctr[i].x[1] - ctr[id1].x[1]
    dxR = ctr[id2].x[1] - ctr[i].x[1]

    KitBase.reconstruct3!(
        sw,
        wL,
        ctr[i].w,
        wR,
        dxL,
        dxR,
        Symbol(ks.set.limiter),
    )
    
end






for iter = 1:10
    if iter == 7
        continue
    end
    @show iter
end