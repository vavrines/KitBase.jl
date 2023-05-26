using KitBase, Plots
using Base.Threads: @threads
using KitBase.ProgressMeter: @showprogress

set = Setup(
    case = "cylinder",
    space = "2d0f0v",
    boundary = ["fix", "extra", "mirror", "extra"],
    limiter = "minmod",
    cfl = 0.1,
    maxTime = 2.0, # time
    flux = "hll",
    hasForce = true,
)
ps = PSpace2D(0, 3, 60, 0, 2, 40, 1, 1)
vs = VSpace2D(-3, 3, 28, -3, 3, 28)
gas = Gas(Kn = 1e-3, Ma = 0.8, K = 1.0)

prim0 = [1.0, 0.0, 0.0, 1.0]
prim1 = [1.0, gas.Ma * sound_speed(1.0, gas.γ), 0.0, 1.0]
fw = (args...) -> prim_conserve(prim1, gas.γ)
ff = function (args...)
    prim = conserve_prim(fw(args...), gas.γ)
    h = maxwellian(vs.u, vs.v, prim)
    b = h .* gas.K / 2 / prim[end]
    return h, b
end
bc = function (x, y, args...)
    if abs((x - 3)^2 + y^2 - 1) < 1e-3
        return prim0
    else
        return prim1
    end
end
ib = IB2F(fw, ff, bc, NamedTuple())

ks = SolverSet(set, ps, nothing, gas, ib)
ctr, a1face, a2face = init_fvm(ks; structarray = true)

radius = 1
θs = linspace(0, π / 2, 30)
xs = @. ps.x1 - radius * cos(θs)
ys = @. sin(θs)
lps = hcat(xs, ys)
h = ps.dx[1]
ΔV = 2π * radius / 4 / length(θs) * h

#scatter(xs, ys, xticks = 0:0.05:3, yticks = 0:0.05:1, ratio = 1)
#savefig("cylinder.pdf")

# 1: fluid
# 0: solid
# -1: outer
# -2: ghost
flags = ones(Int, axes(ps.x))
for i in axes(flags, 1), j in axes(flags, 2)
    if (ps.x[i, j] - 3)^2 + ps.y[i, j]^2 < 1 # (x-3)^2+y^2=1
        flags[i, j] = 0
    end
end
flags[0, :] .= -1
flags[ks.ps.nx+1, :] .= -1
flags[:, 0] .= -1
flags[:, ks.ps.ny+1] .= -1

function ib_flag!(flags, ps)
    for j = 1:ps.ny, i = 1:ps.nx
        if flags[i, j] == 0
            if 1 in [flags[i-1, j], flags[i+1, j], flags[i, j-1], flags[i, j+1]]
                flags[i, j] = -2
            end
        end
    end

    for j = 1:ps.ny, i = 1:ps.nx
        if flags[i, j] == 1
            @assert 0 ∉ [flags[i-1, j], flags[i+1, j], flags[i, j-1], flags[i, j+1]] @show i j
        end
    end
end

ib_flag!(flags, ks.ps)

ghost_ids = findall(flags .== -2)

xbis = zeros(length(ghost_ids), 2)
nbis = zero(xbis)
for iter in axes(xbis, 1)
    idx = ghost_ids[iter]

    θ = atan(ks.ps.y[idx] / (3 - ks.ps.x[idx]))
    xbis[iter, 1] = 3 - radius * cos(θ)
    xbis[iter, 2] = radius * sin(θ)

    nbis[iter, :] .= [- radius * cos(θ), radius * sin(θ)]
end

xips = zero(xbis)
for iter in axes(xips, 1)
    idx = ghost_ids[iter]

    xips[iter, 1] = 2 * xbis[iter, 1] - ks.ps.x[idx]
    xips[iter, 2] = 2 * xbis[iter, 2] - ks.ps.y[idx]
end

ip_nids = [CartesianIndex[] for i in axes(xips, 1)]
for iter in axes(xips, 1)
    idx = Int(round(xips[iter, 1] / ps.dx[1]) + 1)
    idy = Int(round(xips[iter, 2] / ps.dy[1]) + 1)

    if flags[idx - 1, idy - 1] == 1
        push!(ip_nids[iter], CartesianIndex(idx - 1, idy - 1))
    end
    if flags[idx, idy - 1] == 1
        push!(ip_nids[iter], CartesianIndex(idx, idy - 1))
    end
    if flags[idx, idy] == 1
        push!(ip_nids[iter], CartesianIndex(idx, idy))
    end
    if flags[idx - 1, idy] == 1
        push!(ip_nids[iter], CartesianIndex(idx - 1, idy))
    end
end

ip_bids = [Int[] for i in axes(xips, 1)]
for iter in axes(xips, 1)
    if length(ip_nids[iter]) == 3
        push!(ip_bids[iter], iter)
    elseif length(ip_nids[iter]) == 2
        push!(ip_bids[iter], iter)
        if iter + 1 <= size(xips, 1)
            push!(ip_bids[iter], iter + 1)
        else
            push!(ip_bids[iter], iter - 1)
        end
    elseif length(ip_nids[iter]) == 1
        @warn "no enough interpolation points"
    end
end

function interp_coeffs(ps, xbis, nbis, nids, bids, W0)
    M = Matrix{Float64}[]
    for id in nids
        xf = ps.x[id]
        yf = ps.y[id]
        push!(M, [1 xf yf xf * yf])
    end

    if length(W0) > length(nids)
        for id in bids
            xb = xbis[id, 1]
            yb = xbis[id, 2]
            push!(M, [1 xb yb xb * yb])
        end
        W = W0
    else
        for id in bids
            nx = nbis[id, 1]
            ny = nbis[id, 2]
            xb = xbis[id, 1]
            yb = xbis[id, 2]
            push!(M, [0 nx ny xb * ny + yb * nx])
        end
        W = [W0; zeros(length(bids))]
    end
    M = vcat(M...)
    
    return M \ W
end



function update_field!(KS, ctr, a1face, a2face, flags, residual)
    nx, ny, dx, dy = KS.ps.nx, KS.ps.ny, KS.ps.dx, KS.ps.dy

    sumRes = zero(ctr[1].w)
    sumAvg = zero(ctr[1].w)

    @inbounds @threads for j ∈ 1:ny
        for i ∈ 1:nx
            if flags[i, j] == 1
                KB.step!(
                    ctr[i, j].w,
                    ctr[i, j].prim,
                    a1face[i, j].fw,
                    a1face[i+1, j].fw,
                    a2face[i, j].fw,
                    a2face[i, j+1].fw,
                    KS.gas.γ,
                    KS.ps.dx[i, j] * KS.ps.dy[i, j],
                    sumRes,
                    sumAvg,
                    :none,
                )
            end
        end
    end

    for i in eachindex(residual)
        residual[i] = sqrt(sumRes[i] * nx * ny) / (sumAvg[i] + 1.e-7)
    end

    KitBase.bc_extra!(ctr; dirc = :xr)
    KitBase.bc_extra!(ctr; dirc = :yr)
    KitBase.bc_mirror!(ctr; dirc = :yl)

    return nothing
end

function update_ghost!(ps, ctr, xbis, nbis, ip_nids, ip_bids, ghost_ids)
    for iter in eachindex(ip_nids)
        nid = ip_nids[iter]
        bid = ip_bids[iter]
        
        xf = xips[iter, 1]
        yf = xips[iter, 2]
        pos = [1, xf, yf, xf * yf]

        # U
        w1 = [ctr[idx].prim[2] for idx in nid]
        w2 = zeros(length(bid))
        w = [w1; w2]
        
        C = interp_coeffs(ps, xbis, nbis, nid, bid, w)
        U1 = C' * pos

        # V
        w1 = [ctr[idx].prim[3] for idx in nid]
        w2 = zeros(length(bid))
        w = [w1; w2]
        
        C = interp_coeffs(ps, xbis, nbis, nid, bid, w)
        V1 = C' * pos

        # T
        w1 = [1/ctr[idx].prim[4] for idx in nid]
        w2 = ones(length(bid))
        w = [w1; w2]
        
        C = interp_coeffs(ps, xbis, nbis, nid, bid, w)
        T1 = C' * pos

        # p
        w1 = [0.5 * ctr[idx].prim[1] / ctr[idx].prim[4] for idx in nid]
        w = w1
        
        C = interp_coeffs(ps, xbis, nbis, nid, bid, w)
        P1 = C' * pos

        ρ1 = 2*P1/T1

        T0 = 2 - T1
        ρ0 = ρ1 * T1 / T0


        idx = ghost_ids[iter]
        ctr[idx].prim .= [ρ0, -U1, -V1, 1/T0]
        ctr[idx].w .= prim_conserve(ctr[idx].prim, ks.gas.γ)
    end
end


#update_ghost!(ks.ps, ctr, xbis, nbis, ip_nids, ip_bids, ghost_ids)



t = 0.0
dt = timestep(ks, ctr, 0.0)
nt = ks.set.maxTime ÷ dt |> Int
res = zeros(4)


@showprogress for iter = 1:50#nt
    evolve!(ks, ctr, a1face, a2face, dt)
    update_ghost!(ks.ps, ctr, xbis, nbis, ip_nids, ip_bids, ghost_ids)
    update_field!(ks, ctr, a1face, a2face, flags, res)

    global t += dt
end


plot(ks, ctr)








###

C = interp_coeffs(ps, xbis, nbis, ip_nids[4], ip_bids[4], [0.0, 2, 1, 0.8])

C' * [1, 2, 0.2, 0.4]
