abstract type AbstractImmersedBoundary end

"""
$(TYPEDEF)

Structure of sharp-interface immersed boundary

## Definition
- `flags`: type of cell -> 1: fluid; 0: solid; -1: outer; -2: ghost
- `idg`: index of ghost cells
- `xb`: location of body intercept points
- `nb`: normal vectors of body intercept points
- `xi`: location of image points
- `idic`: index of cell in which image points fall
- `idin`: index of face intersection of interpolation stencil
- `idib`: index of body intercept points in interpolation stencil
"""
struct SharpIB{TI,TC,TF,T1,T2,T3} <: AbstractImmersedBoundary
    flags::TI
    idg::TC
    xb::TF
    nb::TF
    xi::TF
    idic::T1
    idin::T2
    idib::T3
end

"""
$(SIGNATURES)

Label ghost cells
"""
function ghost_flag!(ps::AbstractPhysicalSpace2D, flags)
    for j = 1:ps.ny
        for i = 1:ps.nx
            if flags[i, j] == 0
                if 1 in [flags[i-1, j], flags[i+1, j], flags[i, j-1], flags[i, j+1]]
                    flags[i, j] = -2
                end
            end
        end
    end

    @threads for j = 1:ps.ny
        for i = 1:ps.nx
            if flags[i, j] == 1
                @assert 0 ∉ [flags[i-1, j], flags[i+1, j], flags[i, j-1], flags[i, j+1]] @show i j
            end
        end
    end
end

"""
$(SIGNATURES)

Compute location of image points
"""
function ip_location(ps::AbstractPhysicalSpace2D, gids, xbis)
    xips = [Vector{Float64}(undef, 2) for iter = 1:length(xbis)]
    for iter in eachindex(xips)
        idx = gids[iter]
        xips[iter][1] = 2 * xbis[iter][1] - ps.x[idx]
        xips[iter][2] = 2 * xbis[iter][2] - ps.y[idx]
    end

    return xips
end

"""
$(SIGNATURES)

Compute connectivity information of image points
"""
function ip_connectivity(ps::AbstractPhysicalSpace2D, xips, flags)
    ip_cids = Vector{CartesianIndex}(undef, 28)
    ip_nids = [CartesianIndex[] for i in eachindex(xips)]
    for iter in eachindex(xips)
        x, y = xips[iter]

        # id of the cell where 
        dxs = abs.(x .- ps.x[:, 1])
        dys = abs.(y .- ps.y[1, :])
        cidx = argmin(dxs)
        cidy = argmin(dys)
        ip_cids[iter] = CartesianIndex(cidx, cidy)

        # id of the center face intersection of the interpolation stencil
        idx = begin
            if x > ps.x[cidx, 1]
                cidx + 1
            else
                cidx
            end
        end
        idy = begin
            if y > ps.y[1, cidy]
                cidy + 1
            else
                cidy
            end
        end

        # id of non-immersed neighboring cell
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

    return ip_cids, ip_nids, ip_bids
end

"""
$(SIGNATURES)

Compute coefficients for bilinear interpolation
"""
function bilinear_coeffs(ps::AbstractPhysicalSpace2D, xbis, nbis, nids, bids, W0)
    M = Matrix{Float64}[]
    for id in nids
        xf = ps.x[id]
        yf = ps.y[id]
        push!(M, [1 xf yf xf * yf])
    end

    if length(W0) > length(nids)
        for id in bids
            xb, yb = xbis[id]
            push!(M, [1 xb yb xb * yb])
        end
        W = W0
    else
        for id in bids
            nx, ny = nbis[id]
            xb, yb = xbis[id]
            push!(M, [0 nx ny xb * ny + yb * nx])
        end
        W = [W0; zeros(length(bids))]
    end
    M = vcat(M...)
    
    return M \ W
end

function bilinear_coeffs(ps::AbstractPhysicalSpace2D, ib::SharpIB, idx, W0)
    return bilinear_coeffs(ps, ib.xb, ib.nb, ib.idin[idx], ib.idib[idx], W0)
end
