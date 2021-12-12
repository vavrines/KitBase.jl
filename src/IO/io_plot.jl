"""
    plot_line(KS, ctr)

Plot solution profiles
"""
plot_line(args...; kwargs...) = @info "plot_line is deprecated, use Plots.plot instead."

@recipe function plot_line_backend(
    KS::X,
    ctr::Y,
) where {X<:AbstractSolverSet,Y<:AA{<:AbstractControlVolume,1}}

    # solution
    pltx = KS.pSpace.x[1:KS.pSpace.nx]
    plty = zeros(KS.pSpace.nx, 3)
    for i in eachindex(pltx)
        for j = 1:2
            plty[i, j] = ctr[i].prim[j]
        end

        plty[i, 3] = 1.0 / ctr[i].prim[end]
    end

    # attributes
    xguide --> "x"
    :linewidth --> 1.5

    @series begin
        label := "ρ"
        pltx, plty[:, 1]
    end

    if ctr[1].w isa AA
        @series begin
            label := "u"
            pltx, plty[:, 2]
        end
        @series begin
            label := "T"
            pltx, plty[:, 3]
        end
    end

    # user-defined
    c = get(plotattributes, :linewidth, :auto)

    return nothing
    
end


"""
    plot_line(KS, ctr)

Plot solution profiles
"""
plot_contour(args...; kwargs...) = @info "plot_line is deprecated, use Plots.plot instead."

@recipe function plot_contour_backend(
    KS::X,
    ctr::Y,
) where {X<:AbstractSolverSet,Y<:AA{<:AbstractControlVolume,2}}

    pltx = KS.ps.x[1:KS.ps.nx, 1]
    plty = KS.ps.y[1, 1:KS.ps.ny]

    sol = zeros(size(ctr[1].w, 1), KS.pSpace.nx, KS.pSpace.ny)
    for i in axes(sol, 2)
        for j in axes(sol, 3)
            for k = 1:size(sol, 1)-1
                sol[k, i, j] = ctr[i, j].prim[k]
            end
            sol[end, i, j] = 1.0 / ctr[i, j].prim[end]
        end
    end

    layout := (2, 2)
    c = get(plotattributes, :inferno, :auto)
    for (i, l) in enumerate(("ρ", "u", "v", "T"))
        @series begin
            subplot := i
            xguide := "x"
            yguide := "y"
            fill := true
            seriescolor := c
            pltx, plty, sol[i, :, :]'
        end
    end

    return nothing

end
