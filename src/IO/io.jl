# ============================================================
# I / O Methods
# ============================================================

export read_dict,
       write_jld,
       write_vtk,
       plot_line,
       plot_contour

"""
    read_dict(filename::T, allowed) where {T<:AbstractString}
    read_dict(filename::T) where {T<:AbstractString}

Read text into dictionary

* @args filename: configuration text file
* @args allowed: keywords in range
* @return vars: dictionary with values of variables

"""
function read_dict(filename::T, allowed) where {T<:AbstractString}

    @info "Reading config from $filename"
    println("--------------------------------------------------------------")
    f = open(filename)
    vars = Dict{String,Any}()

    for line in eachline(f)
        # skip comments
        if length(line) == 0 || line[1] == '#'
            continue
        end

        #print("\t");
        #println(line)
        var, val = split(line, "=")
        stripped = strip(var)
        if stripped in allowed
            println(line)

            #vars[stripped] = parse(Float64, val)
            #vars[stripped] = strip(val)
            tmp = tryparse(Float64, val)
            if isa(tmp, Nothing)
                vars[stripped] = strip(val)
            else
                vars[stripped] = isinteger(tmp) ? Int(tmp) : tmp
            end
        end
    end
    println("--------------------------------------------------------------")
    println("")
    return vars

end

function read_dict(filename::T) where {T<:AbstractString}

    @info "Reading config from $filename"
    println("--------------------------------------------------------------")
    f = open(filename)
    vars = Dict{Symbol,Any}()

    for line in eachline(f)
        if length(line) == 0 || line[1] == '#'
            continue
        end
        println(line)

        var, val = split(line, "=")
        stripped = strip(var)
        stripped = Symbol(stripped)
        val = split(val, "#")[1]

        tmp = tryparse(Float64, val)
        if isa(tmp, Nothing)
            vars[stripped] = strip(val)
        else
            vars[stripped] = isinteger(tmp) ? Int(tmp) : tmp
        end
    end
    println("--------------------------------------------------------------")
    println("")
    return vars

end


"""
    write_jld(KS, ctr, t)

Write data into JLD2
"""
function write_jld(KS::X, ctr::Y, t = 0) where {X<:AbstractSolverSet,Y}

    strIter = string(t)
    fileOut = KS.outputFolder * "data/t=" * strIter * ".jld2"

    save(fileOut, Dict("set" => KS, "ctr" => ctr, "t" => t))

end


"""
    write_vtk(points, cells, cdata, pdata = zeros(axes(points, 1)))

Write data into VTK
"""
function write_vtk(
    points::T,
    cells,
    cdata,
    pdata = zeros(axes(points, 1)),
) where {T<:AbstractMatrix}
    mcells = [MeshCell(VTKCellTypes.VTK_TRIANGLE, cells[i, :]) for i in axes(cells, 1)]
    vtkfile = vtk_grid("sol", permutedims(points), mcells)

    len = size(cdata, 2)
    cname = Array{String}(undef, len)
    for i in eachindex(cname)
        cname[i] = "c" * string(i)
        vtkfile[cname[i], VTKCellData()] = cdata[:, i]
    end

    len = size(pdata, 2)
    pname = Array{String}(undef, len)
    for i in eachindex(pname)
        pname[i] = "p" * string(i)
        vtkfile[pname[i], VTKPointData()] = pdata[:, i]
    end

    outfiles = vtk_save(vtkfile)

    return nothing
end

function write_vtk(ks::T1, ctr) where {T1<:AbstractSolverSet}
    cdata = zeros(length(ctr), length(ctr[1].w))
    for i in eachindex(ctr)
        cdata[i, :] .= ctr[i].prim
        cdata[i, end] = 1.0 / cdata[i, end]
    end
    KitBase.write_vtk(ks.pSpace.points, ks.pSpace.cellid, cdata)

    return nothing
end


"""
    plot_line(KS, ctr; backend = :plots)

Plot solution profiles
"""
function plot_line(
    KS::X,
    ctr::Y;
    backend = :plots::Symbol,
) where {X<:AbstractSolverSet,Y<:AbstractArray{<:AbstractControlVolume,1}}

    pltx = KS.pSpace.x[1:KS.pSpace.nx]
    plty = zeros(KS.pSpace.nx, 6)

    for i in eachindex(pltx)
        for j = 1:2
            plty[i, j] = ctr[i].prim[j]
        end

        plty[i, 3] = 1.0 / ctr[i].prim[end]
    end

    if backend == :plots
        p1 = plot(pltx, plty[:, 1], label = "Density", lw = 2, xlabel = "x")
        p1 = plot!(pltx, plty[:, 2], label = "Velocity", lw = 2)
        p1 = plot!(pltx, plty[:, 3], label = "Temperature", lw = 2)
        display(p1)
    elseif backend == :gr
        xlabel("x")
        ylabel("Density")
        legend("n")
        p1 = plot(pltx, plty[:, 1])
        display(p1)
        xlabel("x")
        ylabel("Velocity")
        legend("U")
        p2 = plot(pltx, plty[:, 2])
        display(p2)
        xlabel("x")
        ylabel("Temperature")
        legend("T")
        p3 = plot(pltx, plty[:, 3])
        display(p3)
    else
        throw("undefined plotting backend")
    end

end


"""
    plot_contour(KS, ctr; backend = :plots)

Plot solution contours
"""
function plot_contour(
    KS::X,
    ctr::Y;
    backend = :plots::Symbol,
) where {X<:AbstractSolverSet,Y<:AbstractArray{<:AbstractControlVolume,2}}

    sol = zeros(size(ctr[1].w, 1), KS.pSpace.nx, KS.pSpace.ny)
    for i in axes(sol, 2)
        for j in axes(sol, 3)
            for k = 1:size(sol, 1)-1
                sol[k, i, j] = ctr[i, j].prim[k]
            end
            sol[end, i, j] = 1.0 / ctr[i, j].prim[end]
        end
    end

    p1 = contourf(
        KS.pSpace.x[1:KS.pSpace.nx, 1],
        KS.pSpace.y[1, 1:KS.pSpace.ny],
        sol[1, :, :]',
    )
    p2 = contourf(
        KS.pSpace.x[1:KS.pSpace.nx, 1],
        KS.pSpace.y[1, 1:KS.pSpace.ny],
        sol[2, :, :]',
    )
    p3 = contourf(
        KS.pSpace.x[1:KS.pSpace.nx, 1],
        KS.pSpace.y[1, 1:KS.pSpace.ny],
        sol[3, :, :]',
    )
    p4 = contourf(
        KS.pSpace.x[1:KS.pSpace.nx, 1],
        KS.pSpace.y[1, 1:KS.pSpace.ny],
        sol[4, :, :]',
    )

    plot(p1, p2, p3, p4, layout = (2, 2))

end
