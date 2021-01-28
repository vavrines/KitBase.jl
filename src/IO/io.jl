# ============================================================
# I / O Methods
# ============================================================

export read_dict,
       write_jld,
       plot_line,
       plot_contour

"""
    read_dict(filename::String, allowed)
    read_dict(filename::String)

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

        var, val = split(line, "=")
        stripped = strip(var)
        stripped = Symbol(stripped)
        println(line)

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
function write_jld(
    KS::X,
    ctr::Y,
    t = 0,
) where {X<:AbstractSolverSet,Y}

    strIter = string(t)
    fileOut = KS.outputFolder * "data/t=" * strIter * ".jld2"

    save(fileOut, Dict("set" => KS, "ctr" => ctr, "t" => t))

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

Plot solution contour
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

    p1 = contourf(KS.pSpace.x[1:KS.pSpace.nx, 1], KS.pSpace.y[1, 1:KS.pSpace.ny], sol[1, :, :]')
    p2 = contourf(KS.pSpace.x[1:KS.pSpace.nx, 1], KS.pSpace.y[1, 1:KS.pSpace.ny], sol[2, :, :]')
    p3 = contourf(KS.pSpace.x[1:KS.pSpace.nx, 1], KS.pSpace.y[1, 1:KS.pSpace.ny], sol[3, :, :]')
    p4 = contourf(KS.pSpace.x[1:KS.pSpace.nx, 1], KS.pSpace.y[1, 1:KS.pSpace.ny], sol[4, :, :]')

    plot(p1, p2, p3, p4, layout = (2, 2))

end