"""
$(SIGNATURES)

Write solution data

## Arguments
* `mode`: data format (`:bson`, `:jld`, `:vtk`, `:tec`)
"""
function write_sol(args...; mode=:bson)
    fn = eval(Symbol("write_" * string(mode)))
    fn(args...)

    return nothing
end

"""
$(SIGNATURES)

Write data into BSON
"""
function write_bson(KS::AbstractSolverSet, ctr, t=0)
    strIter = string(t)
    fileOut = KS.outputFolder * "data/t=" * strIter * ".bson"
    save(fileOut, Dict("set" => KS, "ctr" => ctr, "t" => t))

    return nothing
end

"""
$(SIGNATURES)

Write data into JLD2
"""
function write_jld(KS::AbstractSolverSet, ctr, t=0)
    strIter = string(t)
    fileOut = KS.outputFolder * "data/t=" * strIter * ".jld2"
    save(fileOut, Dict("set" => KS, "ctr" => ctr, "t" => t))

    return nothing
end

"""
$(SIGNATURES)

Write data into VTK
"""
function write_vtk(points::T, cells, cdata, pdata=zeros(axes(points, 1))) where {T<:AM}
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

"""
$(SIGNATURES)

Write data into tecplot data file
"""
function write_tec(x::AV, sol)
    open("sol.dat", "w") do f
        nx = length(x)
        len, type = begin
            if sol isa AM
                size(sol, 2), "matrix"
            elseif length(sol) == nx
                1, "scalar"
            else
                length(sol), "vector"
            end
        end

        varnm = ""
        for i in 1:len-1
            varnm = varnm * "V" * string(i) * ", "
        end
        varnm = varnm * "V" * string(len)

        println(f, "VARIABLES = X, " * varnm)
        println(f, "ZONE I = $nx")

        dp = "DATAPACKING = BLOCK"
        println(f, dp)

        write_num(f, x)

        if type == "scalar"
            write_num(f, sol)
        elseif type == "vector"
            for sv in sol
                write_num(f, sv)
            end
        elseif type == "matrix"
            for i in axes(sol)[end]
                sv = @view sol[:, i]
                write_num(f, sv)
            end
        end
    end
end

"""
$(SIGNATURES)

Write data into tecplot data file
"""
function write_tec(x::AM, y::AM, sol)
    open("sol.dat", "w") do f
        nx = size(x, 1)
        ny = size(y, 2)
        len, type = begin
            if ndims(sol) == ndims(x)
                1, "scalar"
            elseif sol isa AbstractVector
                length(sol), "vector"
            else
                size(sol)[end], "matrix"
            end
        end
        isCenter = begin
            if size(sol) == size(x) || size(sol)[1:end-1] == size(x)
                true
            else
                false
            end
        end

        varnm = ""
        for i in 1:len-1
            varnm = varnm * "V" * string(i) * ", "
        end
        varnm = varnm * "V" * string(len)

        println(f, "VARIABLES = X, Y, " * varnm)
        println(f, "ZONE I = $nx, J = $ny")

        dp = "DATAPACKING = BLOCK, "
        vl = begin
            if isCenter
                ""
            elseif len == 1
                "VARLOCATION = ([3]=CELLCENTERED)"
            else
                "VARLOCATION = ([3-$(len+2)]=CELLCENTERED)"
            end
        end
        println(f, dp * vl)

        write_num(f, x)
        write_num(f, y)

        if type == "scalar"
            write_num(f, sol)
        elseif type == "vector"
            for sv in sol
                write_num(f, sv)
            end
        elseif type == "matrix"
            for i in axes(sol)[end]
                sv = @view sol[:, :, i]
                write_num(f, sv)
            end
        end
    end

    return nothing
end

"""
$(SIGNATURES)
"""
function write_num(io, xs::AA)
    for x in xs
        @printf(io, "%.6f ", x)
    end
    print(io, "\n")

    return nothing
end
