symbolize(x::AbstractString) = Symbol(x)
symbolize(x::AbstractArray{T}) where {T<:AbstractString} = [symbolize(y) for y in x]


function static_array(x::AbstractVector)
    y = MVector{length(x)}(collect(x))

    if x isa OffsetArray
        idx0 = eachindex(x) |> first
        idx1 = eachindex(x) |> last

        y = OffsetArray(y, idx0:idx1)
    end

    return y
end

function static_array(x::AbstractMatrix)
    y = MMatrix{size(x, 1),size(x, 2)}(collect(x))

    if x isa OffsetArray
        idx0 = axes(x, 1) |> first
        idx1 = axes(x, 1) |> last
        idy0 = axes(x, 2) |> first
        idy1 = axes(x, 2) |> last

        y = OffsetArray(y, idx0:idx1, idy0:idy1)
    end

    return y
end

function static_array(x::AbstractArray{<:Number,3})
    y = MArray{Tuple{size(x, 1),size(x, 2),size(x, 3)}}(collect(x))

    if x isa OffsetArray
        idx0 = axes(x, 1) |> first
        idx1 = axes(x, 1) |> last
        idy0 = axes(x, 2) |> first
        idy1 = axes(x, 2) |> last
        idz0 = axes(x, 3) |> first
        idz1 = axes(x, 3) |> last

        y = OffsetArray(y, idx0:idx1, idy0:idy1, idz0:idz1)
    end

    return y
end

function static_array(x::AbstractArray{<:Number,4})
    y = MArray{Tuple{size(x, 1),size(x, 2),size(x, 3),size(x, 4)}}(collect(x))

    if x isa OffsetArray
        ida0 = axes(x, 1) |> first
        ida1 = axes(x, 1) |> last
        idb0 = axes(x, 2) |> first
        idb1 = axes(x, 2) |> last
        idc0 = axes(x, 3) |> first
        idc1 = axes(x, 3) |> last
        idd0 = axes(x, 4) |> first
        idd1 = axes(x, 4) |> last

        y = OffsetArray(y, ida0:ida1, idb0:idb1, idc0:idc1, idd0:idd1)
    end

    return y
end


dynamic_array(x::AbstractArray) = x
dynamic_array(x::Number) = x
