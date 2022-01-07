const difffn = Dict(
    :Central => :CenteredDifference,
    :Center => :CenteredDifference,
    :Centered => :CenteredDifference,
    :central => :CenteredDifference,
    :center => :CenteredDifference,
    :centered => :CenteredDifference,
    :CenteredDifference => :CenteredDifference,
    :Upwind => :UpwindDifference,
    :upwind => :UpwindDifference,
    :UpwindDifference => :UpwindDifference,
)

const bcfn = Dict(
    :Period => :PeriodicBC,
    :Periodic => :PeriodicBC,
    :period => :PeriodicBC,
    :Period => :PeriodicBC,
    :PeriodBC => :PeriodicBC,
    :Dirichlet => :DirichletBC,
    :dirichlet => :DirichletBC,
    :DirichletBC => :DirichletBC,
    :none => nothing,
)


*(Δ::AbstractDerivativeOperator, B::Nothing) = Δ


"""
$(SIGNATURES)

Finite difference operation

Note that for upwind difference, the coeff = 1/-1 can't take the default value as nothing.
"""
function finite_difference(
    y::AA,
    x::AA,
    coeff = 1,
    uL = first(y),
    uR = last(y);
    method = :central,
    dimension = 1,
    derivative = 1,
    truncation = 2,
    bc = :dirichlet,
)

    oprtfn = eval(difffn[method])

    i0 = firstindex(x)
    i1 = lastindex(x)

    if bc in (nothing, :none)
        dx = x[i0+1:i1] .- x[i0:i1-1]
        Δ = oprtfn{dimension}(derivative, truncation, dx, length(y) - 2, coeff)
        B = nothing
    else
        dx = [x[i0+1] - x[i0]; x[i0+1:i1] .- x[i0:i1-1]; x[i1] - x[i1-1]]
        Δ = oprtfn{dimension}(derivative, truncation, dx, length(y), coeff)
    end

    bcnm = bcfn[bc]
    bcrtfn = eval(bcnm)

    if bcnm == :PeriodicBC
        B = bcrtfn(eltype(y))
    elseif bcnm == :DirichletBC
        B = bcrtfn(uL, uR)
    end

    return Δ * B * y

end

function finite_difference(y, dx::Real, args...; kwargs...)
    x = linspace(0, length(y) - 1, length(y)) .* dx
    return finite_difference(y, x, args...; kwargs...)
end


"""
$(SIGNATURES)

Central difference
"""
function central_diff(y::T, x::T) where {T<:AV{<:Number}}
    dy = zeros(eltype(y), axes(y))

    idx = eachindex(y) |> collect
    i0 = idx[1]
    i1 = idx[end]

    dy[i0] = (y[i0+1] - y[i0]) / (x[i0+1] - x[i0] + 1e-7)
    dy[i1] = (y[i1] - y[i1-1]) / (x[i1] - x[i1-1] + 1e-7)
    for i = i0+1:i1-1
        dy[i] = (y[i+1] - y[i-1]) / (x[i+1] - x[i-1] + 1e-7)
    end

    return dy
end

"""
$(SIGNATURES)
"""
function central_diff(y::AV{T}, dx) where {T}
    x = linspace(0, length(y) - 1, length(y)) .* dx
    dy = central_diff(y, x)

    return dy
end


"""
$(SIGNATURES)

Central difference for second-order derivatives
"""
function central_diff2(y::T, x::T) where {T<:AV{<:Number}}
    dy = zeros(eltype(y), axes(y))

    i0 = eachindex(y) |> first
    i1 = eachindex(y) |> last

    for i = i0+1:i1-1
        dy[i] = (y[i+1] - 2.0 * y[i] + y[i-1]) / (x[i+1] - x[i-1])^2
    end

    return dy
end

"""
$(SIGNATURES)
"""
function central_diff2(y::AV{T}, dx) where {T}
    x = linspace(0, length(y) - 1, length(y)) .* dx
    dy = central_diff2(y, x)

    return dy
end


"""
$(SIGNATURES)

Central difference
"""
function central_diff!(dy::AV{T1}, y::T2, x::T2) where {T1,T2<:AV{<:Number}}
    @assert axes(dy) == axes(y) == axes(x)

    idx = eachindex(y) |> collect
    i0 = idx[1]
    i1 = idx[end]

    dy[i0] = (y[i0+1] - y[i0]) / (x[i0+1] - x[i0] + 1e-7)
    dy[i1] = (y[i1] - y[i1-1]) / (x[i1] - x[i1-1] + 1e-7)
    for i = i0+1:i1-1
        dy[i] = (y[i+1] - y[i-1]) / (x[i+1] - x[i-1] + 1e-7)
    end

    return nothing
end

"""
$(SIGNATURES)
"""
function central_diff!(dy::AV{T1}, y::AV{T2}, dx) where {T1,T2}
    x = linspace(0, length(y) - 1, length(y)) .* dx
    central_diff!(dy, y, x)

    return nothing
end


"""
$(SIGNATURES)
    
Central difference
"""
function central_diff2!(dy::AV{T1}, y::T2, x::T2) where {T1,T2<:AV{<:Number}}
    @assert axes(dy) == axes(y) == axes(x)

    i0 = eachindex(y) |> first
    i1 = eachindex(y) |> last

    dy[i0] = 0.0
    dy[i1] = 0.0
    for i = i0+1:i1-1
        dy[i] = (y[i+1] - 2.0 * y[i] + y[i-1]) / (x[i+1] - x[i-1])^2
    end

    return nothing
end

"""
$(SIGNATURES)
"""
function central_diff2!(dy::AV{T1}, y::AV{T2}, dx) where {T1,T2}
    x = ones(eltype(y), axes(y)) .* dx
    central_diff2!(dy, y, x)

    return nothing
end


"""
$(SIGNATURES)

Upwind difference
"""
function upwind_diff(y::AV{T1}, x::AV{T2}; stream = :right::Symbol) where {T1,T2}
    dy = zeros(eltype(y), axes(y))

    idx = eachindex(y) |> collect
    i0 = idx[1]
    i1 = idx[end]

    if stream == :right
        dy[i0] = 0.0
        for i = i0+1:i1
            dy[i] = (y[i] - y[i-1]) / (x[i] - x[i-1] + 1e-7)
        end
    elseif stream == :left
        dy[i1] = 0.0
        for i = i0:i1-1
            dy[i] = (y[i+1] - y[i]) / (x[i+1] - x[i] + 1e-7)
        end
    else
        throw("streaming direction should be :left or :right")
    end

    return dy
end

"""
$(SIGNATURES)
"""
function upwind_diff(y::AV{T}, dx; stream = :right::Symbol) where {T}
    x = linspace(0, length(y) - 1, length(y)) .* dx
    dy = upwind_diff(y, x, stream = stream)

    return dy
end


"""
$(SIGNATURES)

Upwind difference
"""
function upwind_diff!(
    dy::AV{T1},
    y::AV{T2},
    x::AV{T3};
    stream = :right::Symbol,
) where {T1,T2,T3}

    @assert axes(dy) == axes(y) == axes(x)

    idx = eachindex(y) |> collect
    i0 = idx[1]
    i1 = idx[end]

    if stream == :right
        dy[i0] = 0.0
        for i = i0+1:i1
            dy[i] = (y[i] - y[i-1]) / (x[i] - x[i-1] + 1e-7)
        end
    elseif stream == :left
        dy[i1] = 0.0
        for i = i0:i1-1
            dy[i] = (y[i+1] - y[i]) / (x[i+1] - x[i] + 1e-7)
        end
    else
        throw("streaming direction should be :left or :right")
    end

    return nothing

end

"""
$(SIGNATURES)
"""
function upwind_diff!(dy::AV{T1}, y::AV{T2}, dx; stream = :right::Symbol) where {T1,T2}
    x = linspace(0, length(y) - 1, length(y)) .* dx
    upwind_diff!(dy, y, x, stream = stream)

    return nothing
end


"""
$(SIGNATURES)

Finite difference for pseudo-unstructured mesh
"""
function unstruct_diff(
    u::AV{T1},
    x::AV{T2},
    nx::Integer;
    mode = :central::Symbol,
) where {T1,T2}

    uu = reshape(u, (nx, :))
    xx = reshape(x, (nx, :))

    dux = similar(xx)
    for i = 1:nx
        if mode == :central
            dux[i, :] .= central_diff(uu[i, :], xx[i, :])
        elseif mode == :upwind
            dux[i, :] .= upwind_diff(uu[i, :], xx[i, :])
        else
            throw("difference mode should be central or upwind")
        end
    end

    return reshape(dux, (1, :))

end

"""
$(SIGNATURES)
"""
function unstruct_diff(
    u::Function,
    x::AV{T},
    nx::Integer,
    dim::Integer;
    mode = :central::Symbol,
) where {T}

    uu = reshape(u.(x), (nx, :))
    xx = reshape(x, (nx, :))
    dux = zeros(eltype(x), axes(xx))

    if dim == 1
        for j in axes(xx, 2)
            if mode == :central
                dux[:, j] .= central_diff(uu[:, j], xx[:, j])
            elseif mode == :upwind
                dux[:, j] .= upwind_diff(uu[:, j], xx[:, j])
            end
        end
    elseif dim == 2
        for i in axes(xx, 1)
            if mode == :central
                dux[i, :] .= central_diff(uu[i, :], xx[i, :])
            elseif mode == :upwind
                dux[i, :] .= upwind_diff(uu[i, :], xx[i, :])
            end
        end
    end

    return dux[:]

end
