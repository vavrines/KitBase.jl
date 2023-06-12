*(Δ::AbstractDerivativeOperator, B::Nothing) = Δ

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

# test
u = Float64[0, 1, 2, 3, 2, 1, 0]
dx = 1.0
KitBase.finite_difference(u, dx; method = :central, bc = :period)
KitBase.finite_difference(u, dx; method = :central, bc = :none)
KitBase.finite_difference(u, dx, -1; method = :upwind, bc = :none)
