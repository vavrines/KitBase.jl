"""
$(SIGNATURES)

Create basis functions for moment closure hierarchies

_C. D. Levermore. J. Stat. Phys. 83(5): 1021-1065, 1996._
"""
function moment_basis(u, n::Integer)
    m = zeros(typeof(u), n + 1)
    for i = 1:n+1
        m[i] = u^(i - 1)
    end

    return m
end

"""
$(SIGNATURES)
"""
moment_basis(u::AV, n::Integer) = hcat(moment_basis.(u, n)...)

"""
$(SIGNATURES)

2D closure
"""
function moment_basis(u, v, n::Integer)
    idx = moment_index(n, Dimension{2})
    m = zeros(length(idx))
    for i in eachindex(m)
        idu, idv = idx[i]
        m[i] = u^idu * v^idv
    end

    return m
end

"""
$(SIGNATURES)
"""
moment_basis(u::AA, v::AA, n::Integer) = hcat(moment_basis.(u, v, n)...)

"""
$(SIGNATURES)

3D closure

The dimensions of moment basis are 10 and 35 for n = 2 and n = 4 respectively.
"""
function moment_basis(u, v, w, n::Integer)
    idx = moment_index(n, Dimension{3})
    m = Array{typeof(u)}(undef, length(idx))
    for i in eachindex(m)
        idu, idv, idw = idx[i]
        m[i] = u^idu * v^idv * w^idw
    end

    return m
end

"""
$(SIGNATURES)
"""
moment_basis(u::AA, v::AA, w::AA, n::Integer) = hcat(moment_basis.(u, v, w, n)...)


"""
$(SIGNATURES)

Create index of basis function in moment closure
"""
moment_index(n, ::Type{Dimension{1}}) = collect(0:n)

"""
$(SIGNATURES)
"""
function moment_index(n, ::Type{Dimension{2}})
    idx = Tuple[]
    for i = 1:n+1
        for idu = i-1:-1:0
            idv = i - 1 - idu
            push!(idx, (idu, idv))
        end
    end

    return idx
end

"""
$(SIGNATURES)
"""
function moment_index(n, ::Type{Dimension{3}})
    idx = Tuple[]
    for i = 1:n+1
        for idu = i-1:-1:0
            for idv = i-1-idu:-1:0
                idw = i - 1 - idu - idv
                push!(idx, (idu, idv, idw))
            end
        end
    end

    return idx
end


"""
$(SIGNATURES)

Optimizer for the entropy closure: `argmin(<η(α*m)> - α*u)`
"""
function optimize_closure(α, m, ω, u, η::Function; optimizer = Newton())
    loss(_α) = sum(η.(_α' * m)[:] .* ω) - dot(_α, u)
    res = Optim.optimize(loss, α, optimizer)

    return res
end


"""
$(SIGNATURES)

Reconstruct realizable moments based on Lagrange multiplier α
"""
function realizable_reconstruct(α, m, ω, η_dual_prime::Function)
    u = zero(α)
    for i in eachindex(u)
        u[i] = sum(η_dual_prime.(α' * m)[:] .* ω .* m[i, :])
    end

    return u
end


"""
$(SIGNATURES)

Sample distribution functions from entropy closure

## Arguments
* `m`: basis function {1, u, u^2, ..., u^n}
* `n`: order
* `prim`: primitive variables
* `pdf`: probability density
"""
function sample_pdf(m, n::Integer, prim, pdf)
    pdf1 = truncated(pdf, -Inf, 0)
    α = zeros(size(m, 1))

    if length(prim) == 3
        α[1] = log((prim[1] / (π / prim[end]))^0.5) - prim[2]^2 * prim[end]
        α[2] = 2.0 * prim[2] * prim[end]
        α[3] = -prim[3]

        cut = 4
        idx = moment_index(n, Dimension{1})
    elseif length(prim) == 4
        α[1] = log((prim[1] / (π / prim[end]))) - (prim[2]^2 + prim[3]^2) * prim[end]
        α[2] = 2.0 * prim[2] * prim[end]
        α[3] = 2.0 * prim[3] * prim[end]
        α[4] = -prim[end]
        α[5] = 0.0
        α[6] = -prim[end]

        cut = 7
        idx = moment_index(n, Dimension{2})
    elseif length(prim) == 5
        α[1] =
            log((prim[1] / (π / prim[end]))^1.5) -
            (prim[2]^2 + prim[3]^2 + prim[4]^2) * prim[end]
        α[2] = 2.0 * prim[2] * prim[end]
        α[3] = 2.0 * prim[3] * prim[end]
        α[4] = 2.0 * prim[4] * prim[end]
        α[5] = -prim[end]
        α[6] = 0.0
        α[7] = 0.0
        α[8] = -prim[end]
        α[9] = 0.0
        α[10] = -prim[end]

        cut = 11
        idx = moment_index(n, Dimension{3})
    end

    for i = cut:lastindex(α)
        if iseven(sum(idx[i]))
            α[i] = rand(pdf1)
        else
            α[i] = rand(pdf)
        end
    end

    return exp.(α' * m)[:]
end
