"""
$(SIGNATURES)

Create basis functions for entropy closure
"""
function moment_basis(u, n::Integer)
    m = similar(u, n, length(u))
    for j in axes(m, 2)
        for i in axes(m, 1)
            m[i, j] = u[j]^(i - 1)
        end
    end

    return m
end

"""
$(SIGNATURES)
"""
function moment_basis(u, v, n::Integer)
    m = Array{eltype(u)}(undef, 0, length(u))

    for i = 1:n
        for idu = 0:i-1
            idv = i - 1 - idu
            tmp = u.^idu .* v.^idv |> permutedims
            m = vcat(m, tmp)
        end
    end

    return m
end

"""
$(SIGNATURES)

According to Levermore's paper, the dimensions of moment basis are 10 and 
35 for n = 3 and n = 5 respectively.
"""
function moment_basis(u, v, w, n::Integer)
    m = Array{eltype(u)}(undef, 0, length(u))

    for i = 1:n
        for idu = 0:i-1
            for idv = 0:i-1-idu
                idw = i - 1 - idu - idv
                tmp = u.^idu .* v.^idv .* w.^idw |> permutedims
                m = vcat(m, tmp)
            end
        end
    end

    return m
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
"""
function sample_pdf(m, prim, pdf)
    α = zeros(size(m, 1))
    α[1] = log((prim[1] / (π / prim[end]))^0.5) - prim[2]^2 * prim[end]
    α[2] = 2.0 * prim[2] * prim[3]
    α[3] = -prim[3]

    pdf1 = Truncated(pdf, -Inf, 0)
    for i = 4:lastindex(α)
        if iseven(i)
            α[i] = rand(pdf)
        else
            α[i] = rand(pdf1)
        end
    end

    return exp.(α' * m)[:]
end
