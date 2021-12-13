# ------------------------------------------------------------
# Filter functions
# ------------------------------------------------------------

"""
    modal_filter!(u, args...; filter)

Filter of modal solutions

- @arg u: 1D modal solution
- @arg args...: filter parameters including strength, norm, etc.
- @arg filter: symbolic filter options including :l2, l2opt, :l1, :lasso, :exp, :houli
"""
function modal_filter!(u::AA{T}, args...; filter::Symbol) where {T<:FN}
    filtstr = "filter_" * string(filter) * "!"
    filtfunc = Symbol(filtstr) |> eval
    filtfunc(u, args...)

    return nothing
end

function filter_l2!(u::AV{T}, args...) where {T<:FN}
    q0 = eachindex(u) |> first
    q1 = eachindex(u) |> last
    @assert q0 >= 0

    λ = args[1]
    for i = q0+1:q1
        u[i] /= (1.0 + λ * (i - q0 + 1)^2 * (i - q0)^2)
    end

    return nothing
end

function filter_l2!(u::AM{T}, args...) where {T<:FN}
    p0 = axes(u, 1) |> first
    q0 = axes(u, 2) |> first
    @assert p0 >= 0
    @assert q0 >= 0

    λx, λξ = args[1:2]
    for j in axes(u, 2), i in axes(u, 1)
        u[i, j] /=
            (1.0 + λx * (i - p0 + 1)^2 * (i - p0)^2 + λξ * (j - q0 + 1)^2 * (j - q0)^2)
    end

    return nothing
end

function filter_l2opt!(u::AV{T}, args...) where {T<:FN}
    q0 = eachindex(u) |> first
    q1 = eachindex(u) |> last
    @assert q0 >= 0

    λ = args[1]
    η = λ * 2.0
    for i = q0+1:q1
        u[i] /= (1.0 + λ * (i - q0 + 1)^2 * (i - q0)^2 - η)
    end

    return nothing
end

function filter_l2opt!(u::AM{T}, args...) where {T<:FN}
    p0 = axes(u, 1) |> first
    q0 = axes(u, 2) |> first
    @assert p0 >= 0
    @assert q0 >= 0

    λx, λξ = args[1:2]
    η0 = λx * 2.0^2 + λξ * 2.0^2
    for j in axes(u, 2)
        for i in axes(u, 1)
            if i == p0 && j == q0
                continue
            elseif i == 1
                η = λξ * 2.0^2
            elseif j == 1
                η = λx * 2.0^2
            else
                η = η0
            end

            u[i, j] /= (
                1.0 + λx * (i - p0 + 1)^2 * (i - p0)^2 + λξ * (j - q0 + 1)^2 * (j - q0)^2 - η
            )
        end
    end

    return nothing
end

function filter_l1!(u::AV{T}, args...) where {T<:FN}
    q0 = eachindex(u) |> first
    q1 = eachindex(u) |> last
    @assert q0 >= 0

    λ = args[1]
    ℓ = args[2]
    for i = q0+1:q1
        sc = 1.0 - λ * i * (i - 1) * ℓ[i] / (abs(u[i]) + 1e-8)
        if sc < 0.0
            sc = 0.0
        end
        u[i] *= sc
    end

    return nothing
end

function filter_l1!(u::AM{T}, args...) where {T<:FN}
    p0 = axes(u, 1) |> first
    q0 = axes(u, 2) |> first
    @assert p0 >= 0
    @assert q0 >= 0

    λ1, λ2 = args[1:2]
    ℓ = args[3]
    for j in axes(u, 2), i in axes(u, 1)
        sc = 1.0 - (λ1 * i * (i - 1) + λ2 * j * (j - 1)) * ℓ[i, j] / (abs(u[i, j]) + 1e-8)
        if sc < 0.0
            sc = 0.0
        end
        u[i, j] *= sc
    end

    return nothing
end

function filter_lasso!(u::AV{T}, args...) where {T<:FN}
    q0 = eachindex(u) |> first
    @assert q0 >= 0

    ℓ = args[1]
    nr = length(u)
    λ = abs(u[end]) / (nr * (nr - 1) * ℓ[end])
    filter_l1!(u, λ, ℓ)

    return nothing
end

function filter_lasso!(u::AM{T}, args...) where {T<:FN}
    nr, nz = size(u)
    ℓ = args[1]
    λ1 = abs(u[end, 1]) / (nr * (nr - 1) * ℓ[end, 1])
    λ2 = abs(u[1, end]) / (nz * (nz - 1) * ℓ[1, end])

    filter_l1!(u, λ1, λ2, ℓ)

    return nothing
end

function filter_exp!(u::AV{T}, args...) where {T<:FN}
    N = length(u) - 1
    s = args[1]
    Nc = begin
        if length(args) > 1
            args[2]
        else
            0
        end
    end

    σ = filter_exp1d(N, s, Nc)
    u .*= σ

    return nothing
end

function filter_houli!(u::AV{T}, args...) where {T<:FN}
    N = length(u) - 1
    s = args[1]
    Nc = begin
        if length(args) > 1
            args[2]
        else
            0
        end
    end

    σ = filter_exp1d(N, s, Nc)
    for i in eachindex(σ)
        if i / length(σ) < 2 / 3
            σ[i] = 1.0
        end
    end

    u .*= σ

    return nothing
end


"""
    filter_exp1d(N, s, Nc)

Construct exponential filter for modal solution

- @arg N: degree of polynomials
- @arg s: order of filter (must be even)
- @arg Nc: cutoff location
"""
function filter_exp1d(N, s, Nc = 0)
    alpha = -log(eps())

    filterdiag = ones(N + 1)
    for i = Nc:N
        filterdiag[i+1] = exp(-alpha * ((i - Nc) / (N - Nc))^s)
    end

    return filterdiag
end
