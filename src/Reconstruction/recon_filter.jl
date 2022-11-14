# ------------------------------------------------------------
# Filter functions
# ------------------------------------------------------------

"""
$(SIGNATURES)

Filter of modal solutions

## Arguments
- `u`: modal solution
- `args...`: filter parameters including strength, norm, etc.
- `filter`: symbolic filter options (`:l2`, `l2opt`, `:l1`, `:lasso`, `:exp`, `:houli`)
"""
function modal_filter!(u::AA, args...; filter::Symbol)
    filtstr = "filter_" * string(filter) * "!"
    filtfunc = Symbol(filtstr) |> eval
    filtfunc(u, args...)

    return nothing
end


function filter_l2!(u::AV, λ)
    q0 = firstindex(u)
    @assert q0 >= 0

    for i in eachindex(u)
        u[i] /= (1.0 + λ * (i - q0 + 1)^2 * (i - q0)^2)
    end

    return nothing
end

function filter_l2!(u::AM, λx, λξ)
    p0 = axes(u, 1) |> first
    q0 = axes(u, 2) |> first
    @assert p0 >= 0
    @assert q0 >= 0

    for j in axes(u, 2), i in axes(u, 1)
        u[i, j] /=
            (1.0 + λx * (i - p0 + 1)^2 * (i - p0)^2 + λξ * (j - q0 + 1)^2 * (j - q0)^2)
    end

    return nothing
end


function filter_l2opt!(u::AV, λ)
    q0 = firstindex(u)
    @assert q0 >= 0

    η = λ * 2.0
    for i in eachindex(u)
        u[i] /= (1.0 + λ * (i - q0 + 1)^2 * (i - q0)^2 - η)
    end

    return nothing
end

function filter_l2opt!(u::AM, λx, λξ)
    p0 = axes(u, 1) |> first
    q0 = axes(u, 2) |> first
    @assert p0 >= 0
    @assert q0 >= 0

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


function filter_l1!(u::AV, λ, ℓ)
    q0 = firstindex(u)
    @assert q0 >= 0

    for i in eachindex(u)
        sc = 1.0 - λ * i * (i - 1) * ℓ[i] / (abs(u[i]) + 1e-8)
        if sc < 0.0
            sc = 0.0
        end
        u[i] *= sc
    end

    return nothing
end

function filter_l1!(u::AM, λ1, λ2, ℓ)
    p0 = axes(u, 1) |> first
    q0 = axes(u, 2) |> first
    @assert p0 >= 0
    @assert q0 >= 0

    for j in axes(u, 2), i in axes(u, 1)
        sc = 1.0 - (λ1 * i * (i - 1) + λ2 * j * (j - 1)) * ℓ[i, j] / (abs(u[i, j]) + 1e-8)
        if sc < 0.0
            sc = 0.0
        end
        u[i, j] *= sc
    end

    return nothing
end


function filter_lasso!(u::AV, ℓ)
    q0 = eachindex(u) |> first
    @assert q0 >= 0

    nr = length(u)
    λ = abs(u[end]) / (nr * (nr - 1) * ℓ[end])
    filter_l1!(u, λ, ℓ)

    return nothing
end

function filter_lasso!(u::AM, ℓ)
    nr, nz = size(u)
    λ1 = abs(u[end, 1]) / (nr * (nr - 1) * ℓ[end, 1])
    λ2 = abs(u[1, end]) / (nz * (nz - 1) * ℓ[1, end])
    filter_l1!(u, λ1, λ2, ℓ)

    return nothing
end


function filter_exp!(u::AV, s, Nc = 0)
    N = length(u) - 1
    σ = filter_exp1d(N, s, Nc)
    u .*= σ

    return nothing
end

function filter_exp!(u::AM, spx, spy = spx, λ = 1.0, Nc = 0)
    nx, nz = size(u)
    σ = filter_exp2d(nx - 1, nz - 1, spx, spy) .^ λ
    u .*= σ

    return nothing
end


function filter_houli!(u::AV, s, Nc = 0)
    N = length(u) - 1
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
$(SIGNATURES)

Calculate strength for 1D exponential filter

Note that the implementation here, `filterdiag[i+1] = exp(-alpha * (i / (N + 1))^s)`,
is slightly different from Hesthaven's monograph, `filterdiag[i+1] = exp(-alpha * ((i - Nc) / (N - Nc))^s)`
"""
function filter_exp1d(N, s, Nc = 0)
    alpha = -log(eps())

    filterdiag = ones(N + 1)
    for i = Nc:N
        filterdiag[i+1] = exp(-alpha * (i / (N + 1))^s)
    end

    return filterdiag
end


"""
$(SIGNATURES)

Calculate strength for 2D exponential filter
"""
function filter_exp2d(Nx, Ny, spx, spy, Nc = 0)
    alpha = -log(eps())

    filterdiag = ones(Nx + 1, Ny + 1)
    for i = 0:Nx
        for j = 0:Ny
            if i + j >= Nc
                filterdiag[i+1, j+1] =
                    exp(-alpha * ((i / (Nx + 1))^spx + (j / (Ny + 1))^spy))
            end
        end
    end

    return filterdiag
end
