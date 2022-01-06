"""
$(SIGNATURES)

Sample particle velocity from Maxwellian. Unconstraint methods use Box-Muller sampling,
while constraint one depends on Distributions.jl.
"""
function sample_maxwell(λ::Real, v = 0.0::Real)
    rd1 = rand()
    rd2 = rand()

    return sqrt(-log(rd1) / λ) * sin(2.0 * π * rd2) + v
end

"""
$(SIGNATURES)
"""
function sample_maxwell(prim::T) where {T<:AA{<:Real,1}}

    u = prim[2]
    if length(prim) == 3
        v = 0.0
        w = 0.0
    elseif length(prim) == 4
        v = prim[3]
        w = 0.0
    elseif length(prim) == 5
        v = prim[3]
        w = prim[4]
    end

    rd1 = rand(3)
    rd2 = rand(3)

    return @. sqrt(-log(rd1) / prim[end]) * sin(2.0 * π * rd2) + [u, v, w]

end

"""
$(SIGNATURES)

Truncated sampling with Distributions.jl
"""
function sample_maxwell(λ::Real, l, u)
    d = truncated(Normal(0.0, sqrt(1.0 / λ)), l, u)
    return rand(d)
end

"""
$(SIGNATURES)
"""
function sample_maxwell(λ::Real, v, l, u)
    d = truncated(Normal(0.0, sqrt(1.0 / λ)), l, u)
    return rand(d) + v
end

"""
$(SIGNATURES)
"""
function sample_maxwell(prim::T, l, u) where {T<:AA{<:Real,1}}
    u = prim[2]
    if length(prim) == 3
        v = 0.0
        w = 0.0
    elseif length(prim) == 4
        v = prim[3]
        w = 0.0
    elseif length(prim) == 5
        v = prim[3]
        w = prim[4]
    end

    d = truncated(Normal(0.0, sqrt(1.0 / prim[end])), l, u)
    return rand(d, 3) + [u, v, w]
end


"""
$(SIGNATURES)

Calculate the time interval until next collision
"""
next_collision_time(τ) = -τ * log(rand())


"""
$(SIGNATURES)

Calculate traveling time to edges
"""
function boundary_time(x, v::Real, xL, xR)
    if v < 0.0
        tb = (xL - x) / v
    elseif v > 0.0
        tb = (xR - x) / v
    else
        tb = 1e8
    end

    return tb
end

"""
$(SIGNATURES)
"""
function boundary_time(x, v::T, xL, xR) where {T<:AA{<:Real,1}}
    if v[1] < 0.0
        tb = (xL - x) / v[1]
    elseif v[1] > 0.0
        tb = (xR - x) / v[1]
    else
        tb = 1e8
    end

    return tb
end
