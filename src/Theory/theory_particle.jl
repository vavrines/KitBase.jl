"""
    sample_velocity(prim::AbstractArray{<:Real,1})

Sample particle velocity from Maxwellian

"""
function sample_velocity(v, λ)
    rd1 = rand()
    rd2 = rand()

    return sqrt(-log(rd1) / λ) * sin(2.0 * π * rd2) + v
end

function sample_velocity(prim::T) where {T<:AbstractArray{<:Real,1}}
    rd1 = rand(3)
    rd2 = rand(3)

    return @. sqrt(-log(rd1) / prim[end]) * sin(2.0 * π * rd2) + [prim[2], 0., 0.]
end

function sample_velocity(prim::T, l, u) where {T<:AbstractArray{<:Real,1}}
    rd1 = rand(3)
    rd2 = rand(3)

    d = truncated(Normal(0.0, sqrt(1.0/prim[3])), l, u)
    return rand(d, 3) + [prim[2], 0., 0.]
end



"""
    next_collision_time(τ)

Calculate the time interval until next collision

"""
next_collision_time(τ) = -τ * log(rand())


function boundary_time(x, v::T, x0, dx0) where {T<:AbstractArray{<:Real,1}}
    if v[1] < 0
        tb = (x0 - 0.5 * dx0 - x) / v[1]
    elseif v[1] > 0
        tb = (x0 + 0.5 * dx0 - x) / v[1]
    else
        tb = 1e8
    end

    return tb
end