# ============================================================
# General Structs
# ============================================================

"""
$(TYPEDEF)

Config named tuple
"""
Config = @with_kw (
    t0 = 0,
    t1 = 5,
    nt = 16,
    x0 = 0,
    x1 = 1,
    nx = 100,
    y0 = 0,
    y1 = 1,
    ny = 45,
    z0 = 0,
    z1 = 1,
    nz = 45,
    u0 = -5,
    u1 = 5,
    nu = 48,
    v0 = -5,
    v1 = 5,
    nv = 28,
    w0 = -5,
    w1 = 5,
    nw = 28,
    nm = 5,
    Kn = 1,
    K = 0,
    α = 1.0,
    ω = 0.5,
    Pr = 2 / 3,
)


"""
$(SIGNATURES)

Generate config named tuple
"""
function config_ntuple(nt = Config; kwargs...)
    y = nt()
    
    ks = keys(kwargs)
    vs = values(kwargs)
    
    d = Dict()
    d1 = Dict()
    for i in eachindex(ks)
        if haskey(y, ks[i])
            d[ks[i]] = vs[i]
        else
            d1[ks[i]] = vs[i]
        end
    end
        
    merge(nt(; d...), (; d1...))
end


"""
$(TYPEDEF)

Computational setup

# Fields

$(FIELDS)
"""
@kwdef struct Setup{S,I<:Integer,E<:AV,F<:Real,G<:Real} <: AbstractSetup
    matter::S = "gas"
    case::S = "dev"
    space::S = "1d0f0v"
    flux::S = "kfvs"
    collision::S = "bgk"
    nSpecies::I = 1
    interpOrder::I = 2
    limiter::S = "vanleer"
    boundary::E = ["fix", "fix"]
    cfl::F = 0.5
    maxTime::G = 0.1
end

function Setup(
    matter,
    case,
    space,
    flux,
    collision,
    ns,
    order,
    limiter,
    bc::T,
    cfl,
    time,
) where {T<:Union{AbstractString,Symbol}}
    boundary = begin
        if parse(Int, space[1]) == 1
            [bc, bc]
        elseif parse(Int, space[1]) == 2
            [bc, bc, bc, bc]
        end
    end

    return Setup{typeof(matter),typeof(ns),typeof(boundary),typeof(cfl),typeof(time)}(
        matter,
        case,
        space,
        flux,
        collision,
        ns,
        order,
        limiter,
        boundary,
        cfl,
        time,
    )
end
