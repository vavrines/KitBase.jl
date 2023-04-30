# ============================================================
# General Structs
# ============================================================

"""
$(TYPEDEF)

Config named tuple
"""
Config = @with_kw (
    matter = "gas",
    case = "dev",
    space = "1d1f1v",
    flux = "kfvs",
    collision = "bgk",
    nSpecies = 1,
    interpOrder = 1,
    limiter = "vanleer",
    boundary = ["fix", "fix"],
    hasForce = false,
    cfl = 0.5,
    t0 = 0,
    t1 = 1,
    maxTime = t1,
    nt = 16,
    x0 = 0,
    x1 = 1,
    nx = 100,
    nxg = 0,
    y0 = 0,
    y1 = 1,
    ny = 45,
    nyg = 0,
    z0 = 0,
    z1 = 1,
    nz = 45,
    nzg = 0,
    u0 = -5,
    u1 = 5,
    umin = u0,
    umax = u1,
    nu = 48,
    nug = 0,
    v0 = -5,
    v1 = 5,
    vmin = v0,
    vmax = v1,
    nv = 28,
    nvg = 0,
    w0 = -5,
    w1 = 5,
    wmin = w0,
    wmax = w1,
    nw = 28,
    nwg = 0,
    vMeshType = "rectangle",
    nm = 5,
    Kn = 1,
    knudsen = Kn,
    Ma = 0,
    mach = Ma,
    K = 0,
    inK = K,
    α = 1.0,
    alphaRef = α,
    ω = 0.5,
    omega = ω,
    omegaRef = ω,
    Pr = 2 / 3,
    prandtl = Pr,
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

## Fields

$(FIELDS)
"""
@with_kw struct Setup{S,I<:Integer,E<:AV,F<:Real,G<:Real} <: AbstractSetup
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
    hasForce::Bool = false
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
    hasForce,
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
        hasForce,
    )
end

"""
$(SIGNATURES)

Generate Setup
"""
function set_setup(;
    matter,
    case,
    space,
    flux,
    collision,
    nSpecies,
    interpOrder,
    limiter,
    boundary,
    cfl,
    maxTime,
    hasForce,
    kwargs...,
)
    set = Setup(
        matter,
        case,
        space,
        flux,
        collision,
        nSpecies,
        interpOrder,
        limiter,
        boundary,
        cfl,
        maxTime,
        hasForce,
    )

    return set
end

"""
$(SIGNATURES)
"""
set_setup(dict::Union{AbstractDict,NamedTuple}) = set_setup(; dict...)
