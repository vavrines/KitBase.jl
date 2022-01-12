# ============================================================
# General Structs
# ============================================================

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
    boundary::E = "fix"
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
