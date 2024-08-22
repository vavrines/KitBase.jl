# ============================================================
# Structs of Initial & Boundary Conditions
# ============================================================

"""
$(TYPEDEF)

Initial & boundary conditions with no distribution function

## Fields

$(FIELDS)
"""
mutable struct IB{TF,T,NT} <: AbstractCondition
    fw::TF
    bc::T
    p::NT
end

(m::IB)(args...) = begin
    m.fw(args..., m.p), m.bc(args..., m.p)
end


"""
$(TYPEDEF)

Initial & boundary conditions with 1 distribution function

## Fields

$(FIELDS)
"""
mutable struct IB1F{TF1,TF2,T,NT} <: AbstractCondition
    fw::TF1
    ff::TF2
    bc::T
    p::NT
end

(m::IB1F)(args...) = begin
    m.fw(args..., m.p), m.ff(args..., m.p), m.bc(args..., m.p)
end


"""
$(TYPEDEF)

Initial & boundary conditions with 2 distribution functions

## Fields

$(FIELDS)
"""
mutable struct IB2F{TF1,TF2,T,NT} <: AbstractCondition
    fw::TF1
    ff::TF2
    bc::T
    p::NT
end

(m::IB2F)(args...) = begin
    m.fw(args..., m.p), m.ff(args..., m.p)[1], m.ff(args..., m.p)[2], m.bc(args..., m.p)
end


"""
$(TYPEDEF)

Initial & boundary conditions with 3 distribution functions

## Fields

$(FIELDS)
"""
mutable struct IB3F{TF1,TF2,TF3,TF4,TF5,T,NT} <: AbstractCondition
    fw::TF1
    ff::TF2
    fE::TF3
    fB::TF4
    fL::TF5
    bc::T
    p::NT
end


"""
$(TYPEDEF)

Initial & boundary conditions with 4 distribution functions

## Fields

$(FIELDS)
"""
mutable struct IB4F{TF1,TF2,TF3,TF4,TF5,T,NT} <: AbstractCondition
    fw::TF1
    ff::TF2
    fE::TF3
    fB::TF4
    fL::TF5
    bc::T
    p::NT
end


"""
$(SIGNATURES)

Generate initial & boundary conditions
"""
function set_ib(dict::Union{AbstractDict,NamedTuple}, set, ps, vs, gas)
    if haskey(dict, :uLid)
        ib = set_ib(set, ps, vs, gas, dict[:uLid], dict[:vLid], dict[:tLid])
    else
        ib = set_ib(set, ps, vs, gas)
    end

    return ib
end

"""
$(SIGNATURES)
"""
function set_ib(set, pSpace, vSpace, gas, Um = 0.15, Vm = 0.0, Tm = 1.0)
    ib = begin
        if parse(Int, set.space[3]) == 0
            fw, bc, p = config_ib(set, pSpace, vSpace, gas)
            IB(fw, bc, p)
        elseif parse(Int, set.space[3]) in [3, 4] && gas isa AbstractPlasma
            fw, ff, fE, fB, fL, bc, p = config_ib(set, pSpace, vSpace, gas)
            iname = "IB" * set.space[3] * "F"
            eval(Symbol(iname))(fw, ff, fE, fB, fL, bc, p)
        elseif set.case == "cavity"
            fw, ff, bc, p = config_ib(set, pSpace, vSpace, gas, Um, Vm, Tm)
            iname = "IB" * set.space[3] * "F"
            eval(Symbol(iname))(fw, ff, bc, p)
        else
            fw, ff, bc, p = config_ib(set, pSpace, vSpace, gas)
            iname = "IB" * set.space[3] * "F"
            eval(Symbol(iname))(fw, ff, bc, p)
        end
    end

    return ib
end
