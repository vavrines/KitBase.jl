# ============================================================
# Structs of Initial & Boundary Conditions
# ============================================================

"""
$(TYPEDEF)

Initial & boundary conditions with no distribution function

# Fields

$(FIELDS)
"""
mutable struct IB{T} <: AbstractCondition
    fw::Function
    bc::T
    p::NamedTuple
end

IB(fw, bc = nothing) = IB{typeof(bc)}(fw, bc, NamedTuple())

(m::IB)(args...) = begin
    m.fw(args..., m.p), m.bc(args..., m.p)
end


"""
$(TYPEDEF)

Initial & boundary conditions with 1 distribution function

# Fields

$(FIELDS)
"""
mutable struct IB1F{T} <: AbstractCondition
    fw::Function
    ff::Function
    bc::T
    p::NamedTuple
end

(m::IB1F)(args...) = begin
    m.fw(args..., m.p), m.ff(args..., m.p), m.bc(args..., m.p)
end

#=
function IB1F(fw, vs::AbstractVelocitySpace, gas::AbstractProperty)
    bc = function (args...)
        w = fw(args...)
        prim = begin
            if ndims(w) == 1
                conserve_prim(w, gas.γ)
            else
                mixture_conserve_prim(w, gas.γ)
            end
        end

        return prim
    end

    ff = function (args...)
        prim = bc(args...)

        M = begin
            if !isdefined(vs, :v)
                if ndims(prim) == 1
                    maxwellian(vs.u, prim)
                else
                    mixture_maxwellian(vs.u, prim)
                end
            elseif !isdefined(vs, :w)
                if ndims(prim) == 1
                    maxwellian(vs.u, vs.v, prim)
                else
                    mixture_maxwellian(vs.u, vs.v, prim)
                end
            else
                if ndims(prim) == 1
                    maxwellian(vs.u, vs.v, vs.w, prim)
                else
                    mixture_maxwellian(vs.u, vs.v, vs.w, prim)
                end
            end
        end

        return M
    end

    return IB1F{typeof(bc)}(fw, ff, bc)
end
=#

"""
$(TYPEDEF)

Initial & boundary conditions with 2 distribution functions

# Fields

$(FIELDS)
"""
mutable struct IB2F{T} <: AbstractCondition
    fw::Function
    ff::Function
    bc::T
    p::NamedTuple
end

(m::IB2F)(args...) = begin
    m.fw(args..., m.p), m.ff(args..., m.p)[1], m.ff(args..., m.p)[2], m.bc(args..., m.p)
end

#=
function IB2F(fw, vs::AbstractVelocitySpace, gas::AbstractProperty)
    bc = function (args...)
        w = fw(args...)
        prim = begin
            if ndims(w) == 1
                conserve_prim(w, gas.γ)
            else
                mixture_conserve_prim(w, gas.γ)
            end
        end

        return prim
    end

    ff = function (args...)
        prim = bc(args...)

        if !isdefined(vs, :v)
            if ndims(prim) == 1
                H = maxwellian(vs.u, prim)
                B = H .* gas.K / 2 / prim[end]
            else
                H = mixture_maxwellian(vs.u, prim)
                B = zero(H)
                for j in axes(B, 2)
                    B[:, j] = H[:, j] * gas.K / (2.0 * prim[end, j])
                end
            end
        elseif !isdefined(vs, :w)
            if ndims(prim) == 1
                H = maxwellian(vs.u, vs.v, prim)
                B = H .* gas.K / 2 / prim[end]
            else
                H = mixture_maxwellian(vs.u, vs.v, prim)
                B = zero(H)
                for j in axes(B, 3)
                    B[:, :, j] = H[:, :, j] * gas.K / (2.0 * prim[end, j])
                end
            end
        else
            if ndims(prim) == 1
                H = maxwellian(vs.u, vs.v, vs.w, prim)
                B = H .* gas.K / 2 / prim[end]
            else
                H = mixture_maxwellian(vs.u, vs.v, vs.w, prim)
                B = zero(H)
                for j in axes(B, 4)
                    B[:, :, :, j] = H[:, :, :, j] * gas.K / (2.0 * prim[end, j])
                end
            end
        end

        return H, B
    end

    return IB2F{typeof(bc)}(fw, ff, bc)
end
=#

"""
$(TYPEDEF)

Initial & boundary conditions with 3 distribution functions

# Fields

$(FIELDS)
"""
mutable struct IB3F{T} <: AbstractCondition
    fw::Function
    ff::Function
    fE::Function
    fB::Function
    fL::Function
    bc::T
    p::NamedTuple
end


"""
$(TYPEDEF)

Initial & boundary conditions with 4 distribution functions

# Fields

$(FIELDS)
"""
mutable struct IB4F{T} <: AbstractCondition
    fw::Function
    ff::Function
    fE::Function
    fB::Function
    fL::Function
    bc::T
    p::NamedTuple
end
