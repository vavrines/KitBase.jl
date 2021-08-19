# ============================================================
# Structs of Initial & Boundary Conditions
# ============================================================

"""
    mutable struct IB{T} <: AbstractCondition
        fw::Function
        bc::T
    end

Initial & boundary conditions with no distribution function

"""
mutable struct IB{T<:Union{Function,AbstractArray}} <: AbstractCondition
    fw::Function
    bc::T
end

function IB(fw, gas::Scalar)
    γ = gas.a

    bc = function(args...)
        w = fw(args...)
        return ifelse(γ == 0, conserve_prim(w), conserve_prim(w, γ))
    end

    return IB{typeof(bc)}(fw, bc)
end

function IB(fw, gas::AbstractGas)
    bc = function(args...)
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

    return IB{typeof(bc)}(fw, bc)
end


"""
    mutable struct IB1F{T} <: AbstractCondition
        fw::Function
        ff::Function
        bc::T
    end

Initial & boundary condition with 1 distribution function

"""
mutable struct IB1F{T} <: AbstractCondition
    fw::Function
    ff::Function
    bc::T
end

function IB1F(fw, vs::AbstractVelocitySpace, gas::AbstractProperty)
    bc = function(args...)
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

    ff = function(args...)
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


"""
    mutable struct IB2F{T} <: AbstractCondition
        fw::Function
        ff::Function
        bc::T
    end

Initial & boundary condition with 2 distribution functions

"""
mutable struct IB2F{T} <: AbstractCondition
    fw::Function
    ff::Function
    bc::T
end

function IB2F(fw, vs::AbstractVelocitySpace, gas::AbstractProperty)
    bc = function(args...)
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

    ff = function(args...)
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


"""
    mutable struct IB3F{T} <: AbstractCondition
        fw::Function
        ff::Function
        fE::Function
        fB::Function
        fL::Function
        bc::T
    end

Initial & boundary condition with 3 distribution functions

"""
mutable struct IB3F{T} <: AbstractCondition
    fw::Function
    ff::Function
    fE::Function
    fB::Function
    fL::Function
    bc::T
end


"""
    mutable struct IB4F{T} <: AbstractCondition
        fw::Function
        ff::Function
        fE::Function
        fB::Function
        fL::Function
        bc::T
    end

Initial & boundary condition with 4 distribution functions

"""
mutable struct IB4F{T} <: AbstractCondition
    fw::Function
    ff::Function
    fE::Function
    fB::Function
    fL::Function
    bc::T
end
