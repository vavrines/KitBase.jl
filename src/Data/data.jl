# ============================================================
# Parameters & Data Transformer
# ============================================================

export static_array, slope_array

const BZ = 1.38065e-23
const M_Ar = 6.63e-26
const D_Ar = 3.66e-10

include("conversion.jl")

symbolize(x::AbstractString) = Symbol(x)
symbolize(x::AbstractArray{T}) where {T<:AbstractString} = [symbolize(y) for y in x]
