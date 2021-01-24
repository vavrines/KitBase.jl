# ============================================================
# Thermodynamics
# ============================================================

"""
    heat_capacity_ratio(K, D::T) where {T<:Integer}

Calculate heat capacity ratio

"""
function heat_capacity_ratio(K, D::T) where {T<:Integer}

    if D == 1
        γ = (K + 3.0) / (K + 1.0)
    elseif D == 2
        γ = (K + 4.0) / (K + 2.0)
    elseif D == 3
        γ = (K + 5.0) / (K + 3.0)
    end

    return γ

end

#--- diatomic ---#
function heat_capacity_ratio(K, Nr, D::T) where {T<:Integer}

    if D == 1
        γ = (K + 3.0 + Nr) / (K + 1.0 + Nr)
    elseif D == 2
        γ = (K + 4.0 + Nr) / (K + 2.0 + Nr)
    elseif D == 3
        γ = (K + 5.0 + Nr) / (K + 3.0 + Nr)
    end

    return γ

end


"""
    sound_speed(λ::Real, γ::Real)
    sound_speed(prim::T, γ) where {T<:AbstractArray{<:Real,1}}
    sound_speed(prim::T, γ) where {T<:AbstractArray{<:Real,2}}

Calculate speed of sound

"""
sound_speed(λ::Real, γ::Real) = (0.5 * γ / λ)^0.5

sound_speed(prim::T, γ) where {T<:AbstractArray{<:Real,1}} = sound_speed(prim[end], γ)

# mixture
function sound_speed(prim::T, γ) where {T<:AbstractArray{<:Real,2}}
    c = similar(prim, axes(prim, 2))
    for j in eachindex(c)
        c[j] = sound_speed(prim[end, j], γ)
    end

    return maximum(c)
end
