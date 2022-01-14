# ============================================================
# Thermodynamics
# ============================================================

"""
$(SIGNATURES)

Calculate heat capacity ratio (monatomic gas)
"""
function heat_capacity_ratio(K, D::T) where {T<:Integer}
    γ = begin
        if D == 1
            (K + 3.0) / (K + 1.0)
        elseif D == 2
            (K + 4.0) / (K + 2.0)
        elseif D == 3
            (K + 5.0) / (K + 3.0)
        end
    end

    return γ
end

"""
$(SIGNATURES)

Calculate heat capacity ratio (diatomic gas)
"""
function heat_capacity_ratio(K, Nr, D::T) where {T<:Integer}
    γ = begin
        if D == 1
            (K + 3.0 + Nr) / (K + 1.0 + Nr)
        elseif D == 2
            (K + 4.0 + Nr) / (K + 2.0 + Nr)
        elseif D == 3
            (K + 5.0 + Nr) / (K + 3.0 + Nr)
        end
    end

    return γ
end


"""
$(SIGNATURES)

Calculate speed of sound
"""
sound_speed(λ::Real, γ::Real) = (0.5 * γ / λ)^0.5

"""
$(SIGNATURES)
"""
sound_speed(prim::AV, γ) = sound_speed(prim[end], γ)

"""
$(TYPEDSIGNATURES)

Calculate sound speed in mixture
"""
function sound_speed(prim::AM, γ)
    c = similar(prim, axes(prim, 2))
    for j in eachindex(c)
        c[j] = sound_speed(prim[end, j], γ)
    end

    return maximum(c)
end
