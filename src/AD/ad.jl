# ============================================================
# Automatic Differentiation Methods
# ============================================================

"""
$(SIGNATURES)

Calculate derivatives of Maxwellian
"""
function ∂maxwellian(u::Real, ρ, U, λ)
    Mu = u -> maxwellian(u, ρ, U, λ)
    return ForwardDiff.derivative(Mu, u)
end

"""
$(SIGNATURES)
"""
∂maxwellian(u::Real, prim::AV) = ∂maxwellian(u, prim[1], prim[2], prim[end])

"""
$(SIGNATURES)
"""
function ∂maxwellian(u::AV, ρ, U, λ)
    ∂maxwellian.(u::AV, ρ, U, λ)
end

"""
$(SIGNATURES)
"""
∂maxwellian(u::AV, prim::AV) = ∂maxwellian(u, prim[1], prim[2], prim[end])
