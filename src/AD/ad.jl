# ============================================================
# Automatic Differentiation Methods
# ============================================================

function ∂maxwellian(u::Real, ρ, U, λ)
    Mu = u -> maxwellian(u, ρ, U, λ)
    return ForwardDiff.derivative(Mu, u)
end

∂maxwellian(u::Real, prim::AV) = ∂maxwellian(u, prim[1], prim[2], prim[end])

function ∂maxwellian(u::AV, ρ, U, λ)
    ∂maxwellian.(u::AV, ρ, U, λ)
end

∂maxwellian(u::AV, prim::AV) = ∂maxwellian(u, prim[1], prim[2], prim[end])
