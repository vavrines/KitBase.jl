"""
    maxwell_boltzmann(f)

Maxwell Boltzmann entropy
"""
maxwell_boltzmann(f) = f * log(f) - f


"""
    maxwell_boltzmann_prime(x)

Prim of Maxwell Boltzmann entropy
"""
maxwell_boltzmann_prime(x) = log(x)


"""
    maxwell_boltzmann_dual(f)

Dual of Maxwell Boltzmann entropy
"""
maxwell_boltzmann_dual(f) = exp(f)


"""
    maxwell_boltzmann_dual_prime(f)

Dual prim of Maxwell Boltzmann entropy
"""
maxwell_boltzmann_dual_prime(f) = exp(f)


"""
    kinetic_entropy(α, m, weights)

Reconstruct mathematical entropy from Legendre dual
"""
function kinetic_entropy(α::AbstractArray, m::AbstractArray, weights::AbstractVector)
    B = KitBase.maxwell_boltzmann_dual_prime.(α' * m)[:]
    return sum(maxwell_boltzmann.(B) .* weights)
end
