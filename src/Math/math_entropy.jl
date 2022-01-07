"""
$(SIGNATURES)

Maxwell Boltzmann entropy
"""
maxwell_boltzmann(f) = f * log(f) - f


"""
$(SIGNATURES)

Prim of Maxwell Boltzmann entropy
"""
maxwell_boltzmann_prime(x) = log(x)


"""
$(SIGNATURES)

Dual of Maxwell Boltzmann entropy
"""
maxwell_boltzmann_dual(f) = exp(f)


"""
$(SIGNATURES)

Dual prim of Maxwell Boltzmann entropy
"""
maxwell_boltzmann_dual_prime(f) = exp(f)


"""
$(SIGNATURES)

Reconstruct mathematical entropy from Legendre dual
"""
function kinetic_entropy(α::AA, m::AA, weights::AV)
    B = KitBase.maxwell_boltzmann_dual_prime.(α' * m)[:]
    return sum(maxwell_boltzmann.(B) .* weights)
end
