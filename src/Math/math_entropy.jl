"""
    maxwell_boltzmann(f)

Maxwell Boltzmann entropy
"""
maxwell_boltzmann(f::T) where {T<:Real} = f * log(f) - f 


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
