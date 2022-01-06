# ============================================================
# Mathematical Methods
# ============================================================

export linspace, heaviside, fortsign, mat_split
export finite_difference
export central_diff, central_diff!, central_diff2, central_diff2!
export upwind_diff, upwind_diff!, unstruct_diff
export maxwell_boltzmann, maxwell_boltzmann_prime
export maxwell_boltzmann_dual, maxwell_boltzmann_dual_prime
export kinetic_entropy
export optimize_closure, realizable_reconstruct
export ylm, rlylm, eval_spherharmonic
export basis_size, eval_sphermonomial

include("math_general.jl")
include("math_difference.jl")
include("math_entropy.jl")
include("math_closure.jl")
include("sphere_harmonics.jl")
include("sphere_monomials.jl")
