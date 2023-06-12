# ============================================================
# Mathematical Methods
# ============================================================

export linspace, heaviside, fortsign, mat_split, dirac_delta
export central_diff, central_diff!, central_diff2, central_diff2!
export upwind_diff, upwind_diff!, unstruct_diff
export ylm, rlylm, eval_spherharmonic
export basis_size, eval_sphermonomial
export polylog

include("math_general.jl")
include("math_difference.jl")
include("sphere_harmonics.jl")
include("sphere_monomials.jl")
include("polylog.jl")
