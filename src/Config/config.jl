# ============================================================
# Initial & Boundary Conditions of Specific Problems
# ============================================================

export ib_rh,
       ib_sod,
       ib_briowu,
       ib_cavity

include("cfg_rh.jl")
include("cfg_sod.jl")
include("cfg_briowu.jl")
include("cfg_cavity.jl")
