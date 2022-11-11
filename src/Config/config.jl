# ============================================================
# Initial & Boundary Conditions of Specific Problems
# ============================================================

export config_ib

include("cfg_advection.jl")
include("cfg_rh.jl")
include("cfg_sod.jl")
include("cfg_briowu.jl")
include("cfg_cavity.jl")

"""
$(SIGNATURES)

Config initial and boundary conditions
"""
function config_ib(args...; case = args[1].case)
    func = begin
        if case in ("shock", :shock)
            eval(Symbol("ib_" * "rh"))
        elseif case in ("brio-wu", Symbol("brio-wu"))
            eval(Symbol("ib_" * "briowu"))
        else
            eval(Symbol("ib_" * string(case)))
        end
    end

    return func(args...)
end
