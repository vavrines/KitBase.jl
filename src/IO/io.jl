# ============================================================
# I / O Methods
# ============================================================

export read_cfg, read_dict
export write_sol, write_jld, write_vtk, write_tec

include("io_read.jl")
include("io_write.jl")
include("io_plot.jl")
