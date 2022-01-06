# ============================================================
# Geometrical Methods
# ============================================================

export global_frame, local_frame, unit_normal, point_distance
export PSpace1D, PSpace2D, CSpace2D
export uniform_mesh, ndgrid, meshgrid, find_idx
export UnstructPSpace

include("geo_general.jl")
include("geo_struct.jl")
include("geo_unstruct.jl")
