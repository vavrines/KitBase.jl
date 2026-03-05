# ============================================================
# Geometrical Methods
# ============================================================

export PSpace1D, PSpace2D, CSpace2D, PSpace3D, UnstructPSpace
export global_frame, local_frame, uniform_mesh
export SharpIB

include("geo_general.jl")
include("geo_struct.jl")
include("geo_unstruct.jl")
include("geo_ib.jl")

"""
$(SIGNATURES)

Generate PhysicalSpace
"""
function set_geometry(;
    space,
    x0=nothing,
    x1=nothing,
    nx=nothing,
    nxg=nothing,
    y0=nothing,
    y1=nothing,
    ny=nothing,
    nyg=nothing,
    z0=nothing,
    z1=nothing,
    nz=nothing,
    nzg=nothing,
    mesh=nothing,
    kwargs...,
)
    try
        return UnstructPSpace(mesh)
    catch
        Dx = parse(Int, space[1])
        if Dx == 1
            return PSpace1D(x0, x1, nx, nxg)
        elseif Dx == 2
            return PSpace2D(x0, x1, nx, y0, y1, ny, nxg, nyg)
        elseif Dx == 3
            return PSpace3D(x0, x1, nx, y0, y1, ny, z0, z1, nz, nxg, nyg, nzg)
        else
            throw("No preset available for 3D simulation, please set it up manually.")
        end
    end
end

set_geometry(dict::Union{AbstractDict,NamedTuple}) = set_geometry(; dict...)
