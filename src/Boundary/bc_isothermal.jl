# ============================================================
# Isothermal Wall Functions
# ============================================================

"""
$(SIGNATURES)

Isothermal boundary condition
"""
function bc_isothermal!(
    ctr::T,
    ctr1::T,
    γ::RN,
    λw = 1.0::RN,
) where {T<:Union{ControlVolume1D,ControlVolume2D,ControlVolumeUS}}

    ctr.prim[2] = -ctr1.prim[2]
    if length(ctr.prim) == 4
        ctr.prim[3] = -ctr1.prim[3]
    elseif length(ctr.prim) == 5
        ctr.prim[4] = -ctr1.prim[4]
    end

    ctr.prim[end] = 2.0 * λw - ctr1.prim[end]
    fac = (ctr1.prim[end] - λw) / λw
    ctr.prim[1] = (1.0 - fac) / (1.0 + fac) * ctr1.prim[1]

    ctr.w .= prim_conserve(ctr.prim, γ)

    return nothing
end
