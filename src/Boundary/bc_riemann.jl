"""
Inflow/outflow boundary condition using Riemann invariants

## Arguments

  - `ctr`: ghost control volume
  - `ctr0`: interior control volume
  - `u`,`v`: velocity collocation points
  - `prim`: primitive variables in the reference state
  - `γ`: heat capacity ratio
  - `K`: internal degrees of freedom
  - `n`: pointing ctr0 -> ctr
"""
function bc_riemann!(
    ctr::T,
    ctr0::T,
    u,
    v,
    prim,
    γ,
    K,
    n,
) where {T<:Union{ControlVolume2F,ControlVolume2D2F,ControlVolumeUS2F}}
    # interior cell
    bci = local_frame(ctr0.prim, n[1], n[2])
    ci = sound_speed(bci, γ)

    # outerior cell
    bco = local_frame(prim, n[1], n[2])
    co = sound_speed(bco, γ)

    # Riemann invariants
    Rplus = bci[2] + 2.0 * ci / (γ - 1.0)
    Rminus = bco[2] - 2.0 * co / (γ - 1.0)

    # primitive variables
    Ub = 0.5 * (Rplus + Rminus)
    cb = 0.25 * (γ - 1.0) * (Rplus - Rminus)

    Vb, sb = begin
        if Ub > 0 # outflow
            bci[3], fluid_entropy(bci[1], ci, γ)
        elseif Ub <= 0 # inflow
            bco[3], fluid_entropy(bco[1], co, γ)
        end
    end

    ρb = (cb^2 / γ / sb)^(1.0 / γ - 1.0)
    pb = ρb * cb^2 / γ

    ctr.prim .= global_frame([ρb, Ub, Vb, 0.5 * ρb / pb], n[1], n[2])
    ctr.w .= prim_conserve(ctr.prim, γ)
    ctr.sw .= 0.0
    ctr.h .= maxwellian(u, v, ctr.prim)
    ctr.b .= energy_maxwellian(ctr.h, ctr.prim, K)
    ctr.sh .= 0.0
    ctr.sb .= 0.0

    return nothing
end
