# ============================================================
# Balance Functions
# ============================================================

function bc_balance!(
    ctr::T,
    ctr0::T,
    ctr1::T,
) where {T<:Union{ControlVolume,ControlVolume1D,ControlVolume2D,ControlVolumeUS}}
    ctr.w = 2.0 * ctr0.w - ctr1.w
    @. ctr.prim = 2.0 * ctr0.prim - ctr1.prim

    return nothing
end

function bc_balance!(
    ctr::T,
    ctr0::T,
    ctr1::T,
) where {T<:Union{ControlVolume1F,ControlVolume1D1F,ControlVolume2D1F,ControlVolumeUS1F}}
    @. ctr.w = 2.0 * ctr0.w - ctr1.w
    @. ctr.prim = 2.0 * ctr0.prim - ctr1.prim
    @. ctr.f = 2.0 * ctr0.f - ctr1.f

    return nothing
end

function bc_balance!(
    ctr::T,
    ctr0::T,
    ctr1::T,
) where {T<:Union{ControlVolume2F,ControlVolume1D2F,ControlVolume2D2F,ControlVolumeUS2F}}
    @. ctr.w = 2.0 * ctr0.w - ctr1.w
    @. ctr.prim = 2.0 * ctr0.prim - ctr1.prim
    @. ctr.h = 2.0 * ctr0.h - ctr1.h
    @. ctr.b = 2.0 * ctr0.b - ctr1.b

    return nothing
end

function bc_balance!(
    ctr::T,
    ctr0::T,
    ctr1::T,
) where {T<:Union{ControlVolume1D3F,ControlVolume2D3F}}
    @. ctr.w = 2.0 * ctr0.w - ctr1.w
    @. ctr.prim = 2.0 * ctr0.prim - ctr1.prim
    @. ctr.h0 = 2.0 * ctr0.h0 - ctr1.h0
    @. ctr.h1 = 2.0 * ctr0.h1 - ctr1.h1
    @. ctr.h2 = 2.0 * ctr0.h2 - ctr1.h2
    @. ctr.E = 2.0 * ctr0.E - ctr1.E
    @. ctr.B = 2.0 * ctr0.B - ctr1.B
    ctr.ϕ = 2.0 * ctr0.ϕ - ctr1.ϕ
    ctr.ψ = 2.0 * ctr0.ψ - ctr1.ψ
    @. ctr.lorenz = 2.0 * ctr0.lorenz - ctr1.lorenz

    return nothing
end

function bc_balance!(ctr::T, ctr0::T, ctr1::T) where {T<:ControlVolume1D4F}
    @. ctr.w = 2.0 * ctr0.w - ctr1.w
    @. ctr.prim = 2.0 * ctr0.prim - ctr1.prim
    @. ctr.h0 = 2.0 * ctr0.h0 - ctr1.h0
    @. ctr.h1 = 2.0 * ctr0.h1 - ctr1.h1
    @. ctr.h2 = 2.0 * ctr0.h2 - ctr1.h2
    @. ctr.E = 2.0 * ctr0.E - ctr1.E
    @. ctr.B = 2.0 * ctr0.B - ctr1.B
    ctr.ϕ = 2.0 * ctr0.ϕ - ctr1.ϕ
    ctr.ψ = 2.0 * ctr0.ψ - ctr1.ψ
    @. ctr.lorenz = 2.0 * ctr0.lorenz - ctr1.lorenz

    return nothing
end
