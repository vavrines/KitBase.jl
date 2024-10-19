"""
$(TYPEDEF)

Fluid property for scalar conservation laws

## Fields

$(FIELDS)
"""
@with_kw mutable struct Scalar{T,F} <: AbstractProperty
    a::T = 1.0
    μᵣ::F = 1e-8
end

"""
$(TYPEDEF)

Radiation property for linear Boltzmann equation

## Fields

$(FIELDS)
"""
@with_kw mutable struct Radiation{T,A,F,I<:Integer} <: AbstractProperty
    Kn::T = 1.0
    σs::A = 1.0
    σa::A = 0.0
    m::F = 1e-3
    np::I = 1000
end

function Radiation(_Kn, _ss::Union{Real,AV}, _sa::Union{Real,AV})
    Kn = deepcopy(_Kn)
    σs = deepcopy(_ss)
    σa = deepcopy(_sa)
    m = 1e-3
    np = 1000

    return Radiation{typeof(Kn),typeof(σs),typeof(m),typeof(np)}(Kn, σs, σa, m, np)
end

"""
$(TYPEDEF)

Gas property

## Fields

$(FIELDS)
"""
@with_kw mutable struct Gas{T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11<:Integer,T12} <: AbstractGas
    Kn::T1 = 1e-2
    Ma::T2 = 0.0
    Pr::T3 = 1.0
    K::T4 = 2.0
    γ::T5 = 5 / 3
    ω::T6 = 0.81
    αᵣ::T7 = 1.0
    ωᵣ::T8 = 0.5
    μᵣ::T9 = ref_vhs_vis(Kn, αᵣ, ωᵣ)
    m::T10 = 1e-3
    np::T11 = 1000
    fsm::T12 = nothing
end

"""
$(TYPEDEF)

Diatomic gas property

## Fields

$(FIELDS)
"""
@with_kw struct DiatomicGas{TA,TI,TF} <: AbstractGas
    Kn::TA = 1e-2
    Ma::TA = 0.0
    Pr::TA = 1.0
    K::TI = 2.0
    Kr::TI = 2.0
    γ::TF = 7 / 5
    ω::TF = 0.74
    αᵣ::TF = 1.0
    ωᵣ::TF = 0.5
    μᵣ::TF = ref_vhs_vis(Kn, αᵣ, ωᵣ)
    T₀::TF = 89.1 / 273
    Z₀::TF = 18.1
    σ::TF = 1 / 1.55
    ω₁::TF = 0.2354
    ω₂::TF = 0.3049
end

"""
$(TYPEDEF)

Multi-component gas property

## Fields

$(FIELDS)
"""
@with_kw struct Mixture{A,B,C,D,E,F,G,H,I} <: AbstractGas
    Kn::A = 1e-2
    Ma::B = 0.0
    Pr::C = 1.0
    K::D = 2.0
    γ::E = 5 / 3
    mi::F = 1.0
    ni::G = 0.5
    me::H = 0.5
    ne::I = 0.5
end

"""
$(TYPEDEF)

Multi-component gas property with internal degrees of freedom

## Fields

$(FIELDS)
"""
@with_kw struct PolyatomicMixture{TA,TI,TF,TV} <: AbstractGas
    Kn::TA = 1e-2
    Ma::TA = 0.0
    Pr::TA = 1.0
    K::TI = 2.0
    Kr::TI = 2.0
    γ::TF = 7 / 5
    T₀::TF = 89.1 / 273
    Z₀::TF = 18.1
    m::TV = [1.0, 0.5]
    n::TV = [0.5, 0.5]
end

"""
$(TYPEDEF)

1D plasma property

## Fields

$(FIELDS)
"""
@with_kw struct Plasma1D{A,B,C,TD,E,F,G,H,I,J,TK,L,M,N,O,P} <: AbstractPlasma
    Kn::A = 1e-2
    Ma::B = 0.0
    Pr::C = 1.0
    K::TD = 2.0
    γ::E = 5 / 3

    mi::F = 1.0
    ni::G = 0.5
    me::H = 0.5
    ne::I = 0.5
    lD::J = 0.01
    rL::TK = 0.01

    sol::L = 100
    χ::M = 1
    ν::N = 1
    Ap::O
    An::O
    D::P
end

# unified consideration of deterministic and stochastic conditions
# {mi, ni, me, ne, sol, χ, ν} keep the same as before
function Plasma1D(
    Kn::Array,
    Ma::Union{Real,AA},
    Pr::Union{Real,AA},
    K::Union{Real,AA},
    γ::Union{Real,AA},
    mi::Real,
    ni::Real,
    me::Real,
    ne::Real,
    lD::Union{Real,AA},
    rL::Union{Real,AA},
    sol::Real,
    χ::Real,
    ν::Real,
)
    # A^+
    Ap = Array{Float64}(undef, 8, 8)
    Ap[1, 1] = (sol * χ) / 2.0
    Ap[7, 1] = χ / 2.0
    Ap[2, 2] = sol / 2.0
    Ap[6, 2] = 0.5
    Ap[3, 3] = sol / 2.0
    Ap[5, 3] = -1.0 / 2.0
    Ap[4, 4] = (sol * ν) / 2.0
    Ap[8, 4] = (sol^2 * ν) / 2.0
    Ap[3, 5] = -sol^2 / 2.0
    Ap[5, 5] = sol / 2.0
    Ap[2, 6] = sol^2 / 2.0
    Ap[6, 6] = sol / 2.0
    Ap[1, 7] = (sol^2 * χ) / 2.0
    Ap[7, 7] = (sol * χ) / 2.0
    Ap[4, 8] = ν / 2.0
    Ap[8, 8] = (sol * ν) / 2.0

    # A^-
    An = Array{Float64}(undef, 8, 8)
    An[1, 1] = -(sol * χ) / 2.0
    An[7, 1] = χ / 2.0
    An[2, 2] = -sol / 2.0
    An[6, 2] = 1.0 / 2.0
    An[3, 3] = -sol / 2.0
    An[5, 3] = -1.0 / 2.0
    An[4, 4] = -(sol * ν) / 2.0
    An[8, 4] = (sol^2 * ν) / 2.0
    An[3, 5] = -sol^2 / 2.0
    An[5, 5] = -sol / 2.0
    An[2, 6] = sol^2 / 2.0
    An[6, 6] = -sol / 2.0
    An[1, 7] = (sol^2 * χ) / 2.0
    An[7, 7] = -(sol * χ) / 2.0
    An[4, 8] = ν / 2.0
    An[8, 8] = -(sol * ν) / 2.0

    # eigenvalues of A
    D = [sol, sol, sol * χ, sol * ν, -sol, -sol, -sol * χ, -sol * ν]

    # inner constructor method
    return Plasma1D{
        typeof(Kn),
        typeof(Ma),
        typeof(Pr),
        typeof(K),
        typeof(γ),
        typeof(mi),
        typeof(ni),
        typeof(me),
        typeof(ne),
        typeof(lD),
        typeof(rL),
        typeof(sol),
        typeof(χ),
        typeof(ν),
        typeof(Ap),
        typeof(D),
    }(
        Kn,
        Ma,
        Pr,
        K,
        γ,
        mi,
        ni,
        me,
        ne,
        lD,
        rL,
        sol,
        χ,
        ν,
        Ap,
        An,
        D,
    )
end

"""
$(TYPEDEF)

2D plasma property

## Fields

$(FIELDS)
"""
@with_kw struct Plasma2D{A,B,C,D,E,F,G,H,I,J,TK,L,M,N,O,P} <: AbstractPlasma
    Kn::A = 1e-2
    Ma::B = 0.0
    Pr::C = 1.0
    K::D = 2
    γ::E = 5 / 3

    mi::F = 1.0
    ni::G = 0.5
    me::H = 0.5
    ne::I = 0.5
    lD::J = 0.01
    rL::TK = 0.01

    sol::L = 100
    χ::M = 1
    ν::N = 1
    A1p::O
    A1n::O
    A2p::O
    A2n::O
    D1::P
    D2::P
end

# unified consideration of deterministic and stochastic conditions
# {mi, ni, me, ne, sol, χ, ν} keep the same as before
function Plasma2D(
    Kn::Array,
    Ma::Union{Real,AA},
    Pr::Union{Real,AA},
    K::Union{Real,AA},
    γ::Union{Real,AA},
    mi::Real,
    ni::Real,
    me::Real,
    ne::Real,
    lD::Union{Real,AA},
    rL::Union{Real,AA},
    sol::Real,
    χ::Real,
    ν::Real,
)
    # A₁+
    A1p = zeros(8, 8)
    A1p[1, 1] = (sol * χ) / 2.0
    A1p[7, 1] = χ / 2.0
    A1p[2, 2] = sol / 2.0
    A1p[6, 2] = 0.5
    A1p[3, 3] = sol / 2.0
    A1p[5, 3] = -1.0 / 2.0
    A1p[4, 4] = (sol * ν) / 2.0
    A1p[8, 4] = (sol^2 * ν) / 2.0
    A1p[3, 5] = -sol^2 / 2.0
    A1p[5, 5] = sol / 2.0
    A1p[2, 6] = sol^2 / 2.0
    A1p[6, 6] = sol / 2.0
    A1p[1, 7] = (sol^2 * χ) / 2.0
    A1p[7, 7] = (sol * χ) / 2.0
    A1p[4, 8] = ν / 2.0
    A1p[8, 8] = (sol * ν) / 2.0

    # A₁-
    A1n = zeros(8, 8)
    A1n[1, 1] = -(sol * χ) / 2.0
    A1n[7, 1] = χ / 2.0
    A1n[2, 2] = -sol / 2.0
    A1n[6, 2] = 1.0 / 2.0
    A1n[3, 3] = -sol / 2.0
    A1n[5, 3] = -1.0 / 2.0
    A1n[4, 4] = -(sol * ν) / 2.0
    A1n[8, 4] = (sol^2 * ν) / 2.0
    A1n[3, 5] = -sol^2 / 2.0
    A1n[5, 5] = -sol / 2.0
    A1n[2, 6] = sol^2 / 2.0
    A1n[6, 6] = -sol / 2.0
    A1n[1, 7] = (sol^2 * χ) / 2.0
    A1n[7, 7] = -(sol * χ) / 2.0
    A1n[4, 8] = ν / 2.0
    A1n[8, 8] = -(sol * ν) / 2.0

    D1 = [sol, sol, sol * χ, sol * ν, -sol, -sol, -sol * χ, -sol * ν]

    # A₂+
    A2p = zeros(8, 8)
    A2p[1, 1] = sol / 2.0
    A2p[6, 1] = -1.0 / 2.0
    A2p[2, 2] = (sol * χ) / 2.0
    A2p[7, 2] = χ / 2.0
    A2p[3, 3] = sol / 2.0
    A2p[4, 3] = 1.0 / 2.0
    A2p[3, 4] = (sol^2) / 2.0
    A2p[4, 4] = sol / 2.0
    A2p[5, 5] = (sol * ν) / 2.0
    A2p[8, 5] = (sol^2 * ν) / 2.0
    A2p[1, 6] = -sol^2 / 2.0
    A2p[6, 6] = sol / 2.0
    A2p[2, 7] = (sol^2 * χ) / 2.0
    A2p[7, 7] = (sol * χ) / 2.0
    A2p[5, 8] = ν / 2.0
    A2p[8, 8] = (sol * ν) / 2.0

    # A₂-
    A2n = zeros(8, 8)
    A2n[1, 1] = -sol / 2.0
    A2n[6, 1] = -1.0 / 2.0
    A2n[2, 2] = -(sol * χ) / 2.0
    A2n[7, 2] = χ / 2.0
    A2n[3, 3] = -sol / 2.0
    A2n[4, 3] = 1.0 / 2.0
    A2n[3, 4] = (sol^2) / 2.0
    A2n[4, 4] = -sol / 2.0
    A2n[5, 5] = -(sol * ν) / 2.0
    A2n[8, 5] = (sol^2 * ν) / 2.0
    A2n[1, 6] = -sol^2 / 2.0
    A2n[6, 6] = -sol / 2.0
    A2n[2, 7] = (sol^2 * χ) / 2.0
    A2n[7, 7] = -(sol * χ) / 2.0
    A2n[5, 8] = ν / 2.0
    A2n[8, 8] = -(sol * ν) / 2.0

    D2 = [sol, sol, sol * χ, sol * ν, -sol, -sol, -sol * χ, -sol * ν]

    # inner constructor method
    return Plasma2D{
        typeof(Kn),
        typeof(Ma),
        typeof(Pr),
        typeof(K),
        typeof(γ),
        typeof(mi),
        typeof(ni),
        typeof(me),
        typeof(ne),
        typeof(lD),
        typeof(rL),
        typeof(sol),
        typeof(χ),
        typeof(ν),
        typeof(A1p),
        typeof(D1),
    }(
        Kn,
        Ma,
        Pr,
        K,
        γ,
        mi,
        ni,
        me,
        ne,
        lD,
        rL,
        sol,
        χ,
        ν,
        A1p,
        A1n,
        A2p,
        A2n,
        D1,
        D2,
    )
end

"""
$(SIGNATURES)

Generate property of matter
"""
function set_property(dict::Union{AbstractDict,NamedTuple})
    matter = dict[:matter]
    fn = eval(Symbol(matter * "_property"))
    return fn(; dict...)
end

function scalar_property(; a, mu, kwargs...)
    return Scalar(a, mu)
end

function radiation_property(; knudsen, sigmaS, sigmaA, kwargs...)
    return Radiation(knudsen, sigmaS, sigmaA)
end

function gas_property(;
    space,
    inK,
    nSpecies,
    knudsen,
    mach,
    prandtl,
    omega=0.5,
    alphaRef=1.0,
    omegaRef=0.5,
    mi=nothing,
    ni=nothing,
    me=nothing,
    ne=nothing,
    kwargs...,
)
    Dx = parse(Int, space[1]) # this is type instable
    γD = map(parse(Int, space[3]), parse(Int, space[5])) do x, y # (x)f(y)v
        if x == 0
            return Dx
        elseif x == 1 # 1f
            if y >= 3 # 3v
                return 3
            else # 1v
                return Dx
            end
        elseif x == 2 # 2f
            return Dx
        elseif x >= 3 # 3f / 4f
            return 3
        else
            return nothing
        end
    end
    γ = heat_capacity_ratio(inK, γD)

    if nSpecies == 1
        μᵣ = ref_vhs_vis(knudsen, alphaRef, omegaRef)

        gas = Gas(;
            Kn=knudsen,
            Ma=mach,
            Pr=prandtl,
            K=inK,
            γ=γ,
            ω=omega,
            αᵣ=alphaRef,
            ωᵣ=omegaRef,
            μᵣ=μᵣ,
        )
    elseif nSpecies == 2
        kne = knudsen * (me / mi)
        gas = Mixture([knudsen, kne], mach, prandtl, inK, γ, mi, ni, me, ne)
    end

    return gas
end

function gas_property(
    vs;
    space,
    inK,
    knudsen,
    mach,
    prandtl,
    omega,
    alphaRef,
    omegaRef,
    collision,
    nm=5,
    kwargs...,
)
    @assert collision == "fsm"

    Dx = parse(Int, space[1])
    γD = map(parse(Int, space[3]), parse(Int, space[5])) do x, y # (x)f(y)v
        if x == 0
            return Dx
        elseif x == 1 # 1f
            if y >= 3 # 3v
                return 3
            else # 1v
                return Dx
            end
        elseif x == 2 # 2f
            return Dx
        elseif x >= 3 # 3f / 4f
            return 3
        else
            return nothing
        end
    end
    γ = heat_capacity_ratio(inK, γD)

    μᵣ = ref_vhs_vis(knudsen, alphaRef, omegaRef)
    fsm = fsm_kernel(vs, μᵣ, nm, alphaRef)

    return Gas(;
        Kn=knudsen,
        Ma=mach,
        Pr=prandtl,
        K=inK,
        γ=γ,
        ω=omega,
        αᵣ=alphaRef,
        ωᵣ=omegaRef,
        μᵣ=μᵣ,
        fsm=fsm,
    )
end

function plasma_property(;
    space,
    inK,
    knudsen,
    mach,
    prandtl,
    mi,
    ni,
    me,
    ne,
    lD,
    rL,
    sol,
    echi,
    bnu,
    kwargs...,
)
    Dx = parse(Int, space[1])
    γD = map(parse(Int, space[3]), parse(Int, space[5])) do x, y # (x)f(y)v
        if x == 0
            return Dx
        elseif x == 1 # 1f
            if y >= 3 # 3v
                return 3
            else # 1v
                return Dx
            end
        elseif x == 2 # 2f
            return Dx
        elseif x >= 3 # 3f / 4f
            return 3
        else
            return nothing
        end
    end
    γ = heat_capacity_ratio(inK, γD)

    kne = knudsen * (me / mi)
    if Dx == 1
        gas = Plasma1D(
            [knudsen, kne],
            mach,
            prandtl,
            inK,
            γ,
            mi,
            ni,
            me,
            ne,
            lD,
            rL,
            sol,
            echi,
            bnu,
        )
    elseif Dx == 2
        gas = Plasma2D(
            [knudsen, kne],
            mach,
            prandtl,
            inK,
            γ,
            mi,
            ni,
            me,
            ne,
            lD,
            rL,
            sol,
            echi,
            bnu,
        )
    end

    return gas
end
