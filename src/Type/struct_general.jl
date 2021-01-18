# ============================================================
# General Structs
# ============================================================

"""
Computational setup

@consts: case, space, flux, collision, nSpecies, interpOrder, limiter, cfl, maxTime

"""
struct Setup{S<:AbstractString,I<:Integer,E<:Real,F<:Real} <: AbstractSetup
    case::S
    space::S
    flux::S
    collision::S
    nSpecies::I
    interpOrder::I
    limiter::S
    boundary::S
    cfl::E
    maxTime::F
end

Setup() = Setup{String,Int,Float64,Float64}("sod", "1d1f1v", "kfvs", "bgk", 1, 1, "vanleer", "fix", 0.5, 2.0)


"""
Particle property

@vars: Kn, Ma, Pr, K, γ, ω, αᵣ, ωᵣ, μᵣ, m, np

"""
mutable struct Gas{A,B,C,D,E,F,G,H,I,J,K<:Integer} <: AbstractProperty

    Kn::A
    Ma::B
    Pr::C
    K::D
    γ::E
    ω::F
    αᵣ::G
    ωᵣ::H
    μᵣ::I
    m::J
    np::K

    function Gas(
        _Kn::Union{Real,AbstractArray}, # unified consideration of
        _Ma::Union{Real,AbstractArray}, # 1. deterministic solution, and
        _Pr::Union{Real,AbstractArray}, # 2. uncertainty quantification
        _K::Union{Real,AbstractArray},
        _γ::Union{Real,AbstractArray},
        _ω::Union{Real,AbstractArray},
        _αᵣ::Union{Real,AbstractArray},
        _ωᵣ::Union{Real,AbstractArray},
        _μᵣ::Union{Real,AbstractArray},
        _m = 1e-3::Union{Real,AbstractArray},
        _np = 1000::Union{Integer,AbstractArray},
    )
        Kn = deepcopy(_Kn)
        Ma = deepcopy(_Ma)
        Pr = deepcopy(_Pr)
        K = deepcopy(_K)
        γ = deepcopy(_γ)
        ω = deepcopy(_ω)
        αᵣ = deepcopy(_αᵣ)
        ωᵣ = deepcopy(_ωᵣ)
        μᵣ = deepcopy(_μᵣ)
        m = deepcopy(_m)
        np = deepcopy(_np)

        new{
            typeof(Kn),
            typeof(Ma),
            typeof(Pr),
            typeof(K),
            typeof(γ),
            typeof(ω),
            typeof(αᵣ),
            typeof(ωᵣ),
            typeof(μᵣ),
            typeof(m),
            typeof(np),
        }(
            Kn,
            Ma,
            Pr,
            K,
            γ,
            ω,
            αᵣ,
            ωᵣ,
            μᵣ,
            m,
            np,
        )
    end

end


"""
Diatomic gas property

@vars: Kn, Ma, Pr, K, γ, ω, αᵣ, ωᵣ, μᵣ, m, np

"""
struct DiatomicGas{TA,TI,TF} <: AbstractProperty

    Kn::TA
    Ma::TA
    Pr::TA
    K::TI
    Kr::TI
    γ::TF
    ω::TF
    αᵣ::TF
    ωᵣ::TF
    μᵣ::TF
    T₀::TF
    Z₀::TF
    σ::TF
    ω₁::TF
    ω₂::TF

    function DiatomicGas(
        Kn::Real,
        Ma::Real,
        Pr::Real,
        K::Real,
        Kr::Real,
        γ::Real,
        ω::Real,
        αᵣ::Real,
        ωᵣ::Real,
        μᵣ::Real,
        T₀::Real,
        Z₀ = 18.1::Real,
        σ = 1 / 1.55::Real,
        ω₁ = 0.2354::Real,
        ω₂ = 0.3049::Real,
    )
        _Ma = typeof(Kn)(Ma)
        _Pr = typeof(Kn)(Pr)
        _Kr = typeof(K)(Kr)
        
        T = typeof(γ)
        _ω = T(ω)
        _αᵣ = T(αᵣ)
        _ωᵣ = T(ωᵣ)
        _μᵣ = T(μᵣ)
        _T₀ = T(T₀)
        _Z₀ = T(Z₀)
        _σ = T(σ)
        _ω₁ = T(ω₁)
        _ω₂ = T(ω₂)

        new{
            typeof(Kn),
            typeof(K),
            T,
        }(
            Kn,
            _Ma,
            _Pr,
            K,
            _Kr,
            γ,
            _ω,
            _αᵣ,
            _ωᵣ,
            _μᵣ,
            _T₀,
            _Z₀,
            _σ,
            _ω₁,
            _ω₂,
        )
    end

end


"""
Multi-component gas property

@consts: Kn, Ma, Pr, K, γ, mi, ni, me, ne

"""
struct Mixture{A,B,C,D,E,F,G,H,I} <: AbstractProperty

    Kn::A
    Ma::B
    Pr::C
    K::D
    γ::E
    mi::F
    ni::G
    me::H
    ne::I

    function Mixture(
        Kn::Array,
        Ma::Union{Real,AbstractArray},
        Pr::Union{Real,AbstractArray},
        K::Union{Real,AbstractArray},
        γ::Union{Real,AbstractArray},
        mi::Real,
        ni::Real,
        me::Real,
        ne::Real,
    )

        # inner constructor method
        new{
            typeof(Kn),
            typeof(Ma),
            typeof(Pr),
            typeof(K),
            typeof(γ),
            typeof(mi),
            typeof(ni),
            typeof(me),
            typeof(ne),
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
        )
    end

end


"""
1D plasma property

@consts: Kn, Ma, Pr, K, γ, mi, ni, me, ne, lD, rL, sol, χ, ν, Ap, An, D

"""
struct Plasma1D{A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P} <: AbstractProperty

    Kn::A
    Ma::B
    Pr::C
    K::D
    γ::E

    mi::F
    ni::G
    me::H
    ne::I
    lD::J
    rL::K

    sol::L
    χ::M
    ν::N
    Ap::O
    An::O
    D::P

    # unified consideration of deterministic and stochastic conditions
    # {mi, ni, me, ne, sol, χ, ν} keep the same as before
    function Plasma1D(
        Kn::Array,
        Ma::Union{Real,AbstractArray},
        Pr::Union{Real,AbstractArray},
        K::Union{Real,AbstractArray},
        γ::Union{Real,AbstractArray},
        mi::Real,
        ni::Real,
        me::Real,
        ne::Real,
        lD::Union{Real,AbstractArray},
        rL::Union{Real,AbstractArray},
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
        new{
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

end


"""
2D plasma property

@consts: Kn, Ma, Pr, K, γ, mi, ni, me, ne, lD, rL, sol, χ, ν, A1p, A1n, A2p, A2n, D1, D2

"""
struct Plasma2D{A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P} <: AbstractProperty

    Kn::A
    Ma::B
    Pr::C
    K::D
    γ::E

    mi::F
    ni::G
    me::H
    ne::I
    lD::J
    rL::K

    sol::L
    χ::M
    ν::N
    A1p::O
    A1n::O
    A2p::O
    A2n::O
    D1::P
    D2::P

    # unified consideration of deterministic and stochastic conditions
    # {mi, ni, me, ne, sol, χ, ν} keep the same as before
    function Plasma2D(
        Kn::Array,
        Ma::Union{Real,AbstractArray},
        Pr::Union{Real,AbstractArray},
        K::Union{Real,AbstractArray},
        γ::Union{Real,AbstractArray},
        mi::Real,
        ni::Real,
        me::Real,
        ne::Real,
        lD::Union{Real,AbstractArray},
        rL::Union{Real,AbstractArray},
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
        new{
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

end
