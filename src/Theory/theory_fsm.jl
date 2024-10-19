"""
$(SIGNATURES)

RHS-ODE of Boltzmann equation
"""
function boltzmann_ode!(df, f::AA{T,3}, p, t) where {T}
    Kn, M, phi, psi, phipsi = p
    return df .= boltzmann_fft(f, Kn, M, phi, psi, phipsi)
end

"""
$(SIGNATURES)

Calculate effective Knudsen number for fast spectral method with hard sphere (HS) model
"""
hs_boltz_kn(mu_ref, alpha) =
    64 * sqrt(2.0)^alpha / 5.0 * gamma((alpha + 3) / 2) * gamma(2.0) * sqrt(pi) * mu_ref

"""
$(TYPEDSIGNATURES)

Calculate collision kernel for fast spectral method...
"""
function kernel_mode(
    M::I,
    umax::R,
    vmax::R,
    wmax::R,
    du::R1,
    dv::R1,
    dw::R1,
    unum::I,
    vnum::I,
    wnum::I,
    alpha::R2;
    quad_num=64::I,
) where {I<:Integer,R,R1,R2}
    supp = sqrt(2.0) * 2.0 * max(umax, vmax, wmax) / (3.0 + sqrt(2.0))

    fre_vx = range(-π / du, (unum ÷ 2 - 1) * 2.0 * π / unum / du; length=unum)
    fre_vy = range(-π / dv, (vnum ÷ 2 - 1) * 2.0 * π / vnum / dv; length=vnum)
    fre_vz = range(-π / dw, (wnum ÷ 2 - 1) * 2.0 * π / wnum / dw; length=wnum)

    # abscissa, gweight = gausslegendre(quad_num)
    # @. abscissa = (0. * (1. - abscissa) + supp * (1. + abscissa)) / 2
    # @. gweight *= (supp - 0.) / 2

    abscissa, gweight = lgwt(quad_num, 0, supp)

    phi = zeros(unum, vnum, wnum, M * (M - 1))
    psi = zeros(unum, vnum, wnum, M * (M - 1))
    phipsi = zeros(unum, vnum, wnum)
    for loop in 1:M-1
        theta = π / M * loop
        for loop2 in 1:M
            theta2 = π / M * loop2
            idx = (loop - 1) * M + loop2
            for k in 1:wnum, j in 1:vnum, i in 1:unum
                s =
                    fre_vx[i] * sin(theta) * cos(theta2) +
                    fre_vy[j] * sin(theta) * sin(theta2) +
                    fre_vz[k] * cos(theta)
                # phi
                int_temp = 0.0
                for id in 1:quad_num
                    int_temp +=
                        2.0 * gweight[id] * cos(s * abscissa[id]) * (abscissa[id]^alpha)
                end
                phi[i, j, k, idx] = int_temp * sin(theta)
                # psi
                s = fre_vx[i]^2 + fre_vy[j]^2 + fre_vz[k]^2 - s^2
                if s <= 0.0
                    psi[i, j, k, idx] = π * supp^2
                else
                    s = sqrt(s)
                    bel = supp * s
                    bessel = besselj(1, bel)
                    psi[i, j, k, idx] = 2.0 * π * supp * bessel / s
                end
                # phipsi
                phipsi[i, j, k] += phi[i, j, k, idx] * psi[i, j, k, idx]
            end
        end
    end

    return phi, psi, phipsi
end

"""
$(TYPEDSIGNATURES)
"""
function kernel_mode(
    M::I,
    umax::R,
    vmax::R,
    wmax::R,
    unum::I,
    vnum::I,
    wnum::I,
    alpha::R1;
    quad_num=64::I,
) where {I<:Integer,R,R1}
    supp = sqrt(2.0) * 2.0 * max(umax, vmax, wmax) / (3.0 + sqrt(2.0))

    fre_vx = [i * π / umax for i in -unum÷2:unum÷2-1]
    fre_vy = [i * π / vmax for i in -vnum÷2:vnum÷2-1]
    fre_vz = [i * π / wmax for i in -wnum÷2:wnum÷2-1]

    abscissa, gweight = KitBase.lgwt(quad_num, 0, supp)

    phi = zeros(unum, vnum, wnum, M * (M - 1))
    psi = zeros(unum, vnum, wnum, M * (M - 1))
    phipsi = zeros(unum, vnum, wnum)
    for loop in 1:M-1
        theta = π / M * loop
        for loop2 in 1:M
            theta2 = π / M * loop2
            idx = (loop - 1) * M + loop2
            for k in 1:wnum, j in 1:vnum, i in 1:unum
                s =
                    fre_vx[i] * sin(theta) * cos(theta2) +
                    fre_vy[j] * sin(theta) * sin(theta2) +
                    fre_vz[k] * cos(theta)
                # phi
                int_temp = 0.0
                for id in 1:quad_num
                    int_temp +=
                        2.0 * gweight[id] * cos(s * abscissa[id]) * (abscissa[id]^alpha)
                end
                phi[i, j, k, idx] = int_temp * sin(theta)
                # psi
                s = fre_vx[i]^2 + fre_vy[j]^2 + fre_vz[k]^2 - s^2
                if s <= 0.0
                    psi[i, j, k, idx] = π * supp^2
                else
                    s = sqrt(s)
                    bel = supp * s
                    bessel = besselj(1, bel)
                    psi[i, j, k, idx] = 2.0 * π * supp * bessel / s
                end
                # phipsi
                phipsi[i, j, k] += phi[i, j, k, idx] * psi[i, j, k, idx]
            end
        end
    end

    return phi, psi, phipsi
end

"""
$(TYPEDSIGNATURES)
"""
function kernel_mode(
    M::TI,
    umax::TR,
    vmax::TR,
    du::TR1,
    dv::TR1,
    unum::TI,
    vnum::TI;
    quad_num=64::TI,
) where {TI<:Integer,TR,TR1}
    supp = sqrt(2.0) * 2.0 * max(umax, vmax) / (3.0 + sqrt(2.0))

    fre_vx = range(-π / du, (unum ÷ 2 - 1) * 2.0 * π / unum / du; length=unum)
    fre_vy = range(-π / dv, (vnum ÷ 2 - 1) * 2.0 * π / vnum / dv; length=vnum)
    Fre_vx, Fre_vy = ndgrid(fre_vx, fre_vy)

    αp = zeros(unum, vnum, M)
    αq = zeros(unum, vnum, M)
    for lp in 1:M
        θ = lp * π / M

        s = Fre_vx * cos(θ + π / 2) + Fre_vy * sin(θ + π / 2)
        αq[:, :, lp] .= 2.0 * supp .* sin.(supp .* s) ./ (supp .* s) * (π / M)

        s = Fre_vx * cos(θ) + Fre_vy * sin(θ)
        αp[:, :, lp] .= 2.0 * supp .* sin.(supp .* s) ./ (supp .* s)
    end
    αp ./= π

    pq = αp .* αq
    α0 = zeros(unum, vnum)
    for j in axes(α0, 2), i in axes(α0, 1)
        α0[i, j] = sum(pq[i, j, :])
    end

    return αp, αq, α0
end

"""
$(TYPEDSIGNATURES)

Create pre-computed kernel for fast spectral method
"""
function fsm_kernel(vs::AbstractVelocitySpace, μ, nm=5, α=1.0)
    kn_bz = hs_boltz_kn(μ, α)

    phi, psi, phipsi = kernel_mode(
        nm,
        vs.u1,
        vs.v1,
        vs.w1,
        vs.du[1, 1, 1],
        vs.dv[1, 1, 1],
        vs.dw[1, 1, 1],
        vs.nu,
        vs.nv,
        vs.nw,
        α,
    )

    return (Kn=kn_bz, nm=nm, ϕ=phi, ψ=psi, χ=phipsi)
end

"""
$(TYPEDSIGNATURES)

Calculate collision operator with FFT-based fast spectral method...
"""
function boltzmann_fft(
    f0::AA{<:Real,3},
    Kn,
    M::Integer,
    ϕ::Y,
    ψ::Y,
    phipsi::AA{<:Real,3},
) where {Y<:AA{<:Real,4}}
    f = begin
        if f0 isa Array
            f0
        else
            Array(f0)
        end
    end

    f_spec = f .+ 0im
    bfft!(f_spec)
    f_spec ./= size(f, 1) * size(f, 2) * size(f, 3)
    f_spec .= fftshift(f_spec)

    #--- gain term ---#
    f_temp = zeros(axes(f_spec)) .+ 0im
    for i in 1:M*(M-1)
        fg1 = f_spec .* ϕ[:, :, :, i]
        fg2 = f_spec .* ψ[:, :, :, i]
        fg11 = fft(fg1)
        fg22 = fft(fg2)
        f_temp .+= fg11 .* fg22
    end

    #--- loss term ---#
    fl1 = f_spec .* phipsi
    fl2 = f_spec
    fl11 = fft(fl1)
    fl22 = fft(fl2)
    f_temp .-= fl11 .* fl22

    Q = @. 4.0 * π^2 / Kn / M^2 * real(f_temp)

    return Q
end

"""
$(SIGNATURES)
"""
boltzmann_fft(f, p) = boltzmann_fft(f, p.Kn, p.nm, p.ϕ, p.ψ, p.χ)

"""
$(TYPEDSIGNATURES)

Calculate collision operator with FFT-based fast spectral method...
"""
function boltzmann_fft!(
    Q::AA{<:Real,3},
    f0::AA{<:Real,3},
    Kn::Real,
    M::Integer,
    ϕ::T,
    ψ::T,
    phipsi::AA{<:Real,3},
) where {T<:AA{<:Real,4}}
    f = begin
        if f0 isa Array
            f0
        else
            Array(f0)
        end
    end

    f_spec = f .+ 0im
    bfft!(f_spec)
    f_spec ./= size(f, 1) * size(f, 2) * size(f, 3)
    f_spec .= fftshift(f_spec)

    #--- gain term ---#
    f_temp = zeros(axes(f_spec)) .+ 0im
    for i in 1:M*(M-1)
        fg1 = f_spec .* ϕ[:, :, :, i]
        fg2 = f_spec .* ψ[:, :, :, i]
        fg11 = fft(fg1)
        fg22 = fft(fg2)
        f_temp .+= fg11 .* fg22
    end

    #--- loss term ---#
    fl1 = f_spec .* phipsi
    fl2 = f_spec
    fl11 = fft(fl1)
    fl22 = fft(fl2)
    f_temp .-= fl11 .* fl22

    @. Q = 4.0 * π^2 / Kn / M^2 * real(f_temp)
end

"""
$(SIGNATURES)
"""
boltzmann_fft!(Q, f, p) = boltzmann_fft!(Q, f, p.Kn, p.nm, p.ϕ, p.ψ, p.χ)
