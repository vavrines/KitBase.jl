"""
    hs_boltz_kn(mu_ref, alpha)

Calculate effective Knudsen number for fast spectral method with hard sphere (HS) model
    
"""
hs_boltz_kn(mu_ref, alpha) =
    64 * sqrt(2.0)^alpha / 5.0 * gamma((alpha + 3) / 2) * gamma(2.0) * sqrt(pi) * mu_ref


"""
    kernel_mode(
        M::TI,
        umax::TR,
        vmax::TR,
        du::TR,
        dv::TR,
        unum::TI,
        vnum::TI;
        quad_num = 64::TI,
    ) where {TI<:Integer,TR<:Real}

    kernel_mode(
        M::I,
        umax::R,
        vmax::R,
        wmax::R,
        du::R,
        dv::R,
        dw::R,
        unum::I,
        vnum::I,
        wnum::I,
        alpha::R;
        quad_num = 64::I,
    ) where {I<:Integer,R<:Real}

Calculate collision kernel for fast spectral method

"""
function kernel_mode(
    M::I,
    umax::R,
    vmax::R,
    wmax::R,
    du::R,
    dv::R,
    dw::R,
    unum::I,
    vnum::I,
    wnum::I,
    alpha::R;
    quad_num = 64::I,
) where {I<:Integer,R<:Real}

    supp = sqrt(2.0) * 2.0 * max(umax, vmax, wmax) / (3.0 + sqrt(2.0))

    fre_vx = range(-π / du, (unum ÷ 2 - 1) * 2.0 * π / unum / du, length = unum)
    fre_vy = range(-π / dv, (vnum ÷ 2 - 1) * 2.0 * π / vnum / dv, length = vnum)
    fre_vz = range(-π / dw, (wnum ÷ 2 - 1) * 2.0 * π / wnum / dw, length = wnum)

    # abscissa, gweight = gausslegendre(quad_num)
    # @. abscissa = (0. * (1. - abscissa) + supp * (1. + abscissa)) / 2
    # @. gweight *= (supp - 0.) / 2

    abscissa, gweight = lgwt(quad_num, 0, supp)

    phi = zeros(unum, vnum, wnum, M * (M - 1))
    psi = zeros(unum, vnum, wnum, M * (M - 1))
    phipsi = zeros(unum, vnum, wnum)
    for loop = 1:M-1
        theta = π / M * loop
        for loop2 = 1:M
            theta2 = π / M * loop2
            idx = (loop - 1) * M + loop2
            for k = 1:wnum, j = 1:vnum, i = 1:unum
                s =
                    fre_vx[i] * sin(theta) * cos(theta2) +
                    fre_vy[j] * sin(theta) * sin(theta2) +
                    fre_vz[k] * cos(theta)
                # phi
                int_temp = 0.0
                for id = 1:quad_num
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

function kernel_mode(
    M::I,
    umax::R,
    vmax::R,
    wmax::R,
    unum::I,
    vnum::I,
    wnum::I,
    alpha::R;
    quad_num = 64::I,
) where {I<:Integer,R<:Real}

    supp = sqrt(2.0) * 2.0 * max(umax, vmax, wmax) / (3.0 + sqrt(2.0))

    fre_vx = [i * π / umax for i = -unum÷2:unum÷2-1]
    fre_vy = [i * π / vmax for i = -vnum÷2:vnum÷2-1]
    fre_vz = [i * π / wmax for i = -wnum÷2:wnum÷2-1]

    abscissa, gweight = KitBase.lgwt(quad_num, 0, supp)

    phi = zeros(unum, vnum, wnum, M * (M - 1))
    psi = zeros(unum, vnum, wnum, M * (M - 1))
    phipsi = zeros(unum, vnum, wnum)
    for loop = 1:M-1
        theta = π / M * loop
        for loop2 = 1:M
            theta2 = π / M * loop2
            idx = (loop - 1) * M + loop2
            for k = 1:wnum, j = 1:vnum, i = 1:unum
                s =
                    fre_vx[i] * sin(theta) * cos(theta2) +
                    fre_vy[j] * sin(theta) * sin(theta2) +
                    fre_vz[k] * cos(theta)
                # phi
                int_temp = 0.0
                for id = 1:quad_num
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

function kernel_mode(
    M::TI,
    umax::TR,
    vmax::TR,
    du::TR,
    dv::TR,
    unum::TI,
    vnum::TI;
    quad_num = 64::TI,
) where {TI<:Integer,TR<:Real}

    supp = sqrt(2.0) * 2.0 * max(umax, vmax) / (3.0 + sqrt(2.0))

    fre_vx = range(-π / du, (unum ÷ 2 - 1) * 2.0 * π / unum / du, length = unum)
    fre_vy = range(-π / dv, (vnum ÷ 2 - 1) * 2.0 * π / vnum / dv, length = vnum)
    Fre_vx, Fre_vy = ndgrid(fre_vx, fre_vy)

    αp = zeros(unum, vnum, M)
    αq = zeros(unum, vnum, M)
    for lp = 1:M
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
    boltzmann_fft(
        f::X,
        Kn,
        M::I,
        ϕ::Y,
        ψ::Y,
        phipsi::Z,
    ) where {
        X<:AA{<:Real,3},
        Y<:AA{<:Real,4},
        Z<:AA{<:Real,3},
        I<:Integer,
    }

Calculate collision operator with FFT-based fast spectral method

"""
function boltzmann_fft(
    f::X,
    Kn,
    M::I,
    ϕ::Y,
    ψ::Y,
    phipsi::Z,
) where {
    X<:AA{<:Real,3},
    Y<:AA{<:Real,4},
    Z<:AA{<:Real,3},
    I<:Integer,
}

    f_spec = f .+ 0im
    bfft!(f_spec)
    f_spec ./= size(f, 1) * size(f, 2) * size(f, 3)
    f_spec .= fftshift(f_spec)

    #--- gain term ---#
    f_temp = zeros(axes(f_spec)) .+ 0im
    for i = 1:M*(M-1)
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
    boltzmann_fft!(
        Q::T1,
        f::T2,
        Kn::TR,
        M::TI,
        ϕ::TY,
        ψ::TY,
        phipsi::TZ,
    ) where {
        T1<:AA{<:Real,3},
        T2<:AA{<:Real,3},
        TR<:Real,
        TI<:Integer,
        TY<:AA{<:Real,4},
        TZ<:AA{<:Real,3},
    }

Calculate collision operator with FFT-based fast spectral method

"""
function boltzmann_fft!(
    Q::T1,
    f::T2,
    Kn::TR,
    M::TI,
    ϕ::TY,
    ψ::TY,
    phipsi::TZ,
) where {
    T1<:AA{<:Real,3},
    T2<:AA{<:Real,3},
    TR<:Real,
    TI<:Integer,
    TY<:AA{<:Real,4},
    TZ<:AA{<:Real,3},
}

    f_spec = f .+ 0im
    bfft!(f_spec)
    f_spec ./= size(f, 1) * size(f, 2) * size(f, 3)
    f_spec .= fftshift(f_spec)

    # --- gain term ---#
    f_temp = zeros(axes(f_spec)) .+ 0im
    for i = 1:M*(M-1)
        fg1 = f_spec .* ϕ[:, :, :, i]
        fg2 = f_spec .* ψ[:, :, :, i]
        fg11 = fft(fg1)
        fg22 = fft(fg2)
        f_temp .+= fg11 .* fg22
    end

    # --- loss term ---#
    fl1 = f_spec .* phipsi
    fl2 = f_spec
    fl11 = fft(fl1)
    fl22 = fft(fl2)
    f_temp .-= fl11 .* fl22

    @. Q = 4.0 * π^2 / Kn / M^2 * real(f_temp)

end
