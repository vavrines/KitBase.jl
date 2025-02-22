# ============================================================
# Hydrodynamic Fluxes
# Reference: http://www.cfdbooks.com/
# ============================================================

"""
$(SIGNATURES)

Upwind flux
"""
function flux_upwind(uL, uR, Ω::T, n::T, dt=1.0) where {T<:AV}
    ip = dot(Ω, n)

    if ip > 0
        return dt * ip * uL
    else
        return dt * ip * uR
    end
end

"""
$(SIGNATURES)

Lax-Friedrichs flux

_P. D. Lax, Weak Solutions of Nonlinear Hyperbolic Equations and Their Numerical Computation,
Commun. Pure and Applied Mathematics, 7, 159-193, 1954._
"""
function flux_lax!(fw::AV, wL::T, wR::T, γ, dt, dx) where {T<:AV}
    fw .= 0.5 * dt .* (euler_flux(wL, γ)[1] + euler_flux(wR, γ)[1] - dx / dt .* (wR - wL))
    return nothing
end

"""
$(SIGNATURES)

HLL flux
"""
function flux_hll!(fw::AV, wL::T, wR::T, γ, dt, len=1.0) where {T<:AV}
    primL = conserve_prim(wL, γ)
    primR = conserve_prim(wR, γ)

    aL = sound_speed(primL, γ)
    aR = sound_speed(primR, γ)

    λmin = primL[2] - aL
    λmax = primR[2] + aR

    if λmin >= 0.0
        fw .= euler_flux(wL, γ)[1]
    elseif λmax <= 0.0
        fw .= euler_flux(wR, γ)[1]
    else
        factor = 1.0 / (λmax - λmin)

        flux1 = euler_flux(wL, γ)[1]
        flux2 = euler_flux(wR, γ)[1]

        @. fw = factor * (λmax * flux1 - λmin * flux2 + λmax * λmin * (wR - wL))
    end

    @. fw *= dt * len

    return nothing
end

flux_hll!(KS::AbstractSolverSet, face, ctrL, ctrR, args...) =
    flux_hll!(face, ctrL, ctrR, KS.gas, args...)

function flux_hll!(
    face::Interface,
    ctrL::T,
    ctrR::T,
    gas::Gas,
    p,
    dt=1.0,
) where {T<:ControlVolume}
    dxL, dxR = p[1:2]

    if size(ctrL.w, 1) == 3
        flux_hll!(face.fw, ctrL.w .+ dxL .* ctrL.sw, ctrR.w .- dxR .* ctrR.sw, gas.γ, dt)
    elseif size(ctrL.w, 1) == 4
        len, n, dirc = p[3:5]
        swL = @view ctrL.sw[:, dirc]
        swR = @view ctrR.sw[:, dirc]
        flux_hll!(
            face.fw,
            local_frame(ctrL.w .+ dxL .* swL, n),
            local_frame(ctrR.w .- dxR .* swR, n),
            gas.γ,
            dt,
            len,
        )
        face.fw .= global_frame(face.fw, n)
    end

    return nothing
end

"""
$(SIGNATURES)

HLLC flux
"""
function flux_hllc!(fw, wL::T, wR::T, γ, dt, len=1.0) where {T}
    primL = conserve_prim(wL, γ)
    primR = conserve_prim(wR, γ)

    aL = sound_speed(primL, γ)
    aR = sound_speed(primR, γ)

    tmp1 = 0.5 * (aL + aR) # avg sos
    tmp2 = 0.5 * (primL[1] + primR[1]) # avg density

    pL = 0.5 * primL[1] / primL[end]
    pR = 0.5 * primR[1] / primR[end]

    pvars = 0.5 * (pL + pR) - 0.5 * (primR[2] - primL[2]) * tmp1 * tmp2
    pstar = max(0.0, pvars)

    qk = begin
        if pstar < pL
            1.0
        else
            tmp1 = (γ + 1.0) / (2.0 * γ)
            tmp2 = (pstar / pL - 1.0)
            sqrt(1.0 + tmp1 * tmp2)
        end
    end
    sL = primL[2] - aL * qk

    qk = begin
        if pstar <= pR
            1.0
        else
            tmp1 = (γ + 1.0) / (2.0 * γ)
            tmp2 = (pstar / pR - 1.0)
            sqrt(1.0 + tmp1 * tmp2)
        end
    end
    sR = primR[2] + aR * qk

    tmp1 =
        pR - pL + primL[1] * primL[2] * (sL - primL[2]) -
        primR[1] * primR[2] * (sR - primR[2])
    tmp2 = primL[1] * (sL - primL[2]) - primR[1] * (sR - primR[2])
    star = tmp1 / tmp2

    if sL >= 0
        fw .= euler_flux(wL, γ)[1]
    elseif sR <= 0
        fw .= euler_flux(wR, γ)[1]
    elseif star >= 0 && sL <= 0
        fw .= euler_flux(wL, γ)[1]
        qstar = hllc_var(primL, pL, sL, star, γ)
        @. fw += sL * (qstar - wL)
    elseif star <= 0 && sR >= 0
        fw .= euler_flux(wR, γ)[1]
        qstar = hllc_var(primR, pR, sR, star, γ)
        @. fw += sR * (qstar - wR)
    end

    @. fw *= dt * len

    return nothing
end

flux_hllc!(KS::AbstractSolverSet, face, ctrL, ctrR, args...) =
    flux_hllc!(face, ctrL, ctrR, KS.gas, args...)

function flux_hllc!(
    face::Interface,
    ctrL::T,
    ctrR::T,
    gas::Gas,
    p,
    dt=1.0,
) where {T<:ControlVolume}
    dxL, dxR = p[1:2]

    if size(ctrL.w, 1) == 3
        flux_hllc!(face.fw, ctrL.w .+ dxL .* ctrL.sw, ctrR.w .- dxR .* ctrR.sw, gas.γ, dt)
    elseif size(ctrL.w, 1) == 4
        len, n, dirc = p[3:5]
        swL = @view ctrL.sw[:, dirc]
        swR = @view ctrR.sw[:, dirc]
        flux_hllc!(
            face.fw,
            local_frame(ctrL.w .+ dxL .* swL, n),
            local_frame(ctrR.w .- dxR .* swR, n),
            gas.γ,
            dt,
            len,
        )
        face.fw .= global_frame(face.fw, n)
    end

    return nothing
end

function hllc_var(d1, u1, p1, s1, star1, γ)
    tmp1 = d1 * (s1 - u1) / (s1 - star1)
    tmp2 = 0.5 * (u1^2) + p1 / ((γ - 1.0) * d1)
    tmp3 = star1 + p1 / (d1 * (s1 - u1))

    ustar = [tmp1, tmp1 * star1, tmp1 * (tmp2 + (star1 - u1) * tmp3)]

    return ustar
end

function hllc_var(d1, u1, v1, p1, s1, star1, γ)
    tmp1 = d1 * (s1 - u1) / (s1 - star1)
    tmp2 = 0.5 * (u1^2 + v1^2) + p1 / ((γ - 1.0) * d1)
    tmp3 = star1 + p1 / (d1 * (s1 - u1))

    ustar = [tmp1, tmp1 * star1, tmp1 * v1, tmp1 * (tmp2 + (star1 - u1) * tmp3)]

    return ustar
end

function hllc_var(prim, p, s, star, γ)
    if size(prim, 1) == 3
        return hllc_var(prim[1], prim[2], p, s, star, γ)
    elseif size(prim, 1) == 4
        return hllc_var(prim[1], prim[2], prim[3], p, s, star, γ)
    end
end

"""
$(SIGNATURES)

Roe's flux with entropy fix    

_P. L. Roe, Approximate Riemann Solvers, Parameter Vectors and Difference Schemes, Journal of Computational Physics, 43, pp. 357-372._
(_cf. http://cfdbooks.com/cfdcodes.html_)
"""
flux_roe!(KS::AbstractSolverSet, face, ctrL, ctrR, args...) =
    flux_roe!(face, ctrL, ctrR, KS.gas, args...)

"""
$(SIGNATURES)
"""
function flux_roe!(
    face::Interface,
    ctrL::T,
    ctrR::T,
    gas::Gas,
    p,
    dt=1.0,
) where {T<:ControlVolume}
    dxL, dxR = p[1:2]

    if size(ctrL.w, 1) == 3
        flux_roe!(face.fw, ctrL.w .+ dxL .* ctrL.sw, ctrR.w .- dxR .* ctrR.sw, gas.γ, dt)
    else
    end

    return nothing
end

"""
$(SIGNATURES)

## Arguments
* `n`: unit face normal (L -> R)
"""
function flux_roe!(fw::AV, wL::T, wR::T, γ, dt, δs=1.0, n=[1.0, 0.0]) where {T<:AV}
    primL = conserve_prim(wL, γ)
    primR = conserve_prim(wR, γ)

    if length(fw) == 3 # 1D

        # left state
        rhoL = primL[1]
        vL = primL[2]
        pL = 0.5 * primL[1] / primL[end]
        aL = sound_speed(primL, γ)
        HL = (wL[end] + pL) / rhoL

        # right state
        rhoR = primR[1]
        vR = primR[2]
        pR = 0.5 * primR[1] / primR[end]
        aR = sound_speed(primR, γ)
        HR = (wR[end] + pR) / rhoR

        # compute Roe averages
        RT = sqrt(rhoR / rhoL)
        rho = RT * rhoL
        v = (vL + RT * vR) / (1.0 + RT)
        H = (HL + RT * HR) / (1.0 + RT)
        a = sqrt((γ - 1.0) * (H - 0.5 * v^2))

        # differences in primitive variables
        drho = rhoR - rhoL
        du = vR - vL
        dP = pR - pL

        # wave strength (characteristic variables)
        dV = [
            0.5 * (dP - rho * a * du) / a^2,
            -(dP / a^2 - drho),
            0.5 * (dP + rho * a * du) / a^2,
        ]

        # absolute values of the wave speeds (eigenvalues)
        ws = [abs(v - a), abs(v), abs(v + a)]

        # modified wave speeds for nonlinear fields (to remove expansion shocks).
        Da = max(0.0, 4.0 * ((vR - aR) - (vL - aL)))
        if ws[1] < 0.5 * Da
            ws[1] = ws[1]^2 / Da + 0.25 * Da
        end
        Da = max(0.0, 4.0 * ((vR + aR) - (vL + aL)))
        if ws[3] < 0.5 * Da
            ws[3] = ws[3]^2 / Da + 0.25 * Da
        end

        # right eigenvectors
        R = zeros(3, 3)
        R[1, 1] = 1.0
        R[2, 1] = v - a
        R[3, 1] = H - v * a

        R[1, 2] = 1.0
        R[2, 2] = v
        R[3, 2] = 0.5 * v * v

        R[1, 3] = 1.0
        R[2, 3] = v + a
        R[3, 3] = H + v * a

        # compute average flux
        fw .= 0.5 .* (euler_flux(wL, γ)[1] + euler_flux(wR, γ)[1])

        # add matrix dissipation term to complete Roe flux
        for j in 1:3, k in 1:3
            fw[j] -= 0.5 * ws[k] * dV[k] * R[j, k]
        end

    elseif length(fw) == 4 # 2D

        # normal vector
        nx = n[1]
        ny = n[2]

        # tangent vector
        mx = -ny
        my = nx

        #--- primitive and other variables ---#
        # left state
        rhoL = primL[1]
        uL = primL[2]
        vL = primL[3]
        unL = uL * nx + vL * ny
        umL = uL * mx + vL * my
        pL = 0.5 * primL[1] / primL[4]
        aL = sound_speed(primL[4], γ)
        HL = aL^2 / (γ - 1.0) + 0.5 * (uL^2 + vL^2)

        # right state
        rhoR = primR[1]
        uR = primR[2]
        vR = primR[3]
        unR = uR * nx + vR * ny
        umR = uR * mx + vR * my
        pR = 0.5 * primR[1] / primR[4]
        aR = sound_speed(primR[4], γ)
        HR = aR^2 / (γ - 1.0) + 0.5 * (uR^2 + vR^2)

        # Roe averages
        RT = sqrt(rhoR / rhoL)
        rho = RT * rhoL
        u = (uL + RT * uR) / (1.0 + RT)
        v = (vL + RT * vR) / (1.0 + RT)
        H = (HL + RT * HR) / (1.0 + RT)
        a = sqrt((γ - 1.0) * (H - 0.5 * (u^2 + v^2)))
        un = u * nx + v * ny
        um = u * mx + v * my

        # wave strengths
        drho = rhoR - rhoL
        dp = pR - pL
        dun = unR - unL
        dum = umR - umL

        LdU = [
            (dp - rho * a * dun) / (2.0 * a^2),
            rho * dum,
            drho - dp / (a^2),
            (dp + rho * a * dun) / (2.0 * a^2),
        ]

        # wave speed
        ws = abs.([un - a, un, un, un + a])

        # Harten's entropy fix JCP(1983), 49, pp357-393
        # only for the nonlinear fields.
        if ws[1] < 0.2
            ws[1] = 0.5 * (ws[1]^2 / 0.2 + 0.2)
        end
        if ws[4] < 0.2
            ws[4] = 0.5 * (ws[4]^2 / 0.2 + 0.2)
        end

        # right eigenvectors
        Rv = zeros(4, 4)

        Rv[1, 1] = 1.0
        Rv[2, 1] = u - a * nx
        Rv[3, 1] = v - a * ny
        Rv[4, 1] = H - un * a

        Rv[1, 2] = 0.0
        Rv[2, 2] = mx
        Rv[3, 2] = my
        Rv[4, 2] = um

        Rv[1, 3] = 1.0
        Rv[2, 3] = u
        Rv[3, 3] = v
        Rv[4, 3] = 0.5 * (u * u + v * v)

        Rv[1, 4] = 1.0
        Rv[2, 4] = u + a * nx
        Rv[3, 4] = v + a * ny
        Rv[4, 4] = H + un * a

        # dissipation term
        diss = zeros(4)
        for i in 1:4, j in 1:4
            diss[i] += ws[j] * LdU[j] * Rv[i, j]
        end

        # compute fluxes
        fL = [
            rhoL * unL,
            rhoL * unL * uL + pL * nx,
            rhoL * unL * vL + pL * ny,
            rhoL * unL * HL,
        ]
        fR = [
            rhoR * unR,
            rhoR * unR * uR + pR * nx,
            rhoR * unR * vR + pR * ny,
            rhoR * unR * HR,
        ]

        @. fw = 0.5 * dt * (fL + fR - diss) * δs
    end

    return nothing
end


"""
$(SIGNATURES)

Godunov flux based on the exact solution of Riemann problem
"""
flux_godunov!(KS::AbstractSolverSet, face, ctrL, ctrR, args...) =
    flux_godunov!(face, ctrL, ctrR, KS.gas, args...)

"""
$(SIGNATURES)
"""
function flux_godunov!(
    face::Interface,
    ctrL::T,
    ctrR::T,
    gas::Gas,
    p,
    dt=1.0,
) where {T<:ControlVolume}
    if size(ctrL.w, 1) == 3
        flux_godunov!(face.fw, ctrL.w, ctrR.w, gas.γ, dt)
    else
    end

    return nothing
end

"""
$(SIGNATURES)
"""
function flux_godunov!(fw, wL, wR, γ, dt=1.0, len=1.0; tol=1e-6)
    U_star = ustar(wL, wR, γ; tol=tol)
    fw .= euler_flux(U_star, γ)[1] .* dt .* len

    return nothing
end

function ustar(UL, UR, γ; tol=1e-6)
    # convert left/right conservative variables to primitives
    primL = conserve_prim(UL, γ)
    primR = conserve_prim(UR, γ)
    ρL, uL, _ = primL
    ρR, uR, _ = primR
    pL = 0.5 * primL[1] / primL[3]
    pR = 0.5 * primR[1] / primR[3]
    cL = sound_speed(primL, γ)
    cR = sound_speed(primR, γ)

    # define helper functions f(p) and its derivative for the Newton iteration
    f(p, p_i, ρ_i, c_i) =
        p > p_i ? (p - p_i) * sqrt(2.0 / ((γ + 1) * ρ_i) / (p + (γ - 1) / (γ + 1) * p_i)) :
        (2 * c_i / (γ - 1)) * ((p / p_i)^((γ - 1) / (2γ)) - 1)
    df(p, p_i, ρ_i, c_i) =
        p > p_i ?
        sqrt(2.0 / ((γ + 1) * ρ_i) / (p + (γ - 1) / (γ + 1) * p_i)) *
        (1 - 0.5 * (p - p_i) / (p + (γ - 1) / (γ + 1) * p_i)) :
        (1 / (ρ_i * c_i)) * ((p / p_i)^(-((γ + 1) / (2γ))))

    # initial guess for p_star (using the PVRS approximate solver)
    p_guess = 0.5 * (pL + pR) - 0.125 * (uR - uL) * (ρL + ρR) * (cL + cR)
    p_guess = max(1e-6, p_guess)

    # Newton iteration for p_star
    p_old = p_guess
    for _ in 1:100
        fL = f(p_old, pL, ρL, cL)
        fR = f(p_old, pR, ρR, cR)
        f_total = fL + fR + (uR - uL)
        dfL = df(p_old, pL, ρL, cL)
        dfR = df(p_old, pR, ρR, cR)
        dp = -f_total / (dfL + dfR)
        p_new = p_old + dp
        if abs(dp) < tol * p_old
            p_old = p_new
            break
        end
        p_old = p_new
    end
    p_star = p_old

    # u_star
    u_star = 0.5 * (uL + uR) + 0.5 * (f(p_star, pR, ρR, cR) - f(p_star, pL, ρL, cL))

    # sample state at x/t = 0
    if u_star >= 0
        # left star state
        if p_star > pL
            ρ_star =
                ρL *
                ((p_star / pL + (γ - 1) / (γ + 1)) / ((γ - 1) / (γ + 1) * p_star / pL + 1))
        else
            ρ_star = ρL * (p_star / pL)^(1 / γ)
        end
        U_star = prim_conserve([ρ_star, u_star, 0.5 * ρ_star / p_star], γ)
    else
        # right star state
        if p_star > pR
            ρ_star =
                ρR *
                ((p_star / pR + (γ - 1) / (γ + 1)) / ((γ - 1) / (γ + 1) * p_star / pR + 1))
        else
            ρ_star = ρR * (p_star / pR)^(1 / γ)
        end
        U_star = prim_conserve([ρ_star, u_star, 0.5 * ρ_star / p_star], γ)
    end

    return U_star
end
