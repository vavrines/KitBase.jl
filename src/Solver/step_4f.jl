"""
$(SIGNATURES)

1D4F1V
"""
function step!(
    KS::T,
    cell::ControlVolume1D4F,
    faceL::Interface1D4F,
    faceR::Interface1D4F,
    dx,
    dt,
    RES,
    AVG,
    collision = :bgk::Symbol,
    isMHD = true::Bool,
) where {T<:AbstractSolverSet}

    #--- update conservative flow variables: step 1 ---#
    # w^n
    w_old = deepcopy(cell.w)
    prim_old = deepcopy(cell.prim)

    # flux -> w^{n+1}
    @. cell.w += (faceL.fw - faceR.fw) / dx
    cell.prim .= mixture_conserve_prim(cell.w, KS.gas.γ)

    # temperature protection
    if cell.prim[5, 1] < 0
        @warn ("ion temperature update is negative")
        cell.w .= w_old
        cell.prim .= prim_old
    elseif cell.prim[5, 2] < 0
        @warn ("electron temperature update is negative")
        cell.w .= w_old
        cell.prim .= prim_old
    end

    # source -> w^{n+1}
    if isMHD == false
        #=
        # DifferentialEquations.jl
        tau = get_tau(cell.prim, KS.gas.mi, KS.gas.ni, KS.gas.me, KS.gas.ne, KS.gas.Kn[1])
        for j in axes(wRan, 2)
        prob = ODEProblem( mixture_source,
            vcat(cell.w[1:5,j,1], cell.w[1:5,j,2]),
            dt,
            (tau[1], tau[2], KS.gas.mi, KS.gas.ni, KS.gas.me, KS.gas.ne, KS.gas.Kn[1], KS.gas.γ) )
        sol = solve(prob, Rosenbrock23())

        cell.w[1:5,j,1] .= sol[end][1:5]
        cell.w[1:5,j,2] .= sol[end][6:10]
        for k=1:2
        cell.prim[:,j,k] .= conserve_prim(cell.w[:,j,k], KS.gas.γ)
        end
        end
        =#

        # explicit
        tau = aap_hs_collision_time(
            cell.prim,
            KS.gas.mi,
            KS.gas.ni,
            KS.gas.me,
            KS.gas.ne,
            KS.gas.Kn[1],
        )
        mprim = aap_hs_prim(
            cell.prim,
            tau,
            KS.gas.mi,
            KS.gas.ni,
            KS.gas.me,
            KS.gas.ne,
            KS.gas.Kn[1],
        )
        mw = mixture_prim_conserve(mprim, KS.gas.γ)
        for k = 1:2
            @. cell.w[:, k] += (mw[:, k] - w_old[:, k]) * dt / tau[k]
        end
        cell.prim .= mixture_conserve_prim(cell.w, KS.gas.γ)
    end

    #--- update electromagnetic variables ---#
    # flux -> E^{n+1} & B^{n+1}
    cell.E[1] -= dt * (faceL.femR[1] + faceR.femL[1]) / dx
    cell.E[2] -= dt * (faceL.femR[2] + faceR.femL[2]) / dx
    cell.E[3] -= dt * (faceL.femR[3] + faceR.femL[3]) / dx
    cell.B[1] -= dt * (faceL.femR[4] + faceR.femL[4]) / dx
    cell.B[2] -= dt * (faceL.femR[5] + faceR.femL[5]) / dx
    cell.B[3] -= dt * (faceL.femR[6] + faceR.femL[6]) / dx
    cell.ϕ -= dt * (faceL.femR[7] + faceR.femL[7]) / dx
    cell.ψ -= dt * (faceL.femR[8] + faceR.femL[8]) / dx

    for i = 1:3
        if 1 ∈ vcat(isnan.(cell.E), isnan.(cell.B))
            @warn "NaN electromagnetic update"
        end
    end

    # source -> ϕ
    #@. cell.ϕ += dt * (cell.w[1,:,1] / KS.gas.mi - cell.w[1,:,2] / KS.gas.me) / (KS.gas.lD^2 * KS.gas.rL)

    # source -> U^{n+1}, E^{n+1} and B^{n+1}
    mr = KS.gas.mi / KS.gas.me
    A, b = em_coefficients(cell.prim, cell.E, cell.B, mr, KS.gas.lD, KS.gas.rL, dt)
    x = A \ b

    #--- calculate lorenz force ---#
    cell.lorenz[1, 1] =
        0.5 * (
            x[1] + cell.E[1] + (cell.prim[3, 1] + x[5]) * cell.B[3] -
            (cell.prim[4, 1] + x[6]) * cell.B[2]
        ) / KS.gas.rL
    cell.lorenz[2, 1] =
        0.5 * (
            x[2] + cell.E[2] + (cell.prim[4, 1] + x[6]) * cell.B[1] -
            (cell.prim[2, 1] + x[4]) * cell.B[3]
        ) / KS.gas.rL
    cell.lorenz[3, 1] =
        0.5 * (
            x[3] + cell.E[3] + (cell.prim[2, 1] + x[4]) * cell.B[2] -
            (cell.prim[3, 1] + x[5]) * cell.B[1]
        ) / KS.gas.rL
    cell.lorenz[1, 2] =
        -0.5 *
        (
            x[1] + cell.E[1] + (cell.prim[3, 2] + x[8]) * cell.B[3] -
            (cell.prim[4, 2] + x[9]) * cell.B[2]
        ) *
        mr / KS.gas.rL
    cell.lorenz[2, 2] =
        -0.5 *
        (
            x[2] + cell.E[2] + (cell.prim[4, 2] + x[9]) * cell.B[1] -
            (cell.prim[2, 2] + x[7]) * cell.B[3]
        ) *
        mr / KS.gas.rL
    cell.lorenz[3, 2] =
        -0.5 *
        (
            x[3] + cell.E[3] + (cell.prim[2, 2] + x[7]) * cell.B[2] -
            (cell.prim[3, 2] + x[8]) * cell.B[1]
        ) *
        mr / KS.gas.rL

    cell.E[1] = x[1]
    cell.E[2] = x[2]
    cell.E[3] = x[3]

    #--- update conservative flow variables: step 2 ---#
    cell.prim[2, 1] = x[4]
    cell.prim[3, 1] = x[5]
    cell.prim[4, 1] = x[6]
    cell.prim[2, 2] = x[7]
    cell.prim[3, 2] = x[8]
    cell.prim[4, 2] = x[9]

    cell.w .= mixture_prim_conserve(cell.prim, KS.gas.γ)

    #--- update particle distribution function ---#
    # flux -> f^{n+1}
    @. cell.h0 += (faceL.fh0 - faceR.fh0) / dx
    @. cell.h1 += (faceL.fh1 - faceR.fh1) / dx
    @. cell.h2 += (faceL.fh2 - faceR.fh2) / dx
    @. cell.h3 += (faceL.fh3 - faceR.fh3) / dx

    # force -> f^{n+1} : step 1
    for j in axes(cell.h0, 2)
        _h0 = @view cell.h0[:, j]
        _h1 = @view cell.h1[:, j]
        _h2 = @view cell.h2[:, j]
        _h3 = @view cell.h3[:, j]

        shift_pdf!(_h0, cell.lorenz[1, j], KS.vs.du[1, j], dt)
        shift_pdf!(_h1, cell.lorenz[1, j], KS.vs.du[1, j], dt)
        shift_pdf!(_h2, cell.lorenz[1, j], KS.vs.du[1, j], dt)
        shift_pdf!(_h3, cell.lorenz[1, j], KS.vs.du[1, j], dt)
    end

    # force -> f^{n+1} : step 2
    for k in axes(cell.h1, 3)
        @. cell.h3[:, k] +=
            2.0 * dt * cell.lorenz[2, k] * cell.h1[:, k] +
            (dt * cell.lorenz[2, k])^2 * cell.h0[:, k] +
            2.0 * dt * cell.lorenz[3, k] * cell.h2[:, k] +
            (dt * cell.lorenz[3, k])^2 * cell.h0[:, k]
        @. cell.h2[:, k] += dt * cell.lorenz[3, k] * cell.h0[:, k]
        @. cell.h1[:, k] += dt * cell.lorenz[2, k] * cell.h0[:, k]
    end

    # source -> f^{n+1}
    tau = aap_hs_collision_time(
        cell.prim,
        KS.gas.mi,
        KS.gas.ni,
        KS.gas.me,
        KS.gas.ne,
        KS.gas.Kn[1],
    )

    # interspecies interaction
    if isMHD == true
        prim = deepcopy(cell.prim)
    else
        prim = aap_hs_prim(
            cell.prim,
            tau,
            KS.gas.mi,
            KS.gas.ni,
            KS.gas.me,
            KS.gas.ne,
            KS.gas.Kn[1],
        )
    end
    g = mixture_maxwellian(KS.vs.u, prim)

    # BGK term
    Mu, Mv, Mw, MuL, MuR = mixture_gauss_moments(prim, KS.gas.K)
    for k in axes(cell.h0, 2)
        @. cell.h0[:, k] = (cell.h0[:, k] + dt / tau[k] * g[:, k]) / (1.0 + dt / tau[k])
        @. cell.h1[:, k] =
            (cell.h1[:, k] + dt / tau[k] * Mv[1, k] * g[:, k]) / (1.0 + dt / tau[k])
        @. cell.h2[:, k] =
            (cell.h2[:, k] + dt / tau[k] * Mw[1, k] * g[:, k]) / (1.0 + dt / tau[k])
        @. cell.h3[:, k] =
            (cell.h3[:, k] + dt / tau[k] * (Mv[2, k] + Mw[2, k]) * g[:, k]) /
            (1.0 + dt / tau[k])
    end

    #--- record residuals ---#
    @. RES += (w_old - cell.w)^2
    @. AVG += abs(cell.w)

end
