"""
$(SIGNATURES)

1D3F1V Rykov
"""
function step!(
    w::T3,
    prim::T3,
    h::T4,
    b::T4,
    r::T4,
    fwL::T1,
    fhL::T2,
    fbL::T2,
    frL::T2,
    fwR::T1,
    fhR::T2,
    fbR::T2,
    frR::T2,
    u::T5,
    weights::T5,
    K,
    Kr,
    μᵣ,
    ω,
    Pr,
    T₀,
    Z₀,
    σ,
    ω0,
    ω1,
    dx,
    dt,
    RES,
    AVG,
    collision,
) where {T1<:AA{<:FN,1},T2<:AA{<:FN,1},T3<:AA{<:Real,1},T4<:AA{<:FN,1},T5<:AA{<:FN,1}}

    @assert collision == :rykov

    #--- store W^n and calculate shakhov term ---#
    w_old = deepcopy(w)

    if collision == :rykov
        q = heat_flux(h, b, r, prim, u, weights)
    else
        q = zeros(2)
    end

    #--- update W^{n+1} ---#
    @. w += (fwL - fwR) / dx

    MHT = similar(h)
    MBT = similar(b)
    MRT = similar(r)
    MHR = similar(h)
    MBR = similar(b)
    MRR = similar(r)
    maxwellian!(MHT, MBT, MRT, MHR, MBR, MRR, u, prim, K, Kr)
    τ_old = vhs_collision_time(prim[1:end-1], μᵣ, ω)
    Zr = rykov_zr(1.0 / prim[4], T₀, Z₀)
    Er0_old = 0.5 * sum(@. weights * ((1.0 / Zr) * MRR + (1.0 - 1.0 / Zr) * MRT))

    w[4] += dt * (Er0_old - w_old[4]) / τ_old
    prim .= conserve_prim(w, K, Kr)

    #--- record residuals ---#
    @. RES += (w - w_old)^2
    @. AVG += abs(w)

    #--- calculate M^{n+1} and tau^{n+1} ---#
    maxwellian!(MHT, MBT, MRT, MHR, MBR, MRR, u, prim, K, Kr)

    SHT = similar(h)
    SBT = similar(b)
    SRT = similar(r)
    SHR = similar(h)
    SBR = similar(b)
    SRR = similar(r)
    rykov!(
        SHT,
        SBT,
        SRT,
        SHR,
        SBR,
        SRR,
        u,
        MHT,
        MBT,
        MRT,
        MHR,
        MBR,
        MRR,
        q,
        prim,
        Pr,
        K,
        σ,
        ω0,
        ω1,
    )

    MH = (1.0 - 1.0 / Zr) * (MHT + SHT) + 1.0 / Zr * (MHR + SHR)
    MB = (1.0 - 1.0 / Zr) * (MBT + SBT) + 1.0 / Zr * (MBR + SBR)
    MR = (1.0 - 1.0 / Zr) * (MRT + SRT) + 1.0 / Zr * (MRR + SRR)

    τ = vhs_collision_time(prim[1:end-1], μᵣ, ω)

    #--- update distribution function ---#
    for i in eachindex(u)
        h[i] = (h[i] + (fhL[i] - fhR[i]) / dx + dt / τ * MH[i]) / (1.0 + dt / τ)
        b[i] = (b[i] + (fbL[i] - fbR[i]) / dx + dt / τ * MB[i]) / (1.0 + dt / τ)
        r[i] = (r[i] + (frL[i] - frR[i]) / dx + dt / τ * MR[i]) / (1.0 + dt / τ)
    end

end

"""
$(SIGNATURES)

1D3F2V
"""
function step!(
    KS::T,
    cell::ControlVolume1D3F,
    faceL::Interface1D3F,
    faceR::Interface1D3F,
    dx,
    dt,
    RES,
    AVG,
    collision,
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
    if cell.prim[end, 1] < 0
        @warn ("ion temperature update is negative")
        cell.w .= w_old
        cell.prim .= prim_old
    elseif cell.prim[end, 2] < 0
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
        for k in axes(cell.w, 2)
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
            @warn "electromagnetic update is NaN"
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

    # force -> f^{n+1} : step 1
    for j in axes(cell.h0, 3) # component
        for i in axes(cell.h0, 2) # v
            _h0 = @view cell.h0[:, i, j]
            _h1 = @view cell.h1[:, i, j]
            _h2 = @view cell.h2[:, i, j]

            shift_pdf!(_h0, cell.lorenz[1, j], KS.vs.du[1, i, j], dt)
            shift_pdf!(_h1, cell.lorenz[1, j], KS.vs.du[1, i, j], dt)
            shift_pdf!(_h2, cell.lorenz[1, j], KS.vs.du[1, i, j], dt)
        end
    end

    for j in axes(cell.h0, 3) # component
        for i in axes(cell.h0, 1) # u
            _h0 = @view cell.h0[i, :, j]
            _h1 = @view cell.h1[i, :, j]
            _h2 = @view cell.h2[i, :, j]

            shift_pdf!(_h0, cell.lorenz[2, j], KS.vs.dv[i, 1, j], dt)
            shift_pdf!(_h1, cell.lorenz[2, j], KS.vs.dv[i, 1, j], dt)
            shift_pdf!(_h2, cell.lorenz[2, j], KS.vs.dv[i, 1, j], dt)
        end
    end

    # force -> f^{n+1} : step 2
    for k in axes(cell.h1, 3)
        @. cell.h2[:, :, k] +=
            2.0 * dt * cell.lorenz[3, k] * cell.h1[:, :, k] +
            (dt * cell.lorenz[3, k])^2 * cell.h0[:, :, k]
        @. cell.h1[:, :, k] += dt * cell.lorenz[3, k] * cell.h0[:, :, k]
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

    H0 = similar(KS.vs.u)
    H1 = similar(H0)
    H2 = similar(H0)
    for k in axes(H0, 3)
        H0[:, :, k] .= maxwellian(KS.vs.u[:, :, k], KS.vs.v[:, :, k], prim[:, k])
        @. H1[:, :, k] = H0[:, :, k] * prim[4, k]
        @. H2[:, :, k] = H0[:, :, k] * (prim[4, k]^2 + 1.0 / (2.0 * prim[5, k]))
    end

    # BGK term
    for k in axes(cell.h0, 3)
        @. cell.h0[:, :, k] =
            (cell.h0[:, :, k] + dt / tau[k] * H0[:, :, k]) / (1.0 + dt / tau[k])
        @. cell.h1[:, :, k] =
            (cell.h1[:, :, k] + dt / tau[k] * H1[:, :, k]) / (1.0 + dt / tau[k]) # NOTICE the h1 here is h2 in 1d4f case
        @. cell.h2[:, :, k] =
            (cell.h2[:, :, k] + dt / tau[k] * H2[:, :, k]) / (1.0 + dt / tau[k]) # NOTICE the h2 here is h3 in 1d4f case
    end

    #--- record residuals ---#
    @. RES += (w_old - cell.w)^2
    @. AVG += abs(cell.w)

end
