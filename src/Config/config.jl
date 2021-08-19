# ============================================================
# Initial & Boundary Conditions of Specific Problems
# ============================================================

export ib_rh,
       ib_sod,
       ib_briowu,
       ib_cavity

"""
    1d0f: ib_rh(MaL, gam)
    1d1f1v: ib_rh(MaL, gam, u::T) where {T<:AbstractArray{<:AbstractFloat,1}}
    1d2f1v: ib_rh(MaL, gam, u::T, K) where {T<:AbstractArray{<:AbstractFloat,1}}
    1d1f3v: ib_rh(MaL, gam, u::T, v::T, w::T) where {T<:AbstractArray{<:AbstractFloat,3}}

Initialize Rankine-Hugoniot relation
"""
function ib_rh(
    set::AbstractSetup,
    ps::AbstractPhysicalSpace,
    vs::Union{AbstractVelocitySpace,Nothing},
    gas::AbstractProperty,
)

    gam = gas.γ
    MaL = gas.Ma
    MaR = sqrt((MaL^2 * (gam - 1.0) + 2.0) / (2.0 * gam * MaL^2 - (gam - 1.0)))
    ratioT =
        (1.0 + (gam - 1.0) / 2.0 * MaL^2) * (2.0 * gam / (gam - 1.0) * MaL^2 - 1.0) /
        (MaL^2 * (2.0 * gam / (gam - 1.0) + (gam - 1.0) / 2.0))

    if set.nSpecies == 1

        primL = [1.0, MaL * sqrt(gam / 2.0), 1.0]
        primR = [
            primL[1] * (gam + 1.0) * MaL^2 / ((gam - 1.0) * MaL^2 + 2.0),
            MaR * sqrt(gam / 2.0) * sqrt(ratioT),
            primL[3] / ratioT,
        ]

        primL = ifelse(set.space[5:6] == "3v", [primL[1:2]; zeros(2); primL[end]], primL)
        primR = ifelse(set.space[5:6] == "3v", [primR[1:2]; zeros(2); primR[end]], primR)

        wL = prim_conserve(primL, gam)
        wR = prim_conserve(primR, gam)

        fw = function(x)
            if x <= (ps.x0 + ps.x1) / 2
                return wL
            else
                return wR
            end
        end

        bc = function(x)
            if x <= (ps.x0 + ps.x1) / 2
                return primL
            else
                return primR
            end
        end

        if set.space[1:4] == "1d0f"
            return fw, bc
        elseif set.space == "1d1f1v"
            ff = function(x)
                w = fw(x)
                prim = conserve_prim(w, gas.γ)
                h = maxwellian(vs.u, prim)
                return h
            end

            return fw, ff, bc
        elseif set.space == "1d2f1v"
            ff = function(x)
                w = fw(x)
                prim = conserve_prim(w, gas.γ)
                h = maxwellian(vs.u, prim)
                b = h * gas.K / 2 / prim[end]
                return h, b
            end

            return fw, ff, bc
        elseif set.space == "1d1f3v"
            ff = function(x)
                w = fw(x)
                prim = conserve_prim(w, gas.γ)
                h = maxwellian(vs.u, vs.v, vs.w, prim)
                return h
            end

            return fw, ff, bc
        end

    elseif set.nSpecies == 2

        primL = zeros(3, 2)
        primL[:, 1] .= [mi, MaL * sqrt(gam / 2.0), mi / 1.0]
        primL[:, 2] .= [me, MaL * sqrt(gam / 2.0), me / 1.0]

        primR = zeros(3, 2)
        for j = 1:2
            primR[1, j] = primL[1, j] * (gam + 1.0) * MaL^2 / ((gam - 1.0) * MaL^2 + 2.0)
            primR[2, j] = MaR * sqrt(gam / 2.0 / (mi * ni + me * ne)) * sqrt(ratioT)
            primR[3, j] = primL[3, j] / ratioT
        end

        wL = mixture_prim_conserve(primL, gam)
        wR = mixture_prim_conserve(primR, gam)

        fw = function(x)
            if x <= (ps.x0 + ps.x1) / 2
                return wL
            else
                return wR
            end
        end

        bc = function(x)
            if x <= (ps.x0 + ps.x1) / 2
                return primL
            else
                return primR
            end
        end

        if set.space == "1d2f1v"
            hL = mixture_maxwellian(vs.u, primL)
            hR = mixture_maxwellian(vs.u, primR)
            bL = similar(hL)
            bR = similar(hR)
            for j = 1:2
                bL[:, j] .= hL[:, j] .* gas.K ./ (2.0 .* primL[end, j])
                bR[:, j] .= hR[:, j] .* gas.K ./ (2.0 .* primR[end, j])
            end

            ff = function(x)
                if x <= (ps.x0 + ps.x1) / 2
                    return hL, bL
                else
                    return hR, bR
                end
            end

            return fw, ff, bc
        end

    end

    return nothing

end

# ------------------------------------------------------------
# 1D0F
# ------------------------------------------------------------
function ib_rh(MaL, gam)
    primL = [1.0, MaL * sqrt(gam / 2.0), 1.0]

    MaR = sqrt((MaL^2 * (gam - 1.0) + 2.0) / (2.0 * gam * MaL^2 - (gam - 1.0)))
    ratioT =
        (1.0 + (gam - 1.0) / 2.0 * MaL^2) * (2.0 * gam / (gam - 1.0) * MaL^2 - 1.0) /
        (MaL^2 * (2.0 * gam / (gam - 1.0) + (gam - 1.0) / 2.0))

    primR = [
        primL[1] * (gam + 1.0) * MaL^2 / ((gam - 1.0) * MaL^2 + 2.0),
        MaR * sqrt(gam / 2.0) * sqrt(ratioT),
        primL[3] / ratioT,
    ]

    wL = prim_conserve(primL, gam)
    wR = prim_conserve(primR, gam)

    bcL = deepcopy(primL)
    bcR = deepcopy(primR)

    return wL, primL, bcL, wR, primR, bcR
end

# ------------------------------------------------------------
# 1D1F1V
# ------------------------------------------------------------
function ib_rh(MaL, gam, u::T) where {T<:AbstractArray{<:AbstractFloat,1}}
    primL = [1.0, MaL * sqrt(gam / 2.0), 1.0]

    MaR = sqrt((MaL^2 * (gam - 1.0) + 2.0) / (2.0 * gam * MaL^2 - (gam - 1.0)))
    ratioT =
        (1.0 + (gam - 1.0) / 2.0 * MaL^2) * (2.0 * gam / (gam - 1.0) * MaL^2 - 1.0) /
        (MaL^2 * (2.0 * gam / (gam - 1.0) + (gam - 1.0) / 2.0))

    primR = [
        primL[1] * (gam + 1.0) * MaL^2 / ((gam - 1.0) * MaL^2 + 2.0),
        MaR * sqrt(gam / 2.0) * sqrt(ratioT),
        primL[3] / ratioT,
    ]

    wL = prim_conserve(primL, gam)
    wR = prim_conserve(primR, gam)

    hL = maxwellian(u, primL)
    hR = maxwellian(u, primR)

    bcL = deepcopy(primL)
    bcR = deepcopy(primR)

    return wL, primL, hL, bcL, wR, primR, hR, bcR
end

# ------------------------------------------------------------
# 1D2F1V
# ------------------------------------------------------------
function ib_rh(MaL, gam, K, u::T) where {T<:AbstractArray{<:AbstractFloat,1}}
    primL = [1.0, MaL * sqrt(gam / 2.0), 1.0]

    MaR = sqrt((MaL^2 * (gam - 1.0) + 2.0) / (2.0 * gam * MaL^2 - (gam - 1.0)))
    ratioT =
        (1.0 + (gam - 1.0) / 2.0 * MaL^2) * (2.0 * gam / (gam - 1.0) * MaL^2 - 1.0) /
        (MaL^2 * (2.0 * gam / (gam - 1.0) + (gam - 1.0) / 2.0))

    primR = [
        primL[1] * (gam + 1.0) * MaL^2 / ((gam - 1.0) * MaL^2 + 2.0),
        MaR * sqrt(gam / 2.0) * sqrt(ratioT),
        primL[3] / ratioT,
    ]

    wL = prim_conserve(primL, gam)
    wR = prim_conserve(primR, gam)

    hL = maxwellian(u, primL)
    hR = maxwellian(u, primR)

    bL = hL .* K ./ (2.0 .* primL[end])
    bR = hR .* K ./ (2.0 .* primR[end])

    bcL = deepcopy(primL)
    bcR = deepcopy(primR)

    return wL, primL, hL, bL, bcL, wR, primR, hR, bR, bcR
end

# ------------------------------------------------------------
# 1D1F3V
# ------------------------------------------------------------
function ib_rh(MaL, gam, u::T, v::T, w::T) where {T<:AbstractArray{<:AbstractFloat,3}}
    primL = [1.0, MaL * sqrt(gam / 2.0), 0.0, 0.0, 1.0]

    MaR = sqrt((MaL^2 * (gam - 1.0) + 2.0) / (2.0 * gam * MaL^2 - (gam - 1.0)))
    ratioT =
        (1.0 + (gam - 1.0) / 2.0 * MaL^2) * (2.0 * gam / (gam - 1.0) * MaL^2 - 1.0) /
        (MaL^2 * (2.0 * gam / (gam - 1.0) + (gam - 1.0) / 2.0))

    primR = [
        primL[1] * (gam + 1.0) * MaL^2 / ((gam - 1.0) * MaL^2 + 2.0),
        MaR * sqrt(gam / 2.0) * sqrt(ratioT),
        0.0,
        0.0,
        primL[end] / ratioT,
    ]

    wL = prim_conserve(primL, gam)
    wR = prim_conserve(primR, gam)

    fL = maxwellian(u, v, w, primL)
    fR = maxwellian(u, v, w, primR)

    bcL = deepcopy(primL)
    bcR = deepcopy(primR)

    return wL, primL, fL, bcL, wR, primR, fR, bcR
end

# ------------------------------------------------------------
# 1D2F1V (mixture)
# ------------------------------------------------------------
function ib_rh(
    MaL,
    gam,
    K,
    mi,
    me,
    ni,
    ne,
    u::T,
) where {T<:AbstractArray{<:AbstractFloat,2}}
    MaR = sqrt((MaL^2 * (gam - 1.0) + 2.0) / (2.0 * gam * MaL^2 - (gam - 1.0)))
    ratioT =
        (1.0 + (gam - 1.0) / 2.0 * MaL^2) * (2.0 * gam / (gam - 1.0) * MaL^2 - 1.0) /
        (MaL^2 * (2.0 * gam / (gam - 1.0) + (gam - 1.0) / 2.0))

    primL = zeros(3, 2)
    primL[:, 1] .= [mi, MaL * sqrt(gam / 2.0), mi / 1.0]
    primL[:, 2] .= [me, MaL * sqrt(gam / 2.0), me / 1.0]

    primR = zeros(3, 2)
    for j = 1:2
        primR[1, j] = primL[1, j] * (gam + 1.0) * MaL^2 / ((gam - 1.0) * MaL^2 + 2.0)
        primR[2, j] = MaR * sqrt(gam / 2.0 / (mi * ni + me * ne)) * sqrt(ratioT)
        primR[3, j] = primL[3, j] / ratioT
    end

    wL = mixture_prim_conserve(primL, gam)
    wR = mixture_prim_conserve(primR, gam)

    hL = mixture_maxwellian(u, primL)
    hR = mixture_maxwellian(u, primR)

    bL = similar(hL)
    bR = similar(hR)
    for j = 1:2
        bL[:, j] = hL[:, j] .* K ./ (2.0 .* primL[end, j])
        bR[:, j] = hR[:, j] .* K ./ (2.0 .* primR[end, j])
    end

    bcL = deepcopy(primL)
    bcR = deepcopy(primR)

    return wL, primL, hL, bL, bcL, wR, primR, hR, bR, bcR
end


"""
    1d0f0v: ib_sod(γ)
    1d1f1v: ib_sod(γ, u::T) where {T<:AbstractArray{<:AbstractFloat,1}}
    1d1f3v: ib_sod(γ, u::T, v::T, w::T) where {T<:AbstractArray{<:AbstractFloat,3}}
    1d2f1v: ib_sod(γ, u::T, K) where {T<:AbstractArray{<:AbstractFloat,1}}

Initialize Sod shock tube
"""
function ib_sod(
    set::AbstractSetup,
    ps::AbstractPhysicalSpace,
    vs::Union{AbstractVelocitySpace,Nothing},
    gas::AbstractProperty,
)

    if set.nSpecies == 1

        primL = [1.0, 0.0, 0.5]
        primR = [0.125, 0.0, 0.625]

        wL = prim_conserve(primL, gas.γ)
        wR = prim_conserve(primR, gas.γ)

        fw = function(x)
            if x <= (ps.x0 + ps.x1) / 2
                return wL
            else
                return wR
            end
        end

        bc = function(x)
            if x <= (ps.x0 + ps.x1) / 2
                return primL
            else
                return primR
            end
        end

        if set.space[1:4] == "1d0f"
            return fw, bc
        elseif set.space == "1d1f1v"
            ff = function(x)
                w = fw(x)
                prim = conserve_prim(w, gas.γ)
                h = maxwellian(vs.u, prim)
                return h
            end

            return fw, ff, bc
        elseif set.space == "1d2f1v"
            ff = function(x)
                w = fw(x)
                prim = conserve_prim(w, gas.γ)
                h = maxwellian(vs.u, prim)
                b = h * gas.K / 2 / prim[end]
                return h, b
            end

            return fw, ff, bc
        elseif set.space == "1d1f3v"
            ff = function(x)
                w = fw(x)
                prim = conserve_prim(w, gas.γ)
                h = maxwellian(vs.u, vs.v, vs.w, prim)
                return h
            end

            return fw, ff, bc
        end

    elseif set.nSpecies == 2

        primL = zeros(3, 2)
        primL[:, 1] .= [gas.mi, 0.0, gas.mi / 2]
        primL[:, 2] .= [gas.me, 0.0, gas.me / 2]

        primR = zeros(3, 2)
        primR[:, 1] .= [gas.mi * 0.125, 0.0, gas.mi / 1.6]
        primR[:, 2] .= [gas.me * 0.125, 0.0, gas.me / 1.6]

        wL = mixture_prim_conserve(primL, gas.γ)
        wR = mixture_prim_conserve(primR, gas.γ)

        fw = function(x)
            if x <= (ps.x0 + ps.x1) / 2
                return wL
            else
                return wR
            end
        end

        bc = function(x)
            if x <= (ps.x0 + ps.x1) / 2
                return primL
            else
                return primR
            end
        end

        if set.space[1:4] == "1d0f"
            return fw, bc
        elseif set.space == "1d1f1v"
            hL = mixture_maxwellian(vs.u, primL)
            hR = mixture_maxwellian(vs.u, primR)

            ff = function(x)
                if x <= (ps.x0 + ps.x1) / 2
                    return hL
                else
                    return hR
                end
            end

            return fw, ff, bc
        elseif set.space == "1d2f1v"
            hL = mixture_maxwellian(vs.u, primL)
            hR = mixture_maxwellian(vs.u, primR)
            bL = similar(hL)
            bR = similar(hR)
            for j = 1:2
                bL[:, j] .= hL[:, j] .* K ./ (2.0 .* primL[end, j])
                bR[:, j] .= hR[:, j] .* K ./ (2.0 .* primR[end, j])
            end

            ff = function(x)
                if x <= (ps.x0 + ps.x1) / 2
                    return hL, bL
                else
                    return hR, bR
                end
            end

            return fw, ff, bc
        end

    end

    return nothing

end

function ib_sod(γ) # 1D0F0V

    primL = [1.0, 0.0, 0.5]
    primR = [0.125, 0.0, 0.625]

    wL = prim_conserve(primL, γ)
    wR = prim_conserve(primR, γ)

    bcL = deepcopy(primL)
    bcR = deepcopy(primR)

    return wL, primL, bcL, wR, primR, bcR

end

#--- mixture ---#
function ib_sod(γ::Real, mi::Real, me::Real)
    primL = zeros(3, 2)
    primL[:, 1] .= [mi, 0.0, mi / 2]
    primL[:, 2] .= [me, 0.0, me / 2]

    primR = zeros(3, 2)
    primR[:, 1] .= [mi * 0.125, 0.0, mi / 1.6]
    primR[:, 2] .= [me * 0.125, 0.0, me / 1.6]

    wL = mixture_prim_conserve(primL, γ)
    wR = mixture_prim_conserve(primR, γ)

    bcL = deepcopy(primL)
    bcR = deepcopy(primR)

    return wL, primL, bcL, wR, primR, bcR
end

# ------------------------------------------------------------
# 1D1F1V
# ------------------------------------------------------------
function ib_sod(γ, u::T) where {T<:AbstractArray{<:AbstractFloat,1}}

    primL = [1.0, 0.0, 0.5]
    primR = [0.125, 0.0, 0.625]

    wL = prim_conserve(primL, γ)
    wR = prim_conserve(primR, γ)

    fL = maxwellian(u, primL)
    fR = maxwellian(u, primR)

    bcL = deepcopy(primL)
    bcR = deepcopy(primR)

    return wL, primL, fL, bcL, wR, primR, fR, bcR

end

# ------------------------------------------------------------
# 1D1F3V
# ------------------------------------------------------------
function ib_sod(γ, u::T, v::T, w::T) where {T<:AbstractArray{<:AbstractFloat,3}}

    primL = [1.0, 0.0, 0.0, 0.0, 0.5]
    primR = [0.125, 0.0, 0.0, 0.0, 0.625]

    wL = prim_conserve(primL, γ)
    wR = prim_conserve(primR, γ)

    fL = maxwellian(u, v, w, primL)
    fR = maxwellian(u, v, w, primR)

    bcL = deepcopy(primL)
    bcR = deepcopy(primR)

    return wL, primL, fL, bcL, wR, primR, fR, bcR

end

# ------------------------------------------------------------
# 1D2F1V
# ------------------------------------------------------------
function ib_sod(γ, K, u::T) where {T<:AbstractArray{<:AbstractFloat,1}}

    primL = [1.0, 0.0, 0.5]
    primR = [0.125, 0.0, 0.625]

    wL = prim_conserve(primL, γ)
    wR = prim_conserve(primR, γ)

    hL = maxwellian(u, primL)
    hR = maxwellian(u, primR)

    bL = hL .* K ./ (2.0 .* primL[end])
    bR = hR .* K ./ (2.0 .* primR[end])

    bcL = deepcopy(primL)
    bcR = deepcopy(primR)

    return wL, primL, hL, bL, bcL, wR, primR, hR, bR, bcR

end


"""
    2d0f0v: ib_cavity(gam, Um, Vm, Tm) where {T<:AbstractArray{<:AbstractFloat,2}}
    2d1f2v: ib_cavity(gam, Um, Vm, Tm, u::T, v::T) where {T<:AbstractArray{<:AbstractFloat,2}}
    2d2f2v: ib_cavity(gam, Um, Vm, Tm, u::T, v::T, K) where {T<:AbstractArray{<:AbstractFloat,2}}

Initialize lid-driven cavity
"""
function ib_cavity(
    set::AbstractSetup,
    ps::AbstractPhysicalSpace,
    vs::Union{AbstractVelocitySpace,Nothing},
    gas::AbstractProperty,
    Um = 0.15,
    Vm = 0.0,
    Tm = 1.0,
)

    if set.nSpecies == 1

        prim = [1.0, 0.0, 0.0, 1.0]
        w = prim_conserve(prim, gas.γ)
        h = maxwellian(vs.u, prim)
        b = h .* gas.K / 2.0 / prim[end]

        primU = [1.0, Um, Vm, Tm]
        primD = deepcopy(prim)
        primL = deepcopy(prim)
        primR = deepcopy(prim)

        fw = function(args...)
            return w
        end

        bc = function(x, y)
            if y == ps.y1
                return primU
            else
                return prim
            end
        end

        if set.space[1:4] == "2d0f"
            return fw, bc
        elseif set.space == "2d1f2v"
            ff = function(args...)
                return h
            end

            return fw, ff, bc
        elseif set.space == "2d2f2v"
            ff = function(args...)
                return h, b
            end

            return fw, ff, bc
        end

    end

    return nothing

end

function ib_cavity(gam, Um, Vm, Tm) where {T<:AbstractArray{<:AbstractFloat,2}} # 2D0F

    primL = [1.0, 0.0, 0.0, 1.0]
    primR = deepcopy(primL)

    wL = prim_conserve(primL, gam)
    wR = prim_conserve(primR, gam)

    bcU = [1.0, Um, Vm, Tm]
    bcD = deepcopy(primR)
    bcL = deepcopy(primR)
    bcR = deepcopy(primR)

    return wL, primL, bcL, wR, primR, bcR, bcU, bcD

end


# ------------------------------------------------------------
# 2D1F2V
# ------------------------------------------------------------
function ib_cavity(gam, Um, Vm, Tm, u::T, v::T) where {T<:AbstractArray{<:AbstractFloat,2}}

    primL = [1.0, 0.0, 0.0, 1.0]
    primR = deepcopy(primL)

    wL = prim_conserve(primL, gam)
    wR = prim_conserve(primR, gam)

    fL = maxwellian(u, v, primL)
    fR = maxwellian(u, v, primR)

    bcU = [1.0, Um, Vm, Tm]
    bcD = deepcopy(primR)
    bcL = deepcopy(primR)
    bcR = deepcopy(primR)

    return wL, primL, fL, bcL, wR, primR, fR, bcR, bcU, bcD

end

# ------------------------------------------------------------
# 2D2F2V
# ------------------------------------------------------------
function ib_cavity(
    gam,
    K,
    Um,
    Vm,
    Tm,
    u::T,
    v::T,
) where {T<:AbstractArray{<:AbstractFloat,2}}

    primL = [1.0, 0.0, 0.0, 1.0]
    primR = deepcopy(primL)

    wL = prim_conserve(primL, gam)
    wR = prim_conserve(primR, gam)

    hL = maxwellian(u, v, primL)
    hR = maxwellian(u, v, primR)

    bL = hL .* K ./ (2.0 * primL[end])
    bR = hR .* K ./ (2.0 * primR[end])

    bcU = [1.0, Um, Vm, Tm]
    bcD = deepcopy(primR)
    bcL = deepcopy(primR)
    bcR = deepcopy(primR)

    return wL, primL, hL, bL, bcL, wR, primR, hR, bR, bcR, bcU, bcD

end


"""
    ib_briowu(gam, uspace::T, mi, me) where {T<:AbstractArray{<:AbstractFloat,2}}

Initialize Brio-Wu MHD shock tube

"""
function ib_briowu(
    set::AbstractSetup,
    ps::AbstractPhysicalSpace,
    vs::Union{AbstractVelocitySpace,Nothing},
    gas::AbstractProperty,
)

    # upstream
    primL = zeros(5, 2)
    primL[1, 1] = 1.0 * gas.mi
    primL[2, 1] = 0.0
    primL[3, 1] = 0.0
    primL[4, 1] = 0.0
    primL[5, 1] = gas.mi / 1.0
    primL[1, 2] = 1.0 * gas.me
    primL[2, 2] = 0.0
    primL[3, 2] = 0.0
    primL[4, 2] = 0.0
    primL[5, 2] = gas.me / 1.0

    wL = mixture_prim_conserve(primL, gas.γ)

    EL = zeros(3)
    BL = zeros(3)
    BL[1] = 0.75
    BL[2] = 1.0
    lorenzL = zeros(3, 2)

    # downstream
    primR = zeros(5, 2)
    primR[1, 1] = 0.125 * gas.mi
    primR[2, 1] = 0.0
    primR[3, 1] = 0.0
    primR[4, 1] = 0.0
    primR[5, 1] = gas.mi * 1.25
    primR[1, 2] = 0.125 * gas.me
    primR[2, 2] = 0.0
    primR[3, 2] = 0.0
    primR[4, 2] = 0.0
    primR[5, 2] = gas.me * 1.25

    wR = mixture_prim_conserve(primR, gas.γ)
    ER = zeros(3)
    BR = zeros(3)
    BR[1] = 0.75
    BR[2] = -1.0
    lorenzR = zeros(3, 2)

    fw = function(x)
        if x <= (ps.x0 + ps.x1) / 2
            return wL
        else
            return wR
        end
    end
    fE = function(x)
        if x <= (ps.x0 + ps.x1) / 2
            return EL
        else
            return ER
        end
    end
    fB = function(x)
        if x <= (ps.x0 + ps.x1) / 2
            return BL
        else
            return BR
        end
    end
    fL = function(x)
        if x <= (ps.x0 + ps.x1) / 2
            return lorenzL
        else
            return lorenzR
        end
    end

    bc = function(x)
        if x <= (ps.x0 + ps.x1) / 2
            return primL
        else
            return primR
        end
    end

    if set.space[3:end] == "4f1v"
        h0L = mixture_maxwellian(vs.u, primL)
        h1L = similar(h0L)
        h2L = similar(h0L)
        h3L = similar(h0L)
        for j in axes(h0L, 2)
            h1L[:, j] .= primL[3, j] .* h0L[:, j]
            h2L[:, j] .= primL[4, j] .* h0L[:, j]
            h3L[:, j] .=
                (primL[3, j]^2 + primL[4, j]^2 + 2.0 / (2.0 * primL[end, j])) .* h0L[:, j]
        end

        h0R = mixture_maxwellian(vs.u, primR)
        h1R = similar(h0R)
        h2R = similar(h0R)
        h3R = similar(h0R)
        for j in axes(h0L, 2)
            h1R[:, j] .= primR[3, j] .* h0R[:, j]
            h2R[:, j] .= primR[4, j] .* h0R[:, j]
            h3R[:, j] .=
                (primR[3, j]^2 + primR[4, j]^2 + 2.0 / (2.0 * primR[end, j])) .* h0R[:, j]
        end

        ff = function(x)
            if x <= (ps.x0 + ps.x1) / 2
                return h0L, h1L, h2L, h3L
            else
                return h0R, h1R, h2R, h3R
            end
        end

        return fw, ff, fE, fB, fL, bc
    elseif set.space[3:end] == "3f2v"
        h0L = mixture_maxwellian(vs.u, vs.v, primL)
        h1L = similar(h0L)
        h2L = similar(h0L)
        for j in axes(h0L, 3)
            h1L[:, :, j] .= primL[4, j] .* h0L[:, :, j]
            h2L[:, :, j] .= (primL[4, j]^2 + 1.0 / (2.0 * primL[end, j])) .* h0L[:, :, j]
        end

        h0R = mixture_maxwellian(vs.u, vs.v, primR)
        h1R = similar(h0R)
        h2R = similar(h0R)
        for j in axes(h0R, 3)
            h1R[:, :, j] .= primR[4, j] .* h0R[:, :, j]
            h2R[:, :, j] .= (primR[4, j]^2 + 1.0 / (2.0 * primR[end, j])) .* h0R[:, :, j]
        end

        ff = function(x)
            if x <= (ps.x0 + ps.x1) / 2
                return h0L, h1L, h2L
            else
                return h0R, h1R, h2R
            end
        end

        return fw, ff, fE, fB, fL, bc
    end

    return nothing

end

function ib_briowu(gam, mi, me, uspace::T) where {T<:AbstractArray{<:AbstractFloat,2}}

    # upstream
    primL = zeros(5, 2)
    primL[1, 1] = 1.0 * mi
    primL[2, 1] = 0.0
    primL[3, 1] = 0.0
    primL[4, 1] = 0.0
    primL[5, 1] = mi / 1.0
    primL[1, 2] = 1.0 * me
    primL[2, 2] = 0.0
    primL[3, 2] = 0.0
    primL[4, 2] = 0.0
    primL[5, 2] = me / 1.0

    wL = mixture_prim_conserve(primL, gam)
    h0L = mixture_maxwellian(uspace, primL)

    h1L = similar(h0L)
    h2L = similar(h0L)
    h3L = similar(h0L)
    for j in axes(h0L, 2)
        h1L[:, j] .= primL[3, j] .* h0L[:, j]
        h2L[:, j] .= primL[4, j] .* h0L[:, j]
        h3L[:, j] .=
            (primL[3, j]^2 + primL[4, j]^2 + 2.0 / (2.0 * primL[end, j])) .* h0L[:, j]
    end

    EL = zeros(3)
    BL = zeros(3)
    BL[1] = 0.75
    BL[2] = 1.0

    # downstream
    primR = zeros(5, 2)
    primR[1, 1] = 0.125 * mi
    primR[2, 1] = 0.0
    primR[3, 1] = 0.0
    primR[4, 1] = 0.0
    primR[5, 1] = mi * 1.25
    primR[1, 2] = 0.125 * me
    primR[2, 2] = 0.0
    primR[3, 2] = 0.0
    primR[4, 2] = 0.0
    primR[5, 2] = me * 1.25

    wR = mixture_prim_conserve(primR, gam)
    h0R = mixture_maxwellian(uspace, primR)

    h1R = similar(h0R)
    h2R = similar(h0R)
    h3R = similar(h0R)
    for j in axes(h0L, 2)
        h1R[:, j] .= primR[3, j] .* h0R[:, j]
        h2R[:, j] .= primR[4, j] .* h0R[:, j]
        h3R[:, j] .=
            (primR[3, j]^2 + primR[4, j]^2 + 2.0 / (2.0 * primR[end, j])) .* h0R[:, j]
    end

    ER = zeros(3)
    BR = zeros(3)
    BR[1] = 0.75
    BR[2] = -1.0

    lorenzL = zeros(3, 2)
    lorenzR = zeros(3, 2)
    bcL = zeros(5, 2)
    bcR = zeros(5, 2)

    return wL,
    primL,
    h0L,
    h1L,
    h2L,
    h3L,
    bcL,
    EL,
    BL,
    lorenzL,
    wR,
    primR,
    h0R,
    h1R,
    h2R,
    h3R,
    bcR,
    ER,
    BR,
    lorenzR

end

function ib_briowu(
    gam,
    mi,
    me,
    uspace::T,
    vspace::T,
) where {T<:AbstractArray{<:AbstractFloat,3}}

    # upstream
    primL = zeros(5, 2)
    primL[1, 1] = 1.0 * mi
    primL[2, 1] = 0.0
    primL[3, 1] = 0.0
    primL[4, 1] = 0.0
    primL[5, 1] = mi / 1.0
    primL[1, 2] = 1.0 * me
    primL[2, 2] = 0.0
    primL[3, 2] = 0.0
    primL[4, 2] = 0.0
    primL[5, 2] = me / 1.0

    wL = mixture_prim_conserve(primL, gam)
    h0L = mixture_maxwellian(uspace, vspace, primL)

    h1L = similar(h0L)
    h2L = similar(h0L)
    for j in axes(h0L, 3)
        h1L[:, :, j] .= primL[4, j] .* h0L[:, :, j]
        h2L[:, :, j] .= (primL[4, j]^2 + 1.0 / (2.0 * primL[end, j])) .* h0L[:, :, j]
    end

    EL = zeros(3)
    BL = zeros(3)
    BL[1] = 0.75
    BL[2] = 1.0

    # downstream
    primR = zeros(5, 2)
    primR[1, 1] = 0.125 * mi
    primR[2, 1] = 0.0
    primR[3, 1] = 0.0
    primR[4, 1] = 0.0
    primR[5, 1] = mi * 1.25
    primR[1, 2] = 0.125 * me
    primR[2, 2] = 0.0
    primR[3, 2] = 0.0
    primR[4, 2] = 0.0
    primR[5, 2] = me * 1.25

    wR = mixture_prim_conserve(primR, gam)
    h0R = mixture_maxwellian(uspace, vspace, primR)

    h1R = similar(h0R)
    h2R = similar(h0R)
    for j in axes(h0R, 3)
        h1R[:, :, j] .= primR[4, j] .* h0R[:, :, j]
        h2R[:, :, j] .= (primR[4, j]^2 + 1.0 / (2.0 * primR[end, j])) .* h0R[:, :, j]
    end

    ER = zeros(3)
    BR = zeros(3)
    BR[1] = 0.75
    BR[2] = -1.0

    lorenzL = zeros(3, 2)
    lorenzR = zeros(3, 2)
    bcL = deepcopy(primL)
    bcR = deepcopy(primR)

    return wL,
    primL,
    h0L,
    h1L,
    h2L,
    bcL,
    EL,
    BL,
    lorenzL,
    wR,
    primR,
    h0R,
    h1R,
    h2R,
    bcR,
    ER,
    BR,
    lorenzR

end
