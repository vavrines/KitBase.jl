# ============================================================
# Methods of Particle Velocity Space
# ============================================================

"""
    struct VSpace1D{TR<:Real,TI<:Integer,TA<:AbstractArray,TB<:AbstractArray{<:Real,1}} <: AbstractVelocitySpace
        u0::TR
        u1::TR
        nu::TI
        u::TA
        du::TA
        weights::TB
    end

1D velocity space

"""
struct VSpace1D{TR<:Real,TI<:Integer,TA<:AbstractArray,TB<:AbstractArray{<:Real,1}} <: AbstractVelocitySpace
    u0::TR
    u1::TR
    nu::TI
    u::TA
    du::TA
    weights::TB
end

function VSpace1D(
    U0,
    U1,
    NU::TI,
    TYPE = "rectangle",
    NG = zero(NU)::TI,
    PRECISION = Float64,
) where {TI<:Integer}

    δ = (U1 - U0) / NU
    u = OffsetArray{PRECISION}(undef, 1-NG:NU+NG)
    du = similar(u)
    weights = similar(u)

    if TYPE == "rectangle" # rectangular
        for i in eachindex(u)
            u[i] = U0 + (i - 0.5) * δ
            du[i] = δ
            weights[i] = δ
        end
    elseif TYPE == "newton" # newton-cotes
        for i in eachindex(u)
            u[i] = U0 + (i - 0.5) * δ
            du[i] = δ
            weights[i] = newton_cotes(i + NG, NU + NG * 2) * δ
        end
    else
        throw("No velocity quadrature available")
    end

    return VSpace1D{PRECISION,TI,typeof(u),typeof(weights)}(U0, U1, NU, u, du, weights)

end

VSpace1D() = VSpace1D(-5, 5, 50)
VSpace1D(U0::T, U1::T) where {T<:Real} = VSpace1D(U0, U1, 50)


"""
    struct VSpace2D{TR<:Real,TI<:Integer,TA<:AbstractArray{<:Real,2}} <: AbstractVelocitySpace
        u0::TR
        u1::TR
        nu::TI
        v0::TR
        v1::TR
        nv::TI
        u::TA
        v::TA
        du::TA
        dv::TA
        weights::TA
    end

2D velocity space

"""
struct VSpace2D{TR<:Real,TI<:Integer,TA<:AbstractArray{<:Real,2}} <: AbstractVelocitySpace
    u0::TR
    u1::TR
    nu::TI
    v0::TR
    v1::TR
    nv::TI
    u::TA
    v::TA
    du::TA
    dv::TA
    weights::TA
end

function VSpace2D(
    U0,
    U1,
    NU::TI,
    V0,
    V1,
    NV::TI,
    TYPE = "rectangle",
    NGU = zero(NU)::TI,
    NGV = zero(NV)::TI,
    PRECISION = Float64,
) where {TI<:Integer}

    δu = (U1 - U0) / NU
    δv = (V1 - V0) / NV
    u = OffsetArray{PRECISION}(undef, 1-NGU:NU+NGU, 1-NGV:NV+NGV)
    v = similar(u)
    du = similar(u)
    dv = similar(u)
    weights = similar(u)

    if TYPE == "rectangle" # rectangular
        for j in axes(u, 2)
            for i in axes(u, 1)
                u[i, j] = U0 + (i - 0.5) * δu
                v[i, j] = V0 + (j - 0.5) * δv
                du[i, j] = δu
                dv[i, j] = δv
                weights[i, j] = δu * δv
            end
        end
    elseif TYPE == "newton" # newton-cotes
        for j in axes(u, 2)
            for i in axes(u, 1)
                u[i, j] = U0 + (i - 0.5) * δu
                v[i, j] = V0 + (j - 0.5) * δv
                du[i, j] = δu
                dv[i, j] = δv
                weights[i, j] =
                    newton_cotes(i + NGU, NU + NGU * 2) *
                    δu *
                    newton_cotes(j + NGV, NV + NGV * 2) *
                    δv
            end
        end
    else
        throw("No velocity quadrature available")
    end

    return VSpace2D{PRECISION,TI,typeof(u)}(U0, U1, NU, V0, V1, NV, u, v, du, dv, weights)

end

VSpace2D() = VSpace2D(-5, 5, 28, -5, 5, 28)
VSpace2D(U0::T, U1::T, V0::T, V1::T) where {T<:Real} = VSpace2D(U0, U1, 28, V0, V1, 28)


"""
    struct VSpace3D{TR<:Real,TI<:Integer,TA<:AbstractArray{<:Real,3}} <: AbstractVelocitySpace
        u0::TR
        u1::TR
        nu::TI
        v0::TR
        v1::TR
        nv::TI
        w0::TR
        w1::TR
        nw::TI
        u::TA
        v::TA
        w::TA
        du::TA
        dv::TA
        dw::TA
        weights::TA
    end

3D velocity space

"""
struct VSpace3D{TR<:Real,TI<:Integer,TA<:AbstractArray{<:Real,3}} <: AbstractVelocitySpace
    u0::TR
    u1::TR
    nu::TI
    v0::TR
    v1::TR
    nv::TI
    w0::TR
    w1::TR
    nw::TI
    u::TA
    v::TA
    w::TA
    du::TA
    dv::TA
    dw::TA
    weights::TA
end

function VSpace3D(
    U0,
    U1,
    NU::TI,
    V0,
    V1,
    NV::TI,
    W0,
    W1,
    NW::TI,
    TYPE = "rectangle",
    NGU = zero(NU)::TI,
    NGV = zero(NV)::TI,
    NGW = zero(NW)::TI,
    PRECISION = Float64,
) where {TI<:Integer}

    δu = (U1 - U0) / NU
    δv = (V1 - V0) / NV
    δw = (W1 - W0) / NW
    u = OffsetArray{PRECISION}(undef, 1-NGU:NU+NGU, 1-NGV:NV+NGV, 1-NGW:NW+NGW)
    v = similar(u)
    w = similar(u)
    du = similar(u)
    dv = similar(u)
    dw = similar(u)
    weights = similar(u)

    if TYPE == "rectangle" # rectangular
        for k in axes(u, 3), j in axes(u, 2), i in axes(u, 1)
            u[i, j, k] = U0 + (i - 0.5) * δu
            v[i, j, k] = V0 + (j - 0.5) * δv
            w[i, j, k] = W0 + (k - 0.5) * δw
            du[i, j, k] = δu
            dv[i, j, k] = δv
            dw[i, j, k] = δw
            weights[i, j, k] = δu * δv * δw
        end
    elseif TYPE == "newton" # newton-cotes
        for k in axes(u, 3), j in axes(u, 2), i in axes(u, 1)
            u[i, j, k] = U0 + (i - 0.5) * δu
            v[i, j, k] = V0 + (j - 0.5) * δv
            w[i, j, k] = W0 + (k - 0.5) * δw
            du[i, j, k] = δu
            dv[i, j, k] = δv
            dw[i, j, k] = δw
            weights[i, j, k] =
                newton_cotes(i + NGU, NU + NGU * 2) *
                δu *
                newton_cotes(j + NGV, NV + NGV * 2) *
                δv *
                newton_cotes(k + NGW, NW + NGW * 2) *
                δw
        end
    else
        throw("No velocity quadrature available")
    end

    return VSpace3D{PRECISION,TI,typeof(u)}(U0, U1, NU, V0, V1, NV, W0, W1, NW, u, v, w, du, dv, dw, weights)

end

VSpace3D() = VSpace3D(-5, 5, 28, -5, 5, 28, -5, 5, 28)
VSpace3D(U0::T, U1::T, V0::T, V1::T, W0::T, W1::T) where {T<:Real} =
    VSpace3D(U0, U1, 28, V0, V1, 28, W0, W1, 28)


"""
    struct MVSpace1D{TR<:AbstractArray{<:Real,1},TI<:Integer,TA<:AbstractArray{<:Real,2}} <: AbstractVelocitySpace
        u0::TR
        u1::TR
        nu::TI
        u::TA
        du::TA
        weights::TA
    end

1D multi-component velocity space

"""
struct MVSpace1D{TR<:AbstractArray{<:Real,1},TI<:Integer,TA<:AbstractArray{<:Real,2}} <: AbstractVelocitySpace
    u0::TR
    u1::TR
    nu::TI
    u::TA
    du::TA
    weights::TA
end

function MVSpace1D(
    Ui0,
    Ui1,
    Ue0,
    Ue1,
    NU::TI,
    TYPE = "rectangle",
    NG = zero(NU)::TI,
    PRECISION = Float64,
) where {TI<:Integer}

    u0 = PRECISION.([Ui0, Ue0])
    u1 = PRECISION.([Ui1, Ue1])
    δ = (u1 .- u0) ./ NU
    u = OffsetArray{PRECISION}(undef, 1-NG:NU+NG, 1:2)
    du = similar(u)
    weights = similar(u)

    if TYPE == "rectangle" # rectangular
        for j in axes(u, 2), i in axes(u, 1)
            u[i, j] = u0[j] + (i - 0.5) * δ[j]
            du[i, j] = δ[j]
            weights[i, j] = δ[j]
        end
    elseif TYPE == "newton" # newton-cotes
        for j in axes(u, 2), i in axes(u, 1)
            u[i, j] = u0[j] + (i - 0.5) * δ[j]
            du[i, j] = δ[j]
            weights[i, j] = newton_cotes(i + NG, NU + NG * 2) * δ[j]
        end
    else
        throw("No velocity quadrature available")
    end

    return MVSpace1D{typeof(u0),TI,typeof(u)}(u0, u1, NU, u, du, weights)

end

MVSpace1D() = MVSpace1D(-5, 5, -10, 10, 28)
MVSpace1D(U0::T, U1::T, V0::T, V1::T) where {T<:Real} = MVSpace1D(U0, U1, V0, V1, 28)


"""
    struct MVSpace2D{TR<:AbstractArray{<:Real,1},TI<:Integer,TA<:AbstractArray{Float64,3}} <: AbstractVelocitySpace
        u0::TR
        u1::TR
        nu::TI
        v0::TR
        v1::TR
        nv::TI
        u::TA
        v::TA
        du::TA
        dv::TA
        weights::TA
    end

2D multi-component velocity space

"""
struct MVSpace2D{TR<:AbstractArray{<:Real,1},TI<:Integer,TA<:AbstractArray{Float64,3}} <: AbstractVelocitySpace
    u0::TR
    u1::TR
    nu::TI
    v0::TR
    v1::TR
    nv::TI
    u::TA
    v::TA
    du::TA
    dv::TA
    weights::TA
end

function MVSpace2D(
    Ui0,
    Ui1,
    Ue0,
    Ue1,
    NU::TI,
    Vi0,
    Vi1,
    Ve0,
    Ve1,
    NV::TI,
    TYPE = "rectangle",
    NGU = zero(NU)::TI,
    NGV = zero(NV)::TI,
    PRECISION = Float64,
) where {TI<:Integer}

    u0 = PRECISION.([Ui0, Ue0])
    u1 = PRECISION.([Ui1, Ue1])
    δu = (u1 .- u0) ./ NU
    v0 = PRECISION.([Vi0, Ve0])
    v1 = PRECISION.([Vi1, Ve1])
    δv = (v1 .- v0) ./ NV
    u = OffsetArray{PRECISION}(undef, 1-NGU:NU+NGU, 1-NGV:NV+NGV, 1:2)
    v = similar(u)
    du = similar(u)
    dv = similar(u)
    weights = similar(u)

    if TYPE == "rectangle" # rectangular
        for k in axes(u, 3), j in axes(u, 2), i in axes(u, 1)
            u[i, j, k] = u0[k] + (i - 0.5) * δu[k]
            v[i, j, k] = v0[k] + (j - 0.5) * δv[k]
            du[i, j, k] = δu[k]
            dv[i, j, k] = δv[k]
            weights[i, j, k] = δu[k] * δv[k]
        end
    elseif TYPE == "newton" # newton-cotes
        for k in axes(u, 3), j in axes(u, 2), i in axes(u, 1)
            u[i, j, k] = u0[k] + (i - 0.5) * δu[k]
            v[i, j, k] = v0[k] + (j - 0.5) * δv[k]
            du[i, j, k] = δu[k]
            dv[i, j, k] = δv[k]
            weights[i, j, k] =
                newton_cotes(i + NGU, NU + NGU * 2) *
                δu[k] *
                newton_cotes(j + NGV, NV + NGV * 2) *
                δv[k]
        end
    else
        throw("No velocity quadrature available")
    end

    return MVSpace2D{typeof(u0),TI,typeof(u)}(u0, u1, NU, v0, v1, NV, u, v, du, dv, weights)

end

MVSpace2D() = MVSpace2D(-5, 5, -10, 10, 28, -5, 5, -10, 10, 28)
MVSpace2D(U0::T, U1::T, V0::T, V1::T) where {T<:Real} =
    MVSpace2D(U0, U1, U0, U1, 28, V0, V1, V0, V1, 28)


"""
    struct MVSpace3D{TR<:AbstractArray{<:Real,1},TI<:Integer,TA<:AbstractArray{<:Real,4}} <: AbstractVelocitySpace
        u0::TR
        u1::TR
        nu::TI
        v0::TR
        v1::TR
        nv::TI
        w0::TR
        w1::TR
        nw::TI
        u::TA
        v::TA
        w::TA
        du::TA
        dv::TA
        dw::TA
        weights::TA
    end

3D multi-component velocity space

"""
struct MVSpace3D{TR<:AbstractArray{<:Real,1},TI<:Integer,TA<:AbstractArray{<:Real,4}} <: AbstractVelocitySpace
    u0::TR
    u1::TR
    nu::TI
    v0::TR
    v1::TR
    nv::TI
    w0::TR
    w1::TR
    nw::TI
    u::TA
    v::TA
    w::TA
    du::TA
    dv::TA
    dw::TA
    weights::TA
end

function MVSpace3D(
    Ui0,
    Ui1,
    Ue0,
    Ue1,
    NU::TI,
    Vi0,
    Vi1,
    Ve0,
    Ve1,
    NV::TI,
    Wi0,
    Wi1,
    We0,
    We1,
    NW::TI,
    TYPE = "rectangle",
    NGU = zero(NU)::TI,
    NGV = zero(NV)::TI,
    NGW = zero(NW)::TI,
    PRECISION = Float64,
) where {TI<:Integer}

    u0 = PRECISION.([Ui0, Ue0])
    u1 = PRECISION.([Ui1, Ue1])
    δu = (u1 .- u0) ./ NU
    v0 = PRECISION.([Vi0, Ve0])
    v1 = PRECISION.([Vi1, Ve1])
    δv = (v1 .- v0) ./ NV
    w0 = PRECISION.([Wi0, We0])
    w1 = PRECISION.([Wi1, We1])
    δw = (w1 .- w0) ./ NW

    u = OffsetArray{PRECISION}(undef, 1-NGU:NU+NGU, 1-NGV:NV+NGV, 1-NGW:NW+NGW, 1:2)
    v = similar(u)
    w = similar(u)
    du = similar(u)
    dv = similar(u)
    dw = similar(u)
    weights = similar(u)

    if TYPE == "rectangle" # rectangular
        for l in axes(u, 4), k in axes(u, 3), j in axes(u, 2), i in axes(u, 1)
            u[i, j, k, l] = u0[l] + (i - 0.5) * δu[l]
            v[i, j, k, l] = v0[l] + (j - 0.5) * δv[l]
            w[i, j, k, l] = w0[l] + (k - 0.5) * δw[l]
            du[i, j, k, l] = δu[l]
            dv[i, j, k, l] = δv[l]
            dw[i, j, k, l] = δw[l]
            weights[i, j, k, l] = δu[l] * δv[l] * δw[l]
        end
    elseif TYPE == "newton" # newton-cotes
        for l in axes(u, 4), k in axes(u, 3), j in axes(u, 2), i in axes(u, 1)
            u[i, j, k, l] = u0[l] + (i - 0.5) * δu[l]
            v[i, j, k, l] = v0[l] + (j - 0.5) * δv[l]
            w[i, j, k, l] = w0[l] + (k - 0.5) * δw[l]
            du[i, j, k, l] = δu[l]
            dv[i, j, k, l] = δv[l]
            dw[i, j, k, l] = δw[l]
            weights[i, j, k, l] =
                newton_cotes(i + NGU, NU + NGU * 2) *
                δu[l] *
                newton_cotes(j + NGV, NV + NGV * 2) *
                δv[l] *
                newton_cotes(k + NGW, NW + NGW * 2) *
                δw[l]
        end
    else
        throw("No velocity quadrature available")
    end

    return MVSpace3D{typeof(u0),TI,typeof(u)}(u0, u1, NU, v0, v1, NV, w0, w1, NW, u, v, v, du, dv, dw, weights)

end

MVSpace3D() = MVSpace3D(-5, 5, -10, 10, 20, -5, 5, -10, 10, 20, -5, 5, -10, 10, 20)
MVSpace3D(U0::T, U1::T, V0::T, V1::T, W0::T, W1::T) where {T<:Real} =
    MVSpace3D(U0, U1, U0, U1, 20, V0, V1, V0, V1, 20, W0, W1, W0, W1, 20)


"""
    newton_cotes(idx::T, num::T) where {T<:Integer}

Evaluate quadrature weight from Newton-Cotes rule
"""
function newton_cotes(idx::T, num::T) where {T<:Integer}
    if idx == 1 || idx == num
        nc_coeff = 14.0 / 45.0
    elseif (idx - 5) % 4 == 0
        nc_coeff = 28.0 / 45.0
    elseif (idx - 3) % 4 == 0
        nc_coeff = 24.0 / 45.0
    else
        nc_coeff = 64.0 / 45.0
    end

    return nc_coeff
end

# ------------------------------------------------------------
# Extended Base.show()
# ------------------------------------------------------------
function Base.show(io::IO, vs::VSpace1D{TR,TI,TA,TB}) where {TR,TI,TA,TB}
    print(io, "VelocitySpace1D{$TR,$TI,$TA,$TB}\n",
              "domain: ($(vs.u0),$(vs.u1))\n",
              "resolution: $(vs.nu)\n",
              "ghost: $(1-firstindex(vs.u))\n")
end

function Base.show(io::IO, vs::VSpace2D{TR,TI,TA}) where {TR,TI,TA}
    print(io, "VelocitySpace2D{$TR,$TI,$TA}\n",
              "domain: ($(vs.u0),$(vs.u1)) × ($(vs.v0),$(vs.v1))\n",
              "resolution: $(vs.nu) × $(vs.nv)\n",
              "ghost in u: $(1-firstindex(vs.u[:, 1]))\n",
              "ghost in v: $(1-firstindex(vs.v[1, :]))\n")
end

function Base.show(io::IO, vs::VSpace3D{TR,TI,TA}) where {TR,TI,TA}
    print(io, "VelocitySpace3D{$TR,$TI,$TA}\n",
              "domain: ($(vs.u0),$(vs.u1)) × ($(vs.v0),$(vs.v1)) × ($(vs.w0),$(vs.w1))\n",
              "resolution: $(vs.nu) × $(vs.nv) × $(vs.nw)\n",
              "ghost in u: $(1-firstindex(vs.u[:, 1, 1]))\n",
              "ghost in v: $(1-firstindex(vs.v[1, :, 1]))\n",
              "ghost in w: $(1-firstindex(vs.w[1, 1, :]))\n")
end
