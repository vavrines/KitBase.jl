# ============================================================
# Methods of Particle Velocity Space
# ============================================================

"""
    struct VSpace1D
        u0::TR
        u1::TR
        nu::TI
        u::TA
        du::TA
        weights::TB
    end

1D velocity space

"""
struct VSpace1D{TR<:Real,TI<:Integer,TA<:AA{<:Real},TB<:AV{<:Real}} <:
       AbstractVelocitySpace1D
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
    u = begin
        if NG > 0
            OffsetArray{PRECISION}(undef, 1-NG:NU+NG)
        else
            Array{PRECISION}(undef, NU)
        end
    end
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
    elseif TYPE == "algebra" # algebraic
        _nu = NU + 1
        _u = [U1 / (_nu - 1)^3 * (-_nu + 1 + 2 * (i - 1))^3 for i = 1:_nu]
        u .= (_u[1:end-1] .+ _u[2:end]) ./ 2
        du .= _u[2:end] - _u[1:end-1]
        weights .= du
    else
        throw("No velocity quadrature available")
    end

    return VSpace1D{PRECISION,TI,typeof(u),typeof(weights)}(U0, U1, NU, u, du, weights)

end

VSpace1D() = VSpace1D(-5, 5, 50)
VSpace1D(U0::T, U1::T) where {T<:Real} = VSpace1D(U0, U1, 50)


"""
    struct VSpace2D
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
struct VSpace2D{TR<:Real,TI<:Integer,TA<:AA{<:Real}} <: AbstractVelocitySpace2D
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
    u = begin
        if NGU > 0 || NGV > 0
            OffsetArray{PRECISION}(undef, 1-NGU:NU+NGU, 1-NGV:NV+NGV)
        else
            Array{PRECISION}(undef, NU, NV)
        end
    end
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
    elseif TYPE == "algebra"
        _nu = NU + 1
        _nv = NV + 1
        _u = [U1 / (_nu - 1)^3 * (-_nu + 1 + 2 * (i - 1))^3 for i = 1:_nu]
        _v = [V1 / (_nv - 1)^3 * (-_nv + 1 + 2 * (j - 1))^3 for j = 1:_nv]
        __u = (_u[1:end-1] .+ _u[2:end]) ./ 2
        __v = (_v[1:end-1] .+ _v[2:end]) ./ 2
        u1, v1 = meshgrid(__u, __v)
        u .= u1 |> permutedims
        v .= v1 |> permutedims

        #_wx = [ 3.0 * (-_nu/2 + i - 0.5)^2 / (_nu/2 - 0.5)^3 for i = 1:_nu]
        #_wy = [ 3.0 * (-_nv/2 + j - 0.5)^2 / (_nv/2 - 0.5)^3 for j = 1:_nv]
        #wx = (_wx[1:end-1] .+ _wx[2:end]) ./ 2
        #wy = (_wy[1:end-1] .+ _wy[2:end]) ./ 2

        _du = _u[2:end] - _u[1:end-1]
        _dv = _v[2:end] - _v[1:end-1]
        du1, dv1 = meshgrid(_du, _dv)
        du .= du1 |> permutedims
        dv .= dv1 |> permutedims

        for j in axes(u, 2)
            for i in axes(u, 1)
                weights[i, j] = du[i, j] * dv[i, j]
                #weights[i, j] = wx[i] * wy[j]
            end
        end
    elseif TYPE == "maxwell"
        _u, _uw = maxwell_quadrature(NU, U1 / 5.76)
        _v, _vw = maxwell_quadrature(NV, V1 / 5.76)
        u1, v1 = meshgrid(_u, _v)
        u .= u1 |> permutedims
        v .= v1 |> permutedims
        for j in axes(weights, 2), i in axes(weights, 1)
            weights[i, j] = _uw[i] * _vw[j]
        end

        _du = zero(_u)
        for i in eachindex(_u)
            if i == 1
                _du[i] = _u[2] - _u[1]
            elseif i == length(_u)
                _du[i] = _u[end] - _u[end-1]
            else
                _du[i] = (_u[i+1] - _u[i-1]) / 2
            end
        end
        _dv = zero(_v)
        for i in eachindex(_v)
            if i == 1
                _dv[i] = _v[2] - _v[1]
            elseif i == length(_v)
                _dv[i] = _v[end] - _v[end-1]
            else
                _dv[i] = (_v[i+1] - _v[i-1]) / 2
            end
        end
        du1, dv1 = meshgrid(_du, _dv)
        du .= du1 |> permutedims
        dv .= dv1 |> permutedims
    else
        throw("No velocity quadrature available")
    end

    return VSpace2D{PRECISION,TI,typeof(u)}(U0, U1, NU, V0, V1, NV, u, v, du, dv, weights)

end

VSpace2D() = VSpace2D(-5, 5, 28, -5, 5, 28)
VSpace2D(U0::T, U1::T, V0::T, V1::T) where {T<:Real} = VSpace2D(U0, U1, 28, V0, V1, 28)


"""
    struct VSpace3D
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
struct VSpace3D{TR<:Real,TI<:Integer,TA<:AA{<:Real}} <: AbstractVelocitySpace3D
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
    u = begin
        if NGU > 0 || NGV > 0 || NGW > 0
            OffsetArray{PRECISION}(undef, 1-NGU:NU+NGU, 1-NGV:NV+NGV, 1-NGW:NW+NGW)
        else
            Array{PRECISION}(undef, NU, NV, NW)
        end
    end
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
    elseif TYPE == "algebra"
        _nu = NU + 1
        _nv = NV + 1
        _nw = NW + 1
        _u = [U1 / (_nu - 1)^3 * (-_nu + 1 + 2 * (i - 1))^3 for i = 1:_nu]
        _v = [V1 / (_nv - 1)^3 * (-_nv + 1 + 2 * (j - 1))^3 for j = 1:_nv]
        _w = [W1 / (_nw - 1)^3 * (-_nw + 1 + 2 * (k - 1))^3 for k = 1:_nw]
        __u = (_u[1:end-1] .+ _u[2:end]) ./ 2
        __v = (_v[1:end-1] .+ _v[2:end]) ./ 2
        __w = (_w[1:end-1] .+ _w[2:end]) ./ 2
        u1, v1, w1 = meshgrid(__u, __v, __w)
        u .= permutedims(u1, [3, 2, 1])
        v .= permutedims(v1, [3, 2, 1])
        w .= permutedims(w1, [3, 2, 1])

        _du = _u[2:end] - _u[1:end-1]
        _dv = _v[2:end] - _v[1:end-1]
        _dw = _w[2:end] - _w[1:end-1]
        du1, dv1, dw1 = meshgrid(_du, _dv, _dw)
        du .= permutedims(du1, [3, 2, 1])
        dv .= permutedims(dv1, [3, 2, 1])
        dw .= permutedims(dw1, [3, 2, 1])

        for k in axes(u, 3), j in axes(u, 2), i in axes(u, 1)
            weights[i, j, k] = du[i, j, k] * dv[i, j, k] * dw[i, j, k]
        end
    else
        throw("No velocity quadrature available")
    end

    return VSpace3D{PRECISION,TI,typeof(u)}(
        U0,
        U1,
        NU,
        V0,
        V1,
        NV,
        W0,
        W1,
        NW,
        u,
        v,
        w,
        du,
        dv,
        dw,
        weights,
    )

end

VSpace3D() = VSpace3D(-5, 5, 28, -5, 5, 28, -5, 5, 28)
VSpace3D(U0::T, U1::T, V0::T, V1::T, W0::T, W1::T) where {T<:Real} =
    VSpace3D(U0, U1, 28, V0, V1, 28, W0, W1, 28)


"""
    struct MVSpace1D
        u0::TR
        u1::TR
        nu::TI
        u::TA
        du::TA
        weights::TA
    end

1D multi-component velocity space

"""
struct MVSpace1D{TR<:AA{<:Real,1},TI<:Integer,TA<:AA{<:Real,2}} <: AbstractVelocitySpace1D
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
    u = begin
        if NG > 0
            OffsetArray{PRECISION}(undef, 1-NG:NU+NG, 1:2)
        else
            Array{PRECISION}(undef, NU, 2)
        end
    end
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
    struct MVSpace2D
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
struct MVSpace2D{TR<:AA{<:Real,1},TI<:Integer,TA<:AA{Float64,3}} <: AbstractVelocitySpace2D
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
    u = begin
        if NGU > 0 || NGV > 0
            OffsetArray{PRECISION}(undef, 1-NGU:NU+NGU, 1-NGV:NV+NGV, 1:2)
        else
            Array{PRECISION}(undef, NU, NV, 2)
        end
    end
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
    struct MVSpace3D
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
struct MVSpace3D{TR<:AA{<:Real,1},TI<:Integer,TA<:AA{<:Real,4}} <: AbstractVelocitySpace3D
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

    u = begin
        if NGU > 0 || NGV > 0 || NGW > 0
            OffsetArray{PRECISION}(undef, 1-NGU:NU+NGU, 1-NGV:NV+NGV, 1-NGW:NW+NGW, 1:2)
        else
            Array{PRECISION}(undef, NU, NV, NW, 2)
        end
    end
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

    return MVSpace3D{typeof(u0),TI,typeof(u)}(
        u0,
        u1,
        NU,
        v0,
        v1,
        NV,
        w0,
        w1,
        NW,
        u,
        v,
        v,
        du,
        dv,
        dw,
        weights,
    )

end

MVSpace3D() = MVSpace3D(-5, 5, -10, 10, 20, -5, 5, -10, 10, 20, -5, 5, -10, 10, 20)
MVSpace3D(U0::T, U1::T, V0::T, V1::T, W0::T, W1::T) where {T<:Real} =
    MVSpace3D(U0, U1, U0, U1, 20, V0, V1, V0, V1, 20, W0, W1, W0, W1, 20)


"""
    struct UnstructVSpace
        u0::TR
        u1::TR
        nu::TI
        u::TA
        weights::TB
    end

Unstructured velocity space

"""
struct UnstructVSpace{TR<:Union{Real,AV},TI<:Union{Integer,AV},TA<:AA,TB<:AV{<:Real}} <:
       AbstractVelocitySpace
    u0::TR
    u1::TR
    nu::TI
    u::TA
    weights::TB
end


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
    print(
        io,
        "VelocitySpace1D{$TR,$TI,$TA,$TB}\n",
        "domain: ($(vs.u0),$(vs.u1))\n",
        "resolution: $(vs.nu)\n",
        "ghost: $(1-firstindex(vs.u))\n",
    )
end

function Base.show(io::IO, vs::VSpace2D{TR,TI,TA}) where {TR,TI,TA}
    print(
        io,
        "VelocitySpace2D{$TR,$TI,$TA}\n",
        "domain: ($(vs.u0),$(vs.u1)) × ($(vs.v0),$(vs.v1))\n",
        "resolution: $(vs.nu) × $(vs.nv)\n",
        "ghost in u: $(1-firstindex(vs.u[:, 1]))\n",
        "ghost in v: $(1-firstindex(vs.v[1, :]))\n",
    )
end

function Base.show(io::IO, vs::VSpace3D{TR,TI,TA}) where {TR,TI,TA}
    print(
        io,
        "VelocitySpace3D{$TR,$TI,$TA}\n",
        "domain: ($(vs.u0),$(vs.u1)) × ($(vs.v0),$(vs.v1)) × ($(vs.w0),$(vs.w1))\n",
        "resolution: $(vs.nu) × $(vs.nv) × $(vs.nw)\n",
        "ghost in u: $(1-firstindex(vs.u[:, 1, 1]))\n",
        "ghost in v: $(1-firstindex(vs.v[1, :, 1]))\n",
        "ghost in w: $(1-firstindex(vs.w[1, 1, :]))\n",
    )
end


function maxwell_quadrature(N::Integer, C = 1::Real)
    @assert N <= 33

    py"""
    import numpy as np
    from numpy import linalg as LA

    def dvGH(N2,C):
        N = N2//2

        a = np.zeros(N)
        b = np.zeros(N)
        a[0] = 1.0/np.sqrt(np.pi)
        a[1] = 2.0/np.sqrt(np.pi)/(np.pi-2.0)
        b[1] = a[0]/( a[0] + a[1])/2.0

        for i in range(2,N):
            b[i] = (i-1)+1.0/2.0-b[i-1]-a[i-1]**2
            a[i] = (i**2/4.0/b[i]-b[i-1]-1.0/2)/a[i-1]-a[i-1]

        J = np.diag(a) + np.diag(np.sqrt(b[1:N]),1) \
        + np.diag(np.sqrt(b[1:N]),-1)

        v,V = LA.eig(J)

        w = V[0,:]*V[0,:]*np.sqrt(np.pi)/2.0

        vw = np.transpose(np.vstack((v,w)))
        vw = vw[vw[:,0].argsort()]
        v = vw[:,0]
        w = vw[:,1]

        Xis = np.hstack((-np.flipud(v),v))
        weights = np.hstack((np.flipud(w),w))
        weights = weights*np.exp(Xis**2)*C
        Xis = Xis*C
        return (Xis, weights)
    """

    p, w = py"dvGH"(N, C)
end
