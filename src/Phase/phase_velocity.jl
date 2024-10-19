# ============================================================
# Methods of Particle Velocity Space
# ============================================================

"""
$(TYPEDEF)

1D velocity space

## Fields

$(FIELDS)
"""
struct VSpace1D{TR,TI,TA,TB} <: AbstractVelocitySpace1D
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
    NU::TI;
    type="rectangle",
    ng=zero(NU)::TI,
    precision=Float64,
) where {TI<:Integer}
    δ = (U1 - U0) / NU
    u = begin
        if ng > 0
            OffsetArray{precision}(undef, 1-ng:NU+ng)
        else
            Array{precision}(undef, NU)
        end
    end
    du = similar(u)
    weights = similar(u)

    if type == "rectangle" # rectangular
        for i in eachindex(u)
            u[i] = U0 + (i - 0.5) * δ
            du[i] = δ
            weights[i] = δ
        end
    elseif type == "newton" # newton-cotes
        for i in eachindex(u)
            u[i] = U0 + (i - 0.5) * δ
            du[i] = δ
            weights[i] = newton_cotes(i + ng, NU + ng * 2) * δ
        end
    elseif type == "algebra" # algebraic
        _nu = NU + 1
        _u = [U1 / (_nu - 1)^3 * (-_nu + 1 + 2 * (i - 1))^3 for i in 1:_nu]
        u .= (_u[1:end-1] .+ _u[2:end]) ./ 2
        du .= _u[2:end] - _u[1:end-1]
        weights .= du
    else
        throw("No velocity quadrature available")
    end

    return VSpace1D{precision,TI,typeof(u),typeof(weights)}(U0, U1, NU, u, du, weights)
end

VSpace1D() = VSpace1D(-5, 5, 50)
VSpace1D(U0::T, U1::T) where {T<:Real} = VSpace1D(U0, U1, 50)

# ------------------------------------------------------------
# Multi-component
# ------------------------------------------------------------

function MVSpace1D(
    Ui0,
    Ui1,
    Ue0,
    Ue1,
    NU::TI;
    type="rectangle",
    ng=zero(NU)::TI,
    precision=Float64,
) where {TI<:Integer}
    u0 = precision.([Ui0, Ue0])
    u1 = precision.([Ui1, Ue1])
    δ = (u1 .- u0) ./ NU
    u = begin
        if ng > 0
            OffsetArray{precision}(undef, 1-ng:NU+ng, 1:2)
        else
            Array{precision}(undef, NU, 2)
        end
    end
    du = similar(u)
    weights = similar(u)

    if type == "rectangle" # rectangular
        for j in axes(u, 2), i in axes(u, 1)
            u[i, j] = u0[j] + (i - 0.5) * δ[j]
            du[i, j] = δ[j]
            weights[i, j] = δ[j]
        end
    elseif type == "newton" # newton-cotes
        for j in axes(u, 2), i in axes(u, 1)
            u[i, j] = u0[j] + (i - 0.5) * δ[j]
            du[i, j] = δ[j]
            weights[i, j] = newton_cotes(i + ng, NU + ng * 2) * δ[j]
        end
    else
        throw("No velocity quadrature available")
    end

    return VSpace1D{typeof(u0),TI,typeof(u),typeof(weights)}(u0, u1, NU, u, du, weights)
end

MVSpace1D() = MVSpace1D(-5, 5, -10, 10, 28)
MVSpace1D(U0::T, U1::T, V0::T, V1::T) where {T<:Real} = MVSpace1D(U0, U1, V0, V1, 28)

"""
$(TYPEDEF)

2D velocity space

## Fields

$(FIELDS)
"""
struct VSpace2D{TR,TI,TA} <: AbstractVelocitySpace2D
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
    NV::TI;
    type="rectangle",
    ngu=zero(NU)::TI,
    ngv=zero(NV)::TI,
    precision=Float64,
) where {TI<:Integer}
    δu = (U1 - U0) / NU
    δv = (V1 - V0) / NV
    u = begin
        if ngu > 0 || ngv > 0
            OffsetArray{precision}(undef, 1-ngu:NU+ngu, 1-ngv:NV+ngv)
        else
            Array{precision}(undef, NU, NV)
        end
    end
    v = similar(u)
    du = similar(u)
    dv = similar(u)
    weights = similar(u)

    if type == "rectangle" # rectangular
        for j in axes(u, 2)
            for i in axes(u, 1)
                u[i, j] = U0 + (i - 0.5) * δu
                v[i, j] = V0 + (j - 0.5) * δv
                du[i, j] = δu
                dv[i, j] = δv
                weights[i, j] = δu * δv
            end
        end
    elseif type == "newton" # newton-cotes
        for j in axes(u, 2)
            for i in axes(u, 1)
                u[i, j] = U0 + (i - 0.5) * δu
                v[i, j] = V0 + (j - 0.5) * δv
                du[i, j] = δu
                dv[i, j] = δv
                weights[i, j] =
                    newton_cotes(i + ngu, NU + ngu * 2) *
                    δu *
                    newton_cotes(j + ngv, NV + ngv * 2) *
                    δv
            end
        end
    elseif type == "algebra"
        _nu = NU + 1
        _nv = NV + 1
        _u = [U1 / (_nu - 1)^3 * (-_nu + 1 + 2 * (i - 1))^3 for i in 1:_nu]
        _v = [V1 / (_nv - 1)^3 * (-_nv + 1 + 2 * (j - 1))^3 for j in 1:_nv]
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
    elseif type == "maxwell"
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

    return VSpace2D{precision,TI,typeof(u)}(U0, U1, NU, V0, V1, NV, u, v, du, dv, weights)
end

VSpace2D() = VSpace2D(-5, 5, 28, -5, 5, 28)
VSpace2D(U0::T, U1::T, V0::T, V1::T) where {T<:Real} = VSpace2D(U0, U1, 28, V0, V1, 28)

# ------------------------------------------------------------
# Multi-component
# ------------------------------------------------------------

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
    NV::TI;
    type="rectangle",
    ngu=zero(NU)::TI,
    ngv=zero(NV)::TI,
    precision=Float64,
) where {TI<:Integer}
    u0 = precision.([Ui0, Ue0])
    u1 = precision.([Ui1, Ue1])
    δu = (u1 .- u0) ./ NU
    v0 = precision.([Vi0, Ve0])
    v1 = precision.([Vi1, Ve1])
    δv = (v1 .- v0) ./ NV
    u = begin
        if ngu > 0 || ngv > 0
            OffsetArray{precision}(undef, 1-ngu:NU+ngu, 1-ngv:NV+ngv, 1:2)
        else
            Array{precision}(undef, NU, NV, 2)
        end
    end
    v = similar(u)
    du = similar(u)
    dv = similar(u)
    weights = similar(u)

    if type == "rectangle" # rectangular
        for k in axes(u, 3), j in axes(u, 2), i in axes(u, 1)
            u[i, j, k] = u0[k] + (i - 0.5) * δu[k]
            v[i, j, k] = v0[k] + (j - 0.5) * δv[k]
            du[i, j, k] = δu[k]
            dv[i, j, k] = δv[k]
            weights[i, j, k] = δu[k] * δv[k]
        end
    elseif type == "newton" # newton-cotes
        for k in axes(u, 3), j in axes(u, 2), i in axes(u, 1)
            u[i, j, k] = u0[k] + (i - 0.5) * δu[k]
            v[i, j, k] = v0[k] + (j - 0.5) * δv[k]
            du[i, j, k] = δu[k]
            dv[i, j, k] = δv[k]
            weights[i, j, k] =
                newton_cotes(i + ngu, NU + ngu * 2) *
                δu[k] *
                newton_cotes(j + ngv, NV + ngv * 2) *
                δv[k]
        end
    else
        throw("No velocity quadrature available")
    end

    return VSpace2D{typeof(u0),TI,typeof(u)}(u0, u1, NU, v0, v1, NV, u, v, du, dv, weights)
end

MVSpace2D() = MVSpace2D(-5, 5, -10, 10, 28, -5, 5, -10, 10, 28)
MVSpace2D(U0::T, U1::T, V0::T, V1::T) where {T<:Real} =
    MVSpace2D(U0, U1, U0, U1, 28, V0, V1, V0, V1, 28)

"""
$(TYPEDEF)

3D velocity space

## Fields

$(FIELDS)
"""
struct VSpace3D{TR,TI,TA} <: AbstractVelocitySpace3D
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
    NW::TI;
    type="rectangle",
    ngu=zero(NU)::TI,
    ngv=zero(NV)::TI,
    ngw=zero(NW)::TI,
    precision=Float64,
) where {TI<:Integer}
    δu = (U1 - U0) / NU
    δv = (V1 - V0) / NV
    δw = (W1 - W0) / NW
    u = begin
        if ngu > 0 || ngv > 0 || ngw > 0
            OffsetArray{precision}(undef, 1-ngu:NU+ngu, 1-ngv:NV+ngv, 1-ngw:NW+ngw)
        else
            Array{precision}(undef, NU, NV, NW)
        end
    end
    v = similar(u)
    w = similar(u)
    du = similar(u)
    dv = similar(u)
    dw = similar(u)
    weights = similar(u)

    if type == "rectangle" # rectangular
        for k in axes(u, 3), j in axes(u, 2), i in axes(u, 1)
            u[i, j, k] = U0 + (i - 0.5) * δu
            v[i, j, k] = V0 + (j - 0.5) * δv
            w[i, j, k] = W0 + (k - 0.5) * δw
            du[i, j, k] = δu
            dv[i, j, k] = δv
            dw[i, j, k] = δw
            weights[i, j, k] = δu * δv * δw
        end
    elseif type == "newton" # newton-cotes
        for k in axes(u, 3), j in axes(u, 2), i in axes(u, 1)
            u[i, j, k] = U0 + (i - 0.5) * δu
            v[i, j, k] = V0 + (j - 0.5) * δv
            w[i, j, k] = W0 + (k - 0.5) * δw
            du[i, j, k] = δu
            dv[i, j, k] = δv
            dw[i, j, k] = δw
            weights[i, j, k] =
                newton_cotes(i + ngu, NU + ngu * 2) *
                δu *
                newton_cotes(j + ngv, NV + ngv * 2) *
                δv *
                newton_cotes(k + ngw, NW + ngw * 2) *
                δw
        end
    elseif type == "algebra"
        _nu = NU + 1
        _nv = NV + 1
        _nw = NW + 1
        _u = [U1 / (_nu - 1)^3 * (-_nu + 1 + 2 * (i - 1))^3 for i in 1:_nu]
        _v = [V1 / (_nv - 1)^3 * (-_nv + 1 + 2 * (j - 1))^3 for j in 1:_nv]
        _w = [W1 / (_nw - 1)^3 * (-_nw + 1 + 2 * (k - 1))^3 for k in 1:_nw]
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

    return VSpace3D{precision,TI,typeof(u)}(
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

# ------------------------------------------------------------
# Multi-component
# ------------------------------------------------------------

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
    NW::TI;
    type="rectangle",
    ngu=zero(NU)::TI,
    ngv=zero(NV)::TI,
    ngw=zero(NW)::TI,
    precision=Float64,
) where {TI<:Integer}
    u0 = precision.([Ui0, Ue0])
    u1 = precision.([Ui1, Ue1])
    δu = (u1 .- u0) ./ NU
    v0 = precision.([Vi0, Ve0])
    v1 = precision.([Vi1, Ve1])
    δv = (v1 .- v0) ./ NV
    w0 = precision.([Wi0, We0])
    w1 = precision.([Wi1, We1])
    δw = (w1 .- w0) ./ NW

    u = begin
        if ngu > 0 || ngv > 0 || ngw > 0
            OffsetArray{precision}(undef, 1-ngu:NU+ngu, 1-ngv:NV+ngv, 1-ngw:NW+ngw, 1:2)
        else
            Array{precision}(undef, NU, NV, NW, 2)
        end
    end
    v = similar(u)
    w = similar(u)
    du = similar(u)
    dv = similar(u)
    dw = similar(u)
    weights = similar(u)

    if type == "rectangle" # rectangular
        for l in axes(u, 4), k in axes(u, 3), j in axes(u, 2), i in axes(u, 1)
            u[i, j, k, l] = u0[l] + (i - 0.5) * δu[l]
            v[i, j, k, l] = v0[l] + (j - 0.5) * δv[l]
            w[i, j, k, l] = w0[l] + (k - 0.5) * δw[l]
            du[i, j, k, l] = δu[l]
            dv[i, j, k, l] = δv[l]
            dw[i, j, k, l] = δw[l]
            weights[i, j, k, l] = δu[l] * δv[l] * δw[l]
        end
    elseif type == "newton" # newton-cotes
        for l in axes(u, 4), k in axes(u, 3), j in axes(u, 2), i in axes(u, 1)
            u[i, j, k, l] = u0[l] + (i - 0.5) * δu[l]
            v[i, j, k, l] = v0[l] + (j - 0.5) * δv[l]
            w[i, j, k, l] = w0[l] + (k - 0.5) * δw[l]
            du[i, j, k, l] = δu[l]
            dv[i, j, k, l] = δv[l]
            dw[i, j, k, l] = δw[l]
            weights[i, j, k, l] =
                newton_cotes(i + ngu, NU + ngu * 2) *
                δu[l] *
                newton_cotes(j + ngv, NV + ngv * 2) *
                δv[l] *
                newton_cotes(k + ngw, NW + ngw * 2) *
                δw[l]
        end
    else
        throw("No velocity quadrature available")
    end

    return VSpace3D{typeof(u0),TI,typeof(u)}(
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
$(SIGNATURES)

Generate VelocitySpace
"""
function set_velocity(;
    space,
    nSpecies,
    umin=nothing,
    umax=nothing,
    nu=nothing,
    vMeshType=nothing,
    nug=nothing,
    mi=nothing,
    me=nothing,
    vmin=nothing,
    vmax=nothing,
    nv=nothing,
    nvg=nothing,
    wmin=nothing,
    wmax=nothing,
    nw=nothing,
    nwg=nothing,
    kwargs...,
)
    Dv = parse(Int, space[5])
    if Dv == 0
        vSpace = nothing
    elseif Dv == 1
        if nSpecies == 1
            vSpace = VSpace1D(umin, umax, nu; type=vMeshType, ng=nug)
        elseif nSpecies == 2
            ue0 = umin * sqrt(mi / me)
            ue1 = umax * sqrt(mi / me)
            vSpace = MVSpace1D(umin, umax, ue0, ue1, nu; type=vMeshType, ng=nug)
        else
            throw("The velocity space only supports up to two species.")
        end
    elseif Dv == 2
        if nSpecies == 1
            vSpace =
                VSpace2D(umin, umax, nu, vmin, vmax, nv; type=vMeshType, ngu=nug, ngv=nvg)
        elseif nSpecies == 2
            ue0 = umin * sqrt(mi / me)
            ue1 = umax * sqrt(mi / me)
            ve0 = vmin * sqrt(mi / me)
            ve1 = vmax * sqrt(mi / me)
            vSpace = MVSpace2D(
                umin,
                umax,
                ue0,
                ue1,
                nu,
                vmin,
                vmax,
                ve0,
                ve1,
                nv;
                type=vMeshType,
                ngu=nug,
                ngv=nvg,
            )
        else
            throw("The velocity space only supports up to two species.")
        end
    elseif Dv == 3
        if nSpecies == 1
            vSpace = VSpace3D(
                umin,
                umax,
                nu,
                vmin,
                vmax,
                nv,
                wmin,
                wmax,
                nw;
                type=vMeshType,
                ngu=nug,
                ngv=nvg,
                ngw=nwg,
            )
        elseif nSpecies == 2
            ue0 = umin * sqrt(mi / me)
            ue1 = umax * sqrt(mi / me)
            ve0 = vmin * sqrt(mi / me)
            ve1 = vmax * sqrt(mi / me)
            we0 = wmin * sqrt(mi / me)
            we1 = wmax * sqrt(mi / me)
            vSpace = MVSpace3D(
                umin,
                umax,
                ue0,
                ue1,
                nu,
                vmin,
                vmax,
                ve0,
                ve1,
                nv,
                wmin,
                wmax,
                we0,
                we1,
                nw;
                type=vMeshType,
                ngu=nug,
                ngv=nvg,
                ngw=nwg,
            )
        else
            throw("The velocity space only supports up to two species.")
        end
    end

    return vSpace
end

"""
$(SIGNATURES)
"""
set_velocity(dict::Union{AbstractDict,NamedTuple}) = set_velocity(; dict...)

"""
$(TYPEDEF)

Unstructured velocity space

## Fields

$(FIELDS)
"""
struct UnstructVSpace{TR<:Union{Real,AV},TI<:Union{Integer,AV},TA<:AA,TB<:AV{<:Real}} <:
       AbstractVelocitySpace
    u0::TR
    u1::TR
    nu::TI
    u::TA
    weights::TB
end

# ------------------------------------------------------------
# Extended Base.show()
# ------------------------------------------------------------

function Base.show(io::IO, vs::VSpace1D{TR,TI,TA,TB}) where {TR,TI,TA,TB}
    return print(
        io,
        "VelocitySpace1D{$TR,$TI,$TA,$TB}\n",
        "domain: ($(vs.u0),$(vs.u1))\n",
        "resolution: $(vs.nu)\n",
    )
end

function Base.show(io::IO, vs::VSpace2D{TR,TI,TA}) where {TR,TI,TA}
    return print(
        io,
        "VelocitySpace2D{$TR,$TI,$TA}\n",
        "domain: ($(vs.u0),$(vs.u1)) × ($(vs.v0),$(vs.v1))\n",
        "resolution: $(vs.nu) × $(vs.nv)\n",
    )
end

function Base.show(io::IO, vs::VSpace3D{TR,TI,TA}) where {TR,TI,TA}
    return print(
        io,
        "VelocitySpace3D{$TR,$TI,$TA}\n",
        "domain: ($(vs.u0),$(vs.u1)) × ($(vs.v0),$(vs.v1)) × ($(vs.w0),$(vs.w1))\n",
        "resolution: $(vs.nu) × $(vs.nv) × $(vs.nw)\n",
    )
end
