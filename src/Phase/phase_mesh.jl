"""
$(SIGNATURES)

Generate 1D mesh with quadrature rule

## Arguments
- `U0`: lower bound
- `U1`: upper bound
- `NU`: number of quadrature points
- `type`: quadrature type, default is "rectangle"
- `precision`: precision of the quadrature rule, default is Float64
"""
function mesh_quadrature(U0, U1, NU::Integer; type = "rectangle", precision = Float64)
    δ = (U1 - U0) / NU
    u = Array{precision}(undef, NU)
    du = similar(u)
    weights = similar(u)

    if type == "rectangle" # rectangular
        u .= uniform_mesh(U0, NU, δ)
        du .= δ
        weights .= δ
    elseif type == "newton" # newton-cotes
        @assert isodd(NU) "The number of quadrature points for the Newton-Cotes rule must be odd."
        u .= uniform_mesh(U0, NU, δ)
        du .= δ

        for i in eachindex(u)
            weights[i] = newton_cotes(i, NU) * δ
        end
    elseif type == "algebra" # algebraic
        _nu = NU + 1
        _u = [U1 / (_nu - 1)^3 * (-_nu + 1 + 2 * (i - 1))^3 for i = 1:_nu]
        u .= (_u[1:end-1] .+ _u[2:end]) ./ 2
        du .= _u[2:end] - _u[1:end-1]
        weights .= du
    else
        throw("No velocity quadrature available")
    end

    return u, du, weights
end

"""
$(SIGNATURES)

Generate 2D mesh with quadrature rule

## Arguments
- `U0`,`V0`: lower bound
- `U1`,`V1`: upper bound
- `NU`,`NV`: number of quadrature points
- `type`: quadrature type, default is "rectangle"
- `precision`: precision of the quadrature rule, default is Float64
"""
function mesh_quadrature(
    U0,
    U1,
    NU::TI,
    V0,
    V1,
    NV::TI;
    type = "rectangle",
    precision = Float64,
) where {TI<:Integer}

    δu = (U1 - U0) / NU
    δv = (V1 - V0) / NV
    u = Array{precision}(undef, NU, NV)
    v = similar(u)
    du = similar(u)
    dv = similar(u)
    weights = similar(u)

    if type == "rectangle" # rectangular
        _u = uniform_mesh(U0, NU, δu)
        _v = uniform_mesh(V0, NV, δv)
        u, v = ndgrid(_u, _v)
        du .= δu
        dv .= δv
        weights .= δu * δv
    elseif type == "newton" # newton-cotes
        @assert isodd(NU) && isodd(NV) "The number of quadrature points for the Newton-Cotes rule must be odd."
        _u = uniform_mesh(U0, NU, δu)
        _v = uniform_mesh(V0, NV, δv)
        u, v = ndgrid(_u, _v)
        du .= δu
        dv .= δv
        for j in axes(u, 2)
            for i in axes(u, 1)
                weights[i, j] = newton_cotes(i, NU) * δu * newton_cotes(j, NV) * δv
            end
        end
    elseif type == "algebra"
        _nu = NU + 1
        _nv = NV + 1
        _u = [U1 / (_nu - 1)^3 * (-_nu + 1 + 2 * (i - 1))^3 for i = 1:_nu]
        _v = [V1 / (_nv - 1)^3 * (-_nv + 1 + 2 * (j - 1))^3 for j = 1:_nv]
        __u = (_u[1:end-1] .+ _u[2:end]) ./ 2
        __v = (_v[1:end-1] .+ _v[2:end]) ./ 2
        u1, v1 = meshgrid(__u, __v)
        u .= u1 |> permutedims
        v .= v1 |> permutedims

        _du = _u[2:end] - _u[1:end-1]
        _dv = _v[2:end] - _v[1:end-1]
        du1, dv1 = meshgrid(_du, _dv)
        du .= du1 |> permutedims
        dv .= dv1 |> permutedims

        for j in axes(u, 2)
            for i in axes(u, 1)
                weights[i, j] = du[i, j] * dv[i, j]
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

    return u, v, du, dv, weights

end

"""
$(SIGNATURES)

Generate 3D mesh with quadrature rule

## Arguments
- `U0`,`V0`,`W0`: lower bound
- `U1`,`V1`,`W1`: upper bound
- `NU`,`NV`,`NW`: number of quadrature points
- `type`: quadrature type, default is "rectangle"
- `precision`: precision of the quadrature rule, default is Float64
"""
function mesh_quadrature(
    U0,
    U1,
    NU::TI,
    V0,
    V1,
    NV::TI,
    W0,
    W1,
    NW::TI;
    type = "rectangle",
    precision = Float64,
) where {TI<:Integer}

    δu = (U1 - U0) / NU
    δv = (V1 - V0) / NV
    δw = (W1 - W0) / NW
    u = Array{precision}(undef, NU, NV, NW)
    v = similar(u)
    w = similar(u)
    du = similar(u)
    dv = similar(u)
    dw = similar(u)
    weights = similar(u)

    if type == "rectangle" # rectangular
        _u = uniform_mesh(U0, NU, δu)
        _v = uniform_mesh(V0, NV, δv)
        _w = uniform_mesh(W0, NW, δw)
        u, v, w = ndgrid(_u, _v, _w)
        du .= δu
        dv .= δv
        dw .= δw
        weights .= δu * δv * δw
    elseif type == "newton" # newton-cotes
        @assert isodd(NU) && isodd(NV) && isodd(NW) "The number of quadrature points for the Newton-Cotes rule must be odd."
        _u = uniform_mesh(U0, NU, δu)
        _v = uniform_mesh(V0, NV, δv)
        _w = uniform_mesh(W0, NW, δw)
        u, v, w = ndgrid(_u, _v, _w)
        du .= δu
        dv .= δv
        dw .= δw
        for k in axes(u, 3), j in axes(u, 2), i in axes(u, 1)
            weights[i, j, k] =
                newton_cotes(i, NU) *
                δu *
                newton_cotes(j, NV) *
                δv *
                newton_cotes(k, NW) *
                δw
        end
    elseif type == "algebra"
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

    return u, v, w, du, dv, dw, weights

end
