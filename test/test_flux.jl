begin
    inK = 2.0
    γ = 5 / 3
    wL = [1.13, 0.13, 0.8]
    wR = [0.3, -0.2, 1.5]
    primL = KB.conserve_prim(wL, γ)
    primR = KB.conserve_prim(wR, γ)

    u = collect(-5.0:0.2:5.0)
    ω = ones(length(u)) / length(u)

    hL = KB.maxwellian(u, primL)
    bL = hL .* inK ./ (2.0 * primL[end])
    hR = KB.maxwellian(u, primR)
    bR = hR .* inK ./ (2.0 * primR[end])

    dt = 1e-3
    dx = 1e-2

    fw = similar(wL)
    fh = similar(hL)
    fb = similar(bL)
end

#--- fluid flux ---#
KB.flux_upwind(rand(), rand(), rand(2), rand(2), dt)
KB.flux_lax!(fw, wL, wR, γ, dt, dx)
KB.flux_godunov!(fw, wL, wR, γ, dt)
KB.flux_hll!(fw, wL, wR, γ, dt)
KB.flux_hllc!(fw, wL, wR, γ, dt)
fw1 = zeros(4)
KB.flux_hllc!(fw1, [wL[1:2]; 0.0; wL[end]], [wR[1:2]; 0.0; wR[end]], γ, dt)
KB.flux_roe!(fw, wL, wR, γ, dt)
KB.flux_roe!(zeros(4), [1.0, 0.3, 0.0, 1.0], [0.3, -0.1, 0.0, 2.0], γ, dt)

#--- gks flux ---#
KB.flux_gks(0.3, 1e-3, dt)
KB.flux_gks(0.3, 1e-3, dt, 1e-1, 1.0)
KB.flux_gks(1.0, 0.125, 1e-3, dt, 1e-2, 1e-2)

KB.flux_gks!(zeros(3), [1.0, 0.0, 2.0], inK, γ, 1e-3, 0.72)
KB.flux_gks!(zeros(4), [1.0, 0.0, 0.0, 2.0], inK, γ, 1e-3, 0.72)
KB.flux_gks!(zeros(5), [1.0, 0.0, 0.0, 0.0, 2.0], inK, γ, 1e-3, 0.72)

KB.flux_gks!(fw, wL, wR, inK, γ, 1e-3, 0.72, dt, dx, dx, zeros(3), zeros(3))
KB.flux_gks!(
    zeros(4),
    [1.0, 0.3, 0.0, 1.0],
    [0.3, -0.1, 0.0, 2.0],
    inK,
    γ,
    1e-3,
    0.72,
    dt,
    dx,
    dx,
    dx,
    zeros(4),
    zeros(4),
) # 2D
KB.flux_gks!(
    zeros(4, 2),
    hcat([1.0, 0.3, 0.0, 1.0], [1.0, 0.3, 0.0, 1.0]),
    hcat([0.3, -0.1, 0.0, 2.0], [0.3, -0.1, 0.0, 2.0]),
    inK,
    γ,
    1.0,
    1.0,
    0.5,
    1.0,
    1e-2,
    dt,
    dx,
    dx,
    dx,
    zeros(4, 2),
    zeros(4, 2),
) # mixture
KB.flux_gks!(zeros(4), [1.0, 0.3, 0.0, 1.0], [0.3, -0.1, 0.0, 2.0], inK, γ, 1e-3, 0.72, dt) # FR

# discrete
KB.flux_gks!(fw, fh, wL, wR, u, inK, γ, 1e-3, 0.72, dt, dx, dx, zeros(3), zeros(3))
KB.flux_gks!(fw, fh, fb, wL, wR, u, inK, γ, 1e-3, 0.72, dt, dx, dx, zeros(3), zeros(3))
KB.flux_gks!(
    zeros(4),
    zeros(16, 16),
    [1.0, 0.3, 0.0, 1.0],
    [0.3, -0.1, 0.0, 2.0],
    rand(16, 16),
    rand(16, 16),
    γ,
    inK,
    1e-3,
    0.72,
    dt,
    dx,
    dx,
    dx,
    zeros(4),
    zeros(4),
) # 2D
KB.flux_gks!(
    zeros(4),
    zeros(16, 16),
    zeros(16, 16),
    [1.0, 0.3, 0.0, 1.0],
    [0.3, -0.1, 0.0, 2.0],
    rand(16, 16),
    rand(16, 16),
    γ,
    inK,
    1e-3,
    0.72,
    dt,
    dx,
    dx,
    dx,
    zeros(4),
    zeros(4),
) # 2D

KB.flux_ugks!(
    fw,
    fh,
    fb,
    wL,
    hL,
    bL,
    wR,
    hR,
    bR,
    u,
    ω,
    inK,
    γ,
    1e-3,
    0.81,
    0.72,
    dt,
    dx,
    dx,
)
KB.flux_ugks!(
    zeros(4),
    zeros(28, 28),
    zeros(28, 28),
    [1.0, 0.3, 0.0, 1.0],
    rand(28, 28),
    rand(28, 28),
    [0.3, -0.1, 0.0, 2.0],
    rand(28, 28),
    rand(28, 28),
    rand(28, 28),
    rand(28, 28),
    ones(28, 28),
    inK,
    γ,
    1e-3,
    0.81,
    0.72,
    dt,
    dx,
    dx,
    dx,
)
KB.flux_ugks!(
    zeros(5, 2),
    zeros(16, 16, 2),
    zeros(16, 16, 2),
    zeros(16, 16, 2),
    hcat([1.0, 0.3, 0.0, 0.0, 1.0], [1.0, 0.3, 0.0, 0.0, 1.0]),
    rand(16, 16, 2),
    rand(16, 16, 2),
    rand(16, 16, 2),
    hcat([0.3, -0.1, 0.0, 0.0, 2.0], [0.3, -0.1, 0.0, 0.0, 2.0]),
    rand(16, 16, 2),
    rand(16, 16, 2),
    rand(16, 16, 2),
    rand(16, 16, 2),
    rand(16, 16, 2),
    ones(16, 16, 2),
    inK,
    γ,
    1.0,
    1.0,
    0.5,
    1.0,
    0.1,
    dt,
    dx,
    dx,
    dx,
)

#--- kfvs flux ---#
KB.flux_kfvs(hL, hR, u, dt)
KB.flux_kfvs!(fh, hL, hR, u, dt)
KB.flux_kfvs!(fw, fh, hL, hR, u, ω, dt)
KB.flux_kfvs!(
    zeros(3, 2),
    zeros(16, 2),
    rand(16, 2),
    rand(16, 2),
    randn(16, 2),
    ones(16, 2),
    dt,
)
KB.flux_kfvs!(fw, fh, fb, hL, bL, hR, bR, u, ω, dt)
KB.flux_kfvs!(
    zeros(3, 2),
    zeros(16, 2),
    zeros(16, 2),
    rand(16, 2),
    rand(16, 2),
    rand(16, 2),
    rand(16, 2),
    randn(16, 2),
    ones(16, 2),
    dt,
)
KB.flux_kfvs!(zeros(4), fh, fb, zero(fh), hL, bL, bL, hR, bR, bR, u, ω, dt) # Rykov
KB.flux_kfvs!(
    zeros(5),
    zeros(16, 16, 16),
    rand(16, 16, 16),
    rand(16, 16, 16),
    randn(16, 16, 16),
    randn(16, 16, 16),
    randn(16, 16, 16),
    ones(16, 16, 16),
    dt,
)
KB.flux_kfvs!(
    zeros(5),
    zeros(16, 16, 16),
    rand(16, 16, 16),
    rand(16, 16, 16),
    randn(16, 16, 16),
    randn(16, 16, 16),
    randn(16, 16, 16),
    ones(16, 16, 16),
    dt,
    1.0,
)
KB.flux_kfvs!(
    zeros(5),
    zeros(16),
    zeros(16),
    zeros(16),
    zeros(16),
    rand(16),
    rand(16),
    rand(16),
    rand(16),
    rand(16),
    rand(16),
    rand(16),
    rand(16),
    randn(16),
    ones(16),
    dt,
)
KB.flux_kfvs!(
    zeros(5, 2),
    zeros(16, 2),
    zeros(16, 2),
    zeros(16, 2),
    zeros(16, 2),
    rand(16, 2),
    rand(16, 2),
    rand(16, 2),
    rand(16, 2),
    rand(16, 2),
    rand(16, 2),
    rand(16, 2),
    rand(16, 2),
    randn(16, 2),
    ones(16, 2),
    dt,
)
KB.flux_kfvs!(
    zeros(4),
    zeros(16, 16),
    rand(16, 16),
    rand(16, 16),
    randn(16, 16),
    randn(16, 16),
    ones(16, 16),
    dt,
    dx,
)
KB.flux_kfvs!(
    zeros(4),
    zeros(16, 16),
    zeros(16, 16),
    rand(16, 16),
    rand(16, 16),
    rand(16, 16),
    rand(16, 16),
    randn(16, 16),
    randn(16, 16),
    ones(16, 16),
    dt,
    dx,
)
# 3F2V
KB.flux_kfvs!(
    zeros(5),
    zeros(16, 16),
    zeros(16, 16),
    zeros(16, 16),
    rand(16, 16),
    rand(16, 16),
    rand(16, 16),
    rand(16, 16),
    rand(16, 16),
    rand(16, 16),
    randn(16, 16),
    randn(16, 16),
    ones(16, 16),
    dt,
    dx,
)
KB.flux_kfvs!(
    zeros(5, 2),
    zeros(16, 16, 2),
    zeros(16, 16, 2),
    zeros(16, 16, 2),
    rand(16, 16, 2),
    rand(16, 16, 2),
    rand(16, 16, 2),
    rand(16, 16, 2),
    rand(16, 16, 2),
    rand(16, 16, 2),
    randn(16, 16, 2),
    randn(16, 16, 2),
    ones(16, 16, 2),
    dt,
    dx,
)

#--- kinetic central-upwind flux ---#
# 1F1V
KB.flux_kcu!(fw, fh, wL, hL, wR, hR, u, ω, inK, γ, 1e-3, 0.81, 0.72, dt)
KB.flux_kcu!(
    hcat(fw, fw),
    hcat(fh, fh),
    hcat(wL, wL),
    hcat(hL, hL),
    hcat(wR, wR),
    hcat(hR, hR),
    hcat(u, u),
    hcat(ω, ω),
    inK,
    γ,
    1.0,
    1.0,
    0.5,
    1.0,
    0.01,
    dt,
)
# 2F1V
KB.flux_kcu!(fw, fh, fb, wL, hL, bL, wR, hR, bR, u, ω, inK, γ, 1e-3, 0.81, 0.72, dt)
KB.flux_kcu!(
    zeros(3, 2),
    zeros(16, 2),
    zeros(16, 2),
    zeros(3, 2),
    rand(16, 2),
    rand(16, 2),
    zeros(3, 2),
    rand(16, 2),
    rand(16, 2),
    rand(16, 2),
    ones(16, 2),
    inK,
    γ,
    1.0,
    1.0,
    0.5,
    1.0,
    0.01,
    dt,
)
# 4F1V
KB.flux_kcu!(
    zeros(5),
    zeros(16),
    zeros(16),
    zeros(16),
    zeros(16),
    zeros(5),
    rand(16),
    rand(16),
    rand(16),
    rand(16),
    zeros(5),
    rand(16),
    rand(16),
    rand(16),
    rand(16),
    rand(16),
    ones(16),
    inK,
    γ,
    1e-3,
    0.81,
    0.72,
    dt,
)
KB.flux_kcu!(
    zeros(5, 2),
    zeros(16, 2),
    zeros(16, 2),
    zeros(16, 2),
    zeros(16, 2),
    zeros(5, 2),
    rand(16, 2),
    rand(16, 2),
    rand(16, 2),
    rand(16, 2),
    zeros(5, 2),
    rand(16, 2),
    rand(16, 2),
    rand(16, 2),
    rand(16, 2),
    rand(16, 2),
    ones(16, 2),
    inK,
    γ,
    1.0,
    1.0,
    0.5,
    1.0,
    0.01,
    dt,
)
# 1F2V
KB.flux_kcu!(
    zeros(4),
    zeros(16, 16),
    zeros(4),
    rand(16, 16),
    zeros(4),
    rand(16, 16),
    rand(16, 16),
    rand(16, 16),
    ones(16, 16),
    inK,
    γ,
    1e-3,
    0.81,
    0.72,
    dt,
    dx,
)
# 2F2V
KB.flux_kcu!(
    zeros(4),
    zeros(16, 16),
    zeros(16, 16),
    zeros(4),
    rand(16, 16),
    rand(16, 16),
    zeros(4),
    rand(16, 16),
    rand(16, 16),
    rand(16, 16),
    rand(16, 16),
    ones(16, 16),
    inK,
    γ,
    1e-3,
    0.81,
    0.72,
    dt,
    dx,
)
# 3F2V
KB.flux_kcu!(
    zeros(5),
    zeros(16, 16),
    zeros(16, 16),
    zeros(16, 16),
    zeros(5),
    rand(16, 16),
    rand(16, 16),
    rand(16, 16),
    zeros(5),
    rand(16, 16),
    rand(16, 16),
    rand(16, 16),
    rand(16, 16),
    rand(16, 16),
    ones(16, 16),
    inK,
    γ,
    1e-3,
    0.81,
    0.72,
    dt,
    dx,
)
KB.flux_kcu!(
    zeros(5, 2),
    zeros(16, 16, 2),
    zeros(16, 16, 2),
    zeros(16, 16, 2),
    zeros(5, 2),
    rand(16, 16, 2),
    rand(16, 16, 2),
    rand(16, 16, 2),
    zeros(5, 2),
    rand(16, 16, 2),
    rand(16, 16, 2),
    rand(16, 16, 2),
    rand(16, 16, 2),
    rand(16, 16, 2),
    ones(16, 16, 2),
    inK,
    γ,
    1.0,
    1.0,
    0.5,
    1.0,
    0.01,
    dt,
    dx,
)

#--- boundary flux ---#
KB.flux_boundary_maxwell!(zeros(3), rand(3), [1.0, 0.0, 1.0], 2.0, 5 / 3, 1e-3, 1)
KB.flux_boundary_maxwell!(
    zeros(4),
    rand(4),
    [1.0, 0.0, 0.0, 1.0],
    1.0,
    5 / 3,
    1e-3,
    1e-2,
    1,
)
KB.flux_boundary_maxwell!(
    zeros(3),
    zeros(16),
    zeros(16),
    [1.0, 0.3, 0.89],
    rand(16),
    rand(16),
    randn(16),
    ones(16),
    inK,
    dt,
    1,
)
KB.flux_boundary_maxwell!(
    zeros(4),
    zeros(16, 16),
    zeros(16, 16),
    [1.0, 0.3, -0.2, 0.89],
    rand(16, 16),
    rand(16, 16),
    randn(16, 16),
    randn(16, 16),
    ones(16, 16),
    inK,
    dt,
    dx,
    1,
)
KB.flux_boundary_maxwell!(
    zeros(5),
    zeros(16, 16, 16),
    [1.0, 0.3, 0.1, 0.1, 0.89],
    rand(16, 16, 16),
    randn(16, 16, 16),
    randn(16, 16, 16),
    randn(16, 16, 16),
    ones(16, 16, 16),
    dt,
    1,
    1,
)
# mixture
KB.flux_boundary_maxwell!(
    zeros(3, 2),
    zeros(16, 2),
    zeros(16, 2),
    hcat([1.0, 0.3, 0.89], [1.0, 0.3, 0.89]),
    rand(16, 2),
    rand(16, 2),
    randn(16, 2),
    ones(16, 2),
    2,
    dt,
    1,
)

KB.flux_boundary_specular!(zeros(3), zeros(16), rand(16), randn(16), ones(16), dt)
KB.flux_boundary_specular!(
    zeros(3),
    zeros(16),
    zeros(16),
    rand(16),
    rand(16),
    randn(16),
    ones(16),
    dt,
)

#--- pure equilibrium flux ---#
KB.flux_equilibrium!(fw, wL, wR, inK, γ, 1e-3, 0.81, 0.72, dt, dx, dx, zeros(3), zeros(3))
KB.flux_equilibrium!(
    zeros(4),
    zeros(4),
    zeros(4),
    inK,
    γ,
    1e-3,
    0.81,
    0.72,
    dt,
    dx,
    dx,
    dx,
    zeros(4),
    zeros(4),
)

#--- electromagnetic flux ---#
sol = 100.0
χ = 1.0
ν = 1.0
# A^+
Ap = Array{Float64}(undef, 8, 8)
Ap[1, 1] = (sol * χ) / 2.0
Ap[7, 1] = χ / 2.0
Ap[2, 2] = sol / 2.0
Ap[6, 2] = 0.5
Ap[3, 3] = sol / 2.0
Ap[5, 3] = -1.0 / 2.0
Ap[4, 4] = (sol * ν) / 2.0
Ap[8, 4] = (sol^2 * ν) / 2.0
Ap[3, 5] = -sol^2 / 2.0
Ap[5, 5] = sol / 2.0
Ap[2, 6] = sol^2 / 2.0
Ap[6, 6] = sol / 2.0
Ap[1, 7] = (sol^2 * χ) / 2.0
Ap[7, 7] = (sol * χ) / 2.0
Ap[4, 8] = ν / 2.0
Ap[8, 8] = (sol * ν) / 2.0
# A^-
An = Array{Float64}(undef, 8, 8)
An[1, 1] = -(sol * χ) / 2.0
An[7, 1] = χ / 2.0
An[2, 2] = -sol / 2.0
An[6, 2] = 1.0 / 2.0
An[3, 3] = -sol / 2.0
An[5, 3] = -1.0 / 2.0
An[4, 4] = -(sol * ν) / 2.0
An[8, 4] = (sol^2 * ν) / 2.0
An[3, 5] = -sol^2 / 2.0
An[5, 5] = -sol / 2.0
An[2, 6] = sol^2 / 2.0
An[6, 6] = -sol / 2.0
An[1, 7] = (sol^2 * χ) / 2.0
An[7, 7] = -(sol * χ) / 2.0
An[4, 8] = ν / 2.0
An[8, 8] = -(sol * ν) / 2.0
# eigenvalues
D = [sol, sol, sol * χ, sol * ν, -sol, -sol, -sol * χ, -sol * ν]

KB.flux_em!(
    zeros(8),
    zeros(8),
    randn(3),
    randn(3),
    randn(3),
    randn(3),
    randn(3),
    randn(3),
    randn(3),
    randn(3),
    0.0,
    0.0,
    0.0,
    0.0,
    dx,
    dx,
    Ap,
    An,
    D,
    100.0,
    1.0,
    1.0,
    dt,
)
KB.flux_emx!(
    zeros(8),
    zeros(8),
    zeros(8),
    zeros(8),
    zeros(8),
    zeros(8),
    randn(3),
    randn(3),
    randn(3),
    randn(3),
    randn(3),
    randn(3),
    randn(3),
    randn(3),
    0.0,
    0.0,
    0.0,
    0.0,
    dx,
    dx,
    Ap,
    An,
    Ap,
    An,
    D,
    100.0,
    1.0,
    1.0,
    dt,
)
KB.flux_emy!(
    zeros(8),
    zeros(8),
    zeros(8),
    zeros(8),
    zeros(8),
    zeros(8),
    randn(3),
    randn(3),
    randn(3),
    randn(3),
    randn(3),
    randn(3),
    randn(3),
    randn(3),
    0.0,
    0.0,
    0.0,
    0.0,
    dx,
    dx,
    Ap,
    An,
    Ap,
    An,
    D,
    100.0,
    1.0,
    1.0,
    dt,
)
