import KitBase

inK = 2
γ = 5 / 3
primL = [1.0, 0.0, 1.0]
primR = [0.125, 0.0, 0.625]
wL = KitBase.prim_conserve(primL, γ)
wR = KitBase.prim_conserve(primR, γ)
Δx = 0.01

a = KitBase.pdf_slope(primL, (wR .- wL) ./ Δx, inK)

Mu, Mxi, MuL, MuR = KitBase.gauss_moments(primL, inK)
Mau = KitBase.moments_conserve_slope(a, Mu, Mxi, 1)
A = KitBase.pdf_slope(primL, -primL[1] .* Mau, inK)

u = -5:0.1:5 |> collect
M = KitBase.maxwellian(u, primL)
τ = 0.01
f = @. M * (1 - τ * (a[1] * u + a[2] * u^2 + 0.5 * a[3] * u^3 + A[1] + A[2] * u + 0.5 * A[3] * u^2))
f = KitBase.chapman_enskog(u, primL, a, A, τ)

sum(M)
sum(f)

cd(@__DIR__)
set, ctr, face, t = KitBase.initialize("config.toml")