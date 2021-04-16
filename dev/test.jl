import KitBase

inK = 0
γ = 3

primL = [1.0, 0.0, 1.0]
primR = [0.125, 0.0, 0.625]
wL = KitBase.prim_conserve(primL, γ)
wR = KitBase.prim_conserve(primR, γ)
Δx = 0.01

a = KitBase.pdf_slope(primL, (wR .- wL) ./ Δx, inK)

Mu, Mxi, MuL, MuR = KitBase.gauss_moments(primL, inK)
Mau = KitBase.moments_conserve_slope(a, Mu, Mxi, 1)
A = KitBase.pdf_slope(primL, -primL[1] .* Mau, inK)

vs = KitBase.VSpace1D()
M = KitBase.maxwellian(vs.u, primL)
τ = 0.01
#f = @. M * (1 - τ * (a[1] * vs.u + a[2] * vs.u^2 + 0.5 * a[3] * vs.u^3 + A[1] + A[2] * vs.u + 0.5 * A[3] * vs.u^2))
f = KitBase.chapman_enskog(vs.u, primL, a, A, τ)

KitBase.moments_conserve(f, vs.u, vs.weights)
KitBase.moments_conserve(M, vs.u, vs.weights)

inK = 0
γ = 2

primL = [1.0, 0.0, 0.0, 1.0]
primR = [0.125, 0.0, 0.0, 0.625]
wL = KitBase.prim_conserve(primL, γ)
wR = KitBase.prim_conserve(primR, γ)

a = KitBase.pdf_slope(primL, (wR .- wL) ./ Δx, inK)
b = deepcopy(a)
Mu, Mv, Mxi, MuL, MuR = KitBase.gauss_moments(primL, inK)
Mau = KitBase.moments_conserve_slope(a, Mu, Mv, Mxi, 1, 0)
Mbu = KitBase.moments_conserve_slope(b, Mu, Mv, Mxi, 0, 1)
A = KitBase.pdf_slope(primL, -primL[1] .* (Mau .+ Mbu), inK)

vs = KitBase.VSpace2D()
M = KitBase.maxwellian(vs.u, vs.v, primL)
f = KitBase.chapman_enskog(vs.u, vs.v, primL, a, b, A, τ)

KitBase.moments_conserve(f, vs.u, vs.v, vs.weights)
KitBase.moments_conserve(M, vs.u, vs.v, vs.weights)
