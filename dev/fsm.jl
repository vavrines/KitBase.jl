using OrdinaryDiffEq, Plots
using KitBase

begin
    case = "homogeneous"
    maxTime = 3
    tlen = 16
    u0 = -5
    u1 = 5
    nu = 80
    nug = 0
    v0 = -5
    v1 = 5
    nv = 28
    nvg = 0
    w0 = -5
    w1 = 5
    nw = 28
    nwg = 0
    vMeshType = "rectangle"
    nm = 5
    knudsen = 1
    inK = 0
    alpha = 1.0
    omega = 0.5
    nh = 8
end

tspan = (0.0, maxTime)
tran = linspace(tspan[1], tspan[2], tlen)
γ = heat_capacity_ratio(inK, 3)
vSpace = VSpace3D(u0, u1, nu, v0, v1, nv, w0, w1, nw, vMeshType)

f0 =
    Float32.(
        0.5 * (1 / π)^1.5 .*
        (exp.(-(vSpace.u .- 0.99) .^ 2) .+ exp.(-(vSpace.u .+ 0.99) .^ 2)) .*
        exp.(-vSpace.v .^ 2) .* exp.(-vSpace.w .^ 2),
    ) |> Array
prim0 =
    conserve_prim(moments_conserve(f0, vSpace.u, vSpace.v, vSpace.w, vSpace.weights), γ)
M0 = Float32.(maxwellian(vSpace.u, vSpace.v, vSpace.w, prim0)) |> Array

mu_ref = ref_vhs_vis(knudsen, alpha, omega)
kn_bzm = hs_boltz_kn(mu_ref, 1.0)
τ0 = mu_ref * 2.0 * prim0[end]^(0.5) / prim0[1]

phi, psi, phipsi = kernel_mode(
    nm,
    vSpace.u1,
    vSpace.v1,
    vSpace.w1,
    vSpace.du[1, 1, 1],
    vSpace.dv[1, 1, 1],
    vSpace.dw[1, 1, 1],
    vSpace.nu,
    vSpace.nv,
    vSpace.nw,
    alpha,
)

# Boltzmann
prob = ODEProblem(boltzmann_ode!, f0, tspan, [kn_bzm, nm, phi, psi, phipsi])
data_boltz = solve(prob, Tsit5(), saveat = tran) |> Array

# BGK
prob1 = ODEProblem(bgk_ode!, f0, tspan, [M0, τ0])
data_bgk = solve(prob1, Tsit5(), saveat = tran) |> Array

plot(vSpace.u[:, 1, 1], data_boltz[:, end÷2, end÷2, end], label="FSM")
plot!(vSpace.u[:, 1, 1], data_bgk[:, end÷2, end÷2, end], label="BGK")

###

function KitBase.kernel_mode(
    M::I,
    umax::R,
    vmax::R,
    wmax::R,
    du::AbstractArray,
    dv::AbstractArray,
    dw::AbstractArray,
    unum::I,
    vnum::I,
    wnum::I,
    alpha::R;
    quad_num = 64::I,
) where {I<:Integer,R<:Real}

    supp = sqrt(2.0) * 2.0 * max(umax, vmax, wmax) / (3.0 + sqrt(2.0))

    #fre_vx = collect(range(-π, (unum ÷ 2 - 1) * 2.0 * π / unum, length = unum)) ./ du
    #fre_vy = collect(range(-π, (vnum ÷ 2 - 1) * 2.0 * π / vnum, length = vnum)) ./ dv
    #fre_vz = collect(range(-π, (wnum ÷ 2 - 1) * 2.0 * π / wnum, length = wnum)) ./ dw

    fre_vx = [i * π / umax for i = -unum÷2:unum÷2-1]
    fre_vy = [i * π / vmax for i = -vnum÷2:vnum÷2-1]
    fre_vz = [i * π / wmax for i = -wnum÷2:wnum÷2-1]

    abscissa, gweight = KitBase.lgwt(quad_num, 0, supp)

    phi = zeros(unum, vnum, wnum, M * (M - 1))
    psi = zeros(unum, vnum, wnum, M * (M - 1))
    phipsi = zeros(unum, vnum, wnum)
    for loop = 1:M-1
        theta = π / M * loop
        for loop2 = 1:M
            theta2 = π / M * loop2
            idx = (loop - 1) * M + loop2
            for k = 1:wnum, j = 1:vnum, i = 1:unum
                s =
                    fre_vx[i] * sin(theta) * cos(theta2) +
                    fre_vy[j] * sin(theta) * sin(theta2) +
                    fre_vz[k] * cos(theta)
                # phi
                int_temp = 0.0
                for id = 1:quad_num
                    int_temp +=
                        2.0 * gweight[id] * cos(s * abscissa[id]) * (abscissa[id]^alpha)
                end
                phi[i, j, k, idx] = int_temp * sin(theta)
                # psi
                s = fre_vx[i]^2 + fre_vy[j]^2 + fre_vz[k]^2 - s^2
                if s <= 0.0
                    psi[i, j, k, idx] = π * supp^2
                else
                    s = sqrt(s)
                    bel = supp * s
                    bessel = besselj(1, bel)
                    psi[i, j, k, idx] = 2.0 * π * supp * bessel / s
                end
                # phipsi
                phipsi[i, j, k] += phi[i, j, k, idx] * psi[i, j, k, idx]
            end
        end
    end

    return phi, psi, phipsi

end

vs = VSpace3D(u0, u1, nu, v0, v1, nv, w0, w1, nw, "algebra")
vs = VSpace3D(u0, u1, nu, v0, v1, nv, w0, w1, nw, "rectangle")

phi, psi, phipsi = kernel_mode(
    nm,
    vs.u1,
    vs.v1,
    vs.w1,
    vs.du[:, 1, 1],
    vs.dv[1, :, 1],
    vs.dw[1, 1, :],
    vs.nu,
    vs.nv,
    vs.nw,
    alpha,
)

f0 =
    Float32.(
        0.5 * (1 / π)^1.5 .*
        (exp.(-(vs.u .- 0.99) .^ 2) .+ exp.(-(vs.u .+ 0.99) .^ 2)) .*
        exp.(-vs.v .^ 2) .* exp.(-vs.w .^ 2),
    ) |> Array
prim0 =
    conserve_prim(moments_conserve(f0, vs.u, vs.v, vs.w, vs.weights), γ)
M0 = Float32.(maxwellian(vs.u, vs.v, vs.w, prim0)) |> Array

mu_ref = ref_vhs_vis(knudsen, alpha, omega)
kn_bzm = hs_boltz_kn(mu_ref, 1.0)
τ0 = mu_ref * 2.0 * prim0[end]^(0.5) / prim0[1]

# Boltzmann
prob = ODEProblem(boltzmann_ode!, f0, tspan, [kn_bzm, nm, phi, psi, phipsi])
data_boltz = solve(prob, Tsit5(), saveat = tran) |> Array

plot(vs.u[:, 1, 1], data_boltz[:, end÷2, end÷2, end], label="FSM2")

scatter!(vs.u[:, 1, 1], data_boltz[:, end÷2, end÷2, end], label="FSM2")

plot!(vs.u[:, 1, 1], data_boltz[:, end÷2, end÷2, end], label="FSM2")