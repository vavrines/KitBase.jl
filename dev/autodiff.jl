"""
Develop automatic differentiation (AD) support in Kinetic.jl
"""

using KitBase, ForwardDiff, ReverseDiff, BenchmarkTools, Plots

prim = [1.0, 0.0, 1.0]
vs = VSpace1D(-5, 5, 100)

# We build a unary Maxwellian function for testing.
Mu(u) = maxwellian(u, prim)

@btime fd = ForwardDiff.jacobian(Mu, vs.u)
@btime rd = ReverseDiff.jacobian(Mu, vs.u)
"""
On my workstation, the benchmark results are:
- Forward mode: 30.569 μs (15 allocations: 183.64 KiB)
- Reverse mode: 202.087 ms (6501013 allocations: 215.04 MiB)
So the foward mode seems to be much more efficient.
"""

@btime fd1 = ForwardDiff.derivative.(Mu, vs.u)
"""
15.390 μs (313 allocations: 9.15 KiB)
Using derivative is faster.
"""

# An advantage is that the reverse mode is not restricted to unary functions
@btime rd2 = ReverseDiff.jacobian(maxwellian, (vs.u, prim))
"""
405.675 ms (12996209 allocations: 429.89 MiB)
But it's even slower.
"""

# Check results
fd = ForwardDiff.jacobian(Mu, vs.u)
rd = ReverseDiff.jacobian(Mu, vs.u)
rd2 = ReverseDiff.jacobian(maxwellian, (vs.u, prim))
fd == rd == rd2[1]

dm = [rd2[1][i, i] for i in axes(fd, 1)]
plot(vs.u, dm)

#　These things can be built into a function

function ∂maxwellian(u::Real, ρ, U, λ)
    Mu = u -> maxwellian(u, ρ, U, λ)
    return ForwardDiff.derivative(Mu, u)
end

@btime ∂maxwellian.(vs.u, prim[1], prim[2], prim[3])

function ∂maxwellian(u::AbstractVector, prim)
    Mu = u -> maxwellian(u, prim)
    return ForwardDiff.derivative.(Mu, u)
end

@btime ∂maxwellian(vs.u, prim)
