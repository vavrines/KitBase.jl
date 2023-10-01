# ============================================================
# Ref: Polylogarithms.jl
# M. Roughan. The Polylogarithm Function in Julia. arXiv:2010.09860.
# ============================================================

const ComplexOrReal{T} = Union{T,Complex{T}}
const default_accuracy = 1.0e-12
const default_max_iterations = 1000
const near_int_threshold = 1.0e-6
const series_transition_threshold = 0.25
struct Diagnostics end

"""
$(SIGNATURES)

Calculate the Polylogarithm function `Li_s(z)` defined by

``{Li}_s = \\sum_{n=1}^{\\infty} \\frac{z^n}{n^s}``

It uses double precision complex numbers (not arbitrary precision).
The goal is an relative error bound 10^{-12}.
 
## Arguments
* `s::Complex`: the 'fractional' parameter
* `z::Complex`: the point at which to calculate it
It should also accept input arguments as `Real` or `Rational` or `Integer` but these aren't completely tested.
There are additional keywords, but these are currently intended for testing not use.

## Outputs
* `Li_s(z)`: the result

## Examples
```
julia> polylog(0.35, 0.2)
0.23803890574407033
```
"""
function polylog(
    s::Number,
    z::Number;
    level = 1, # keep track of recursion
    accuracy::Float64 = default_accuracy,
    min_iterations::Integer = 0,
    max_iterations::Integer = default_max_iterations,
)
    return polylog(
        s,
        z,
        Diagnostics();
        level = level,
        accuracy = accuracy,
        min_iterations = min_iterations,
        max_iterations = max_iterations,
    ) |> first # first is the result
end

"""
$(SIGNATURES)

Calculate the Polylogarithm function ``{Li}_s(z)`` defined by

``{Li}_s = \\sum_{n=1}^{\\infty} \\frac{z^n}{n^s}``

It uses double precision complex numbers (not arbitrary precision).
The goal is an relative error bound 10^{-12}.
 
## Arguments
* `s::Complex`: the 'fractional' parameter
* `z::Complex`: the point at which to calculate it
* `::Diagnostics`: use this to indicate that the output should include extra information
It should also accept input arguments as `Real` or `Rational` or `Integer` but these aren't completely tested.
There are additional keywords, but these are currently intended for testing not use.

## Outputs
* `Li_s(z)`: the result
* `n`: the number of elements used in each series
* `series`: the series used to compute results (4 = reciprocal)

## Examples
```
julia> polylog(0.35, 0.2, Diagnostics() )
(0.23803890574407033, 17, 1)
```
"""
function polylog(
    s::Number,
    z::Number,
    ::Diagnostics;
    level = 1, # keep track of recursion
    accuracy::Float64 = default_accuracy,
    min_iterations::Integer = 0,
    max_iterations::Integer = default_max_iterations,
)
    tau_threshold = 1.0e-3
    μ = log(convert(Complex{Float64}, z)) # input z could be an integer or anything
    t = abs(μ / 2 / π)
    T = 0.512 # the duplication formula seems to work if T=0.5 in the small associated wedge, but wyh risk it?
    if abs(μ) < 1.0e-14
        # deal with this case separately or it causes pain in Series 2 and 3, which should be resolvable, but its just easier here
        # there are lots of other special cases to add here eventually
        if real(s) > 1
            return SpecialFunctions.zeta(s), 0, 0
        else
            return typeof(z)(Inf), 0, 0
        end
    elseif abs(z) <= 0.5 && abs(z) < t
        return polylog_series_1(
            s,
            z;
            accuracy = accuracy,
            min_iterations = min_iterations,
            max_iterations = max_iterations,
        )
    elseif t <= T && (abs(round(real(s)) - s) > tau_threshold || real(s) <= 0)
        return polylog_series_2(
            s,
            z;
            accuracy = accuracy,
            min_iterations = min_iterations,
            max_iterations = max_iterations,
        )
    elseif t <= T
        return polylog_series_3(
            s,
            z;
            accuracy = accuracy,
            min_iterations = min_iterations,
            max_iterations = max_iterations,
        )
    else
        return polylog_duplication(
            s,
            z;
            level = level,
            accuracy = accuracy,
            min_iterations = min_iterations,
            max_iterations = max_iterations,
        )
    end
    # we could have a lot more special cases here, particularly for integer input to make the code faster for these cases
end


"""
$(SIGNATURES)

calculate using the duplication formula
"""
function polylog_duplication(
    s::Number,
    z::Number;
    level = 0, # keep track of recursion
    accuracy::Float64 = default_accuracy,
    min_iterations::Integer = 0,
    max_iterations::Integer = default_max_iterations,
)
    z = convert(Complex{Float64}, z)
    f = min(0.5, 2.0^(1 - real(s)))
    # println("  dup level $level, z=$z,  abs(μ)/2π = ", abs(log(z))/twoπ  )
    (Li1, k1, series1) = polylog(
        s,
        sqrt(z),
        Diagnostics();
        level = level + 1,
        accuracy = f * accuracy,
        min_iterations = min_iterations,
        max_iterations = max_iterations,
    )
    (Li2, k2, series2) = polylog(
        s,
        -sqrt(z),
        Diagnostics();
        level = level + 1,
        accuracy = f * accuracy,
        min_iterations = min_iterations,
        max_iterations = max_iterations,
    )
    if typeof(s) <: Real
        s = convert(Float64, s) # convert s into a double
    elseif typeof(s) <: Complex
        s = convert(Complex{Float64}, s) # convert s into doubles
    end
    return (2^(s - 1) * (Li1 + Li2), k1 + k2, 10 + series1 + series2)
end


"""
$(SIGNATURES)

Calculate using direct definition
"""
function polylog_series_1(
    s::Number,
    z::Number;
    accuracy::Float64 = default_accuracy,
    min_iterations::Integer = 0,
    max_iterations::Integer = default_max_iterations,
    existing_total::Number = 0.0,
)
    # prolly should convert z to a double or complex-double
    if abs(z) > 1 || (abs(z) ≈ 1 && real(s) <= 2)
        throw(DomainError(z))
    end
    if abs(z) > 1 / 2
        throw(DomainError(z, "Slow convergence for  |z| > 1/2"))
    end
    total = 0.0
    converged = false
    a = z
    k = 0
    if real(s) < 0
        min_iterations = ceil(real(s) / log(abs(z)))
    end
    while k <= max_iterations && ~converged
        k = k + 1
        total += a
        a *= z * (k / (k + 1.0))^s
        # println("   total = $total")
        if k > min_iterations && abs(a) / abs(total) < 0.5 * accuracy
            converged = true
        end
    end
    return (total, k, 1)
end


"""
$(SIGNATURES)

Calculate using power series around μ = log(z) = 0

This should not be used near positive integer values of s.
"""
function polylog_series_2(
    s::Number,
    z::Number;
    accuracy::Float64 = default_accuracy,
    min_iterations::Integer = 0,
    max_iterations::Integer = default_max_iterations,
)
    μ = log(convert(Complex{Float64}, z)) # input z could be an integer or anything
    if typeof(s) <: Real
        s = convert(Float64, s) # convert s into a double
    elseif typeof(s) <: Complex
        s = convert(Complex{Float64}, s) # convert s into doubles
    end

    if abs(μ) > 2π
        throw(DomainError(z, "we need |log(z)|< 2π for this series"))
    end

    oneminuss = 1.0 - s
    total = SpecialFunctions.gamma(oneminuss) * (-μ)^(-oneminuss)

    converged = false
    tmp = 1
    k = 0
    a = Inf
    a_2 = Inf

    while k <= max_iterations && ~converged
        a_3 = a_2
        a_2 = a
        a = tmp * SpecialFunctions.zeta(s - k)
        total += a
        tmp *= μ / (k + 1)
        if k > min_iterations &&
           abs(a) / abs(total) < 0.5 * accuracy &&
           abs(a_2) / abs(total) < 0.5 * accuracy &&
           abs(a_2) > abs(a)
            # abs( A * (k-real(s))^real(-s) * (μ/twoπ)^k )/abs(total) < accuracy
            # && abs( 2*twoπ^real(s) * (μ/twoπ)^k )/abs(total) < accuracy 
            # the stopping rule should be implemented more efficiently as part of the calculation above
            converged = true
        end
        k = k + 1
    end

    # get correct value along the branch
    if isreal(z) && real(z) >= 1
        # total -= 2*π*im*μ^(s-1)/SpecialFunctions.gamma(s)
        total -= exp(log(2π * im) + (s - 1) * log(μ) - SpecialFunctions.loggamma(s))
    end
    return (total, k, 2)
end

function c_closed(n::Integer, j::Integer, ℒ::Number)
    d2 = SpecialFunctions.digamma(n + 1) - ℒ
    if j == 0
        return harmonic(n) - ℒ
    elseif j == 1
        # Wood:+stieltjes(1) - d2^2/2 + π^2/6 + SpecialFunctions.polygamma(1,n+1)/2 
        return -stieltjes(1) - d2^2 / 2 - π^2 / 6 + SpecialFunctions.polygamma(1, n + 1) / 2
    elseif j == 2
        # Wood:stieltjes(2)   + d2^3/6 + d2*( π^2/6 + SpecialFunctions.polygamma(1,n+1)/2 ) + SpecialFunctions.polygamma(2,n+1)/6
        return stieltjes(2) / 2 +
               d2^3 / 6 +
               d2 * (π^2 / 6 - SpecialFunctions.polygamma(1, n + 1) / 2) +
               SpecialFunctions.polygamma(2, n + 1) / 6
    end
end

function Q_closed(n::Integer, τ::Number, ℒ::Number; n_terms::Integer = 3)
    # τ is the distance from the pole s=n>0, ℒ = log(-μ) = log(-log( z ))
    max_n_terms = 3
    if n_terms < 1 || n_terms > max_n_terms
        throw(DomainError(n_terms))
    end
    return sum(c_closed.(n, 0:n_terms-1, ℒ) .* τ .^ (0:n_terms-1))
end

function Q(n::Integer, τ::Number, ℒ::Number; n_terms::Integer = 5) # Crandall,2012, p.35
    # τ is the distance from the pole s=n>0, ℒ = log(-μ) = log(-log( z ))
    if abs(τ) <= 1.0e-14
        # if really close to the integer, then ignore the extra terms
        return c_closed(n, 0, ℒ)
    else
        max_n_terms = 7
        if n_terms < 1 || n_terms > max_n_terms
            throw(DomainError(n_terms))
        end
        if n_terms <= 3
            # use the direct method in this case
            return Q_closed(n, τ, ℒ; n_terms = n_terms)
        end
        return sum(c_crandall.(n, 0:n_terms-1, ℒ) .* τ .^ (0:n_terms-1))
    end
end


"""
$(SIGNATURES)

Crandall, 2012, p.35
"""
function c_crandall(k::Integer, j::Integer, ℒ)
    return (-1)^j * stieltjes(j) / SpecialFunctions.gamma(j + 1) - b_crandall(k, j + 1, ℒ)
end

"""
$(SIGNATURES)

Crandall, 2012, p.36
"""
function b_crandall(k::Integer, j::Integer, ℒ)
    total = 0
    for q = 0:j
        for t = 0:j-q
            p = j - q - t
            a1 = ℒ^p / SpecialFunctions.gamma(p + 1)
            a2 = (-1)^t * f_crandall(k, q) # Bailey and Borwein, 2015 correct Crandall (t+q - > t)
            a3 = g_crandall(t) / SpecialFunctions.gamma(t + 1)
            total += a1 * a2 * a3
        end
    end
    return total
end

const gamma_t = [
    1.0,
    -0.5772156649015315,
    1.9781119906559432,
    -5.44487445648531,
    23.561474084025583,
    -117.83940826837748,
    715.0673625273184,
    -5019.848872629852,
    40243.62157333573,
    -362526.2891146549,
    3.627042412756892e6,
    -3.990708415143132e7,
    4.7894329176518273e8,
    -6.226641351546061e9,
    8.717563381070836e10,
]
"""
$(SIGNATURES)

Crandall,2012, p.17
"""
function g_crandall(t::Integer)
    # t derivate of Gamma function at 1
    # see "gamma_derivatives.jl" for derivations of these numbers
    if t < 0
        throw(DomainError(t))
    elseif t > 14
        throw(DomainError(t, "only calculate the 1st 14 derivatives"))
    else
        return gamma_t[t+1]
    end
end

"""
$(SIGNATURES)

Crandall,2012, p.36
"""
function f_crandall(k::Integer, q::Integer)
    # Crandall omits the k=0 case. but this is given in Bailey and Borwein and matches other text
    if k == 0 && q == 0
        return 1
    elseif k == 0
        return 0
    elseif q == 0
        return 1
    else
        h = 0:q
        return sum((-1) .^ h .* f_crandall.(k - 1, q .- h) ./ k .^ h)
    end
end


"""
$(SIGNATURES)

For the special case that s is near a postive integer n>0
Calculate in a power series around `z=1`, and `s=n`
"""
function polylog_series_3(
    s::Number,
    z::Number;
    accuracy::Float64 = default_accuracy,
    min_iterations::Integer = 0,
    max_iterations::Integer = default_max_iterations,
    n_terms::Integer = 5,
)
    μ = log(convert(Complex{Float64}, z))
    if abs(μ) > 2π
        throw(DomainError(z, "does not converge for abs(ln(z)) > 2π"))
    end
    if real(s) <= 0.5
        throw(DomainError(s, "for this function s should be near a positive integer"))
    end
    # this series assumes s is near a positive integer
    n = Int(round(real(s)))
    τ = s - n
    # if real(s) > 0
    #     min_iterations = ceil( real(s) )
    # end
    ℒ = log(complex(-μ))  # '\u2112'
    # total = μ^(n-1)*Q(n-1, τ, ℒ; n_terms=n_terms)/SpecialFunctions.gamma(n)
    total = exp(
        (n - 1) * log(μ) + log(Q(n - 1, τ, ℒ; n_terms = n_terms)) -
        SpecialFunctions.loggamma(n),
    )
    converged = false
    a = Inf
    a_2 = Inf
    tmp = 1
    k = 0
    while k <= max_iterations && ~converged
        if n - k != 1
            a_2 = a
            a = tmp * SpecialFunctions.zeta(s - k)
            total += a
        end
        tmp *= μ / (k + 1)
        if k > min_iterations &&
           abs(a) / abs(total) < 0.5 * accuracy &&
           abs(a_2) / abs(total) < 0.5 * accuracy &&
           abs(a_2) > abs(a)
            # abs( (μ/twoπ)^k )/abs(total) < 0.05*accuracy
            converged = true
        end
        k = k + 1
    end
    # get correct value along the branch
    if isreal(z) && real(z) >= 1
        total -= 2 * π * im * μ^(s - 1) / SpecialFunctions.gamma(s)
    end
    return (total, k, 3)
end


"""
$(SIGNATURES)

Calculate harmonic numbers (http://mathworld.wolfram.com/HarmonicNumber.html)

## Arguments
* `n::Integer`: non-negative index of the Harmonic number to calculate

## Examples
```
julia> harmonic(2)
1.5
```
"""
function harmonic(n::Integer)
    if n < 0
        throw(DomainError(n))
    elseif n == 0
        return 0.0 # by convention (see Crandall, p.22)
    elseif n <= 10
        # perform exact sum for small n
        total = 0.0
        for k = 1:n
            total += 1.0 / k
        end
        return total
    else
        return Base.MathConstants.γ + SpecialFunctions.digamma(n + 1) # digamma(m) = ψ(m)
    end
end

"""
$(SIGNATURES)

Calculate harmonic numbers extended to non-integer arguments using the digamma form.

## Arguments
* `x::ComplexOrReal{Float64}`: index of the Harmonic number to calculate

## Examples
```
julia> harmonic(2.0)
1.5000000000000016
```
"""
harmonic(x::ComplexOrReal{Float64}) = Base.MathConstants.γ + SpecialFunctions.digamma(x + 1)

"""
$(SIGNATURES)

Calculate generalized harmonic numbers (e.g., see http://mathworld.wolfram.com/HarmonicNumber.html)

## Arguments
* `n::Integer`: non-negative index 1 of the Harmonic number to calculate
* `r::Real`: index 2 of the Harmonic number to calculate
It should be possible to extend this to complex r, but that requires more testing.

## Examples
```
julia> harmonic(2,1.5)
1.3535533905932737
```
"""
function harmonic(n::Integer, r::Real)
    if n < 0
        throw(DomainError(n))
    end
    if n == 0
        return 0.0
    end
    if r == 1
        return harmonic(n)
    end
    total = 0.0
    for k = 1:n
        total += 1.0 / k^r
    end
    return total
end


"""
$(SIGNATURES)

Calculates generalized harmonic numbers (e.g., see http://mathworld.wolfram.com/HarmonicNumber.html)

Using a better approach which works when both inputs are integers (https://carma.newcastle.edu.au/resources/jon/Preprints/Papers/Published-InPress/Oscillatory%20(Tapas%20II)/Papers/coffey-zeta.pdf, p.341)
 
## Arguments
* `n::Integer`: non-negative index 1 of the Harmonic number to calculate
* `r::Integer`: index 2 of the Harmonic number to calculate

## Examples
```
julia> harmonic(2,1)
1.5000000000000002
```
"""
function harmonic(n::Integer, r::Integer)
    if r < 1
        throw(DomainError(r))
    end
    return (-1)^(r - 1) * (
        SpecialFunctions.polygamma(r - 1, n + 1) - SpecialFunctions.polygamma(r - 1, 1)
    ) / SpecialFunctions.gamma(r)
end


const stieltjes_n = [
    Base.MathConstants.γ,
    -0.0728158454836767248605863758749013191377363383, # A082633
    -0.0096903631928723184845303860352125293590658061, # A086279
    +0.0020538344203033458661600465427533842857158044, # A086280
    +0.0023253700654673000574681701775260680009044694, # A086281
    +0.0007933238173010627017533348774444448307315394, # A086282
    -0.0002387693454301996098724218419080042777837151, # A183141
    -0.0005272895670577510460740975054788582819962534, # A183167
    -0.0003521233538030395096020521650012087417291805, # A183206
    -0.0000343947744180880481779146237982273906207895, # A184853
    +0.0002053328149090647946837222892370653029598537, # A184854
]
"""
$(SIGNATURES)

Provide the first 10 Stieltjes (generalized Euler-Mascheroni) constants (see
Abramowitz and Stegunm, 23.2.5) or (https://en.wikipedia.org/wiki/Stieltjes_constants).

## Arguments
* `n::Integer`: the number of elements to compute.

## Examples
```
julia> stieltjes(0)
0.5772156649015329
```
"""
function stieltjes(n::Integer)
    if n < 0
        throw(DomainError(n))
    elseif n > 10
        throw(DomainError(n, "Only the first 11 Stieltjes numbers are defined so far."))
    end

    return stieltjes_n[n+1]
end
