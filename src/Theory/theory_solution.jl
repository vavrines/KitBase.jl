"""
Exact solution of 1D Riemann problem
https://github.com/ryarazi/ExactRiemannProblemSolver
"""

#--- TypeDefine.jl ---#
"""
$(TYPEDEF)

Hydrodynamic states

## Fields

$(FIELDS)
"""
struct HydroStatus
    rho::Float64 #density
    u::Float64 #velocity
    p::Float64 #pressure
    gamma::Float64 #adiabatic index
    c::Float64 #speed of cound
    A::Float64 #Toro constant (4.8)
    B::Float64 #Toro constant (4.8)

    function HydroStatus(rho, u, p, gamma)
        c = sqrt(gamma * p / rho)
        A = 2.0 / (gamma + 1.0) / rho
        B = (gamma - 1.0) / (gamma + 1.0) * p
        new(rho, u, p, gamma, c, A, B)
    end
end

"""
$(TYPEDEF)

Simplified hydrodynamic states: only saves {density, pressure and velocity}

## Fields

$(FIELDS)
"""
struct HydroStatusSimple
    rho::Float64 #density
    u::Float64 #velocity
    p::Float64 #pressure
end
HydroStatusSimple(status::HydroStatus) = HydroStatusSimple(status.rho, status.u, status.p)

# different types to dispatch the first guess of p_star calculation
# see "PressureGuess.jl"
abstract type FirstGuessScheme end
struct TwoRarefaction <: FirstGuessScheme end
struct PrimitiveValue <: FirstGuessScheme end
struct TwoShock <: FirstGuessScheme end
struct MeanPressure <: FirstGuessScheme end

# the type of the wave that will be created on some side of the contact discontinuity
abstract type WaveType end
struct Rarefaction <: WaveType end
struct Shock <: WaveType end
struct Vaccum <: WaveType end #for vaccum generated in the middle

# the side in which the calculation is happening
abstract type Side end
struct LeftSide <: Side end
struct RightSide <: Side end

# in some equations there is factor difference between the left and right side calculations
side_factor(::LeftSide) = -1
side_factor(::RightSide) = 1

# choose the some variable based on the side given
choose_side(left, right, ::LeftSide) = left
choose_side(left, right, ::RightSide) = right


#--- PressureGuess.jl ---#
# (4.46) from Toro
function p_guess(left::HydroStatus, right::HydroStatus, ::TwoRarefaction)
    @assert left.gamma == right.gamma
    gamma = left.gamma
    gamma_power = (gamma - 1) / 2.0 / gamma
    delta_u = right.u - left.u
    return (
        (left.c + right.c - 0.5 * (gamma - 1) * delta_u) /
        (left.c / left.p^gamma_power + right.c / right.p^gamma_power)
    )^(1.0 / gamma_power)
end

# (4.47) from Toro
function p_guess(left::HydroStatus, right::HydroStatus, ::PrimitiveValue)
    delta_u = right.u - left.u
    mean_pressure = 0.5 * (left.p + right.p)
    mean_density = 0.5 * (left.rho + right.rho)
    mean_speed_of_sound = 0.5 * (left.c + right.c)
    return mean_pressure - delta_u * mean_density * mean_speed_of_sound
end

# (4.48) from Toro
g(p, status::HydroStatus) = sqrt(status.A / (p + status.B))
function p_guess(left::HydroStatus, right::HydroStatus, ::TwoShock)
    p_hat = p_guess(left, right, PrimitiveValue)
    g_left = g(p_hat, left)
    g_right = g(p_hat, right)
    delta_u = right.u - left.u
    return (g_right * left.p + g_left * right.p - delta_u) / (g_right + g_left)
end

# (4.49) from Toro
p_guess(left::HydroStatus, right::HydroStatus, ::MeanPressure) = 0.5 * (left.p + right.p)


"""
$(SIGNATURES)

Main function that samples the Riemann problem and returns the density, pressure and velocity in the space points specified by the array x.

## Arguments
* `x`: an array which saves all the points in space in which we want to sample the problem (@ATTENTION: the contact discontinuity begins at x=0)
* `t`: floating number which tells us in what time we sample to problem
* `left`, `right`: struct from type HydroStatus (see "TypeDefine.jl") to get the initial state of the system
* `guess_scheme`: type of guess which will be used to solve the p_star equation (see "PressureGuess.jl")
* `TOL`: tolerance for the convergance of the p_star finding algorithm
"""
function sample_riemann_solution(
    x,
    t::Float64,
    left::HydroStatus = HydroStatus(1.0, 0.0, 1.0, 7.0 / 5.0),
    right::HydroStatus = HydroStatus(0.125, 0.0, 0.1, 7.0 / 5.0),
    guess_scheme::T = PrimitiveValue(),
    TOL::Float64 = 1.e-6,
) where {T<:FirstGuessScheme}
    @assert left.gamma == right.gamma

    if left.rho == 0 #material is on the right, vaccum on the left
        return sample_riemann_side_vaccum(x, t, right, RightSide())
    elseif right.rho == 0 #material is on the left, vaccum on the right
        return sample_riemann_side_vaccum(x, t, left, LeftSide())
    else
        return sample_riemann_regular(x, t, left, right, guess_scheme, TOL)
    end

end

"""
Sample the riemann problem in the situation where both sides are not vaccum (but vaccum can be generated through the dynamics)
"""
function sample_riemann_regular(
    x,
    t::Float64,
    left::HydroStatus,
    right::HydroStatus,
    guess_scheme::T,
    TOL::Float64,
) where {T<:FirstGuessScheme}
    @assert left.rho != 0.0
    @assert right.rho != 0.0

    profile = similar(x, HydroStatusSimple) #we return the values of rho, u and p in each point x in space

    #the status "far" from the center = like the starting condition
    status_left = HydroStatusSimple(left)
    status_right = HydroStatusSimple(right)

    x_over_t = x ./ t

    if 2.0 / (left.gamma - 1.0) * (left.c + right.c) <= right.u - left.u #vaccum generation
        #if vaccum is generated - the density in the middle will be zero
        p_star = 0.0
        u_star = 0.0
        wave_type_left = Vaccum()
        wave_type_right = Vaccum()
    else
        #if vaccum is not generated we can calculate the star profile in the middle
        p_star = p_star_calc(left, right, guess_scheme, TOL)
        u_star = u_star_calc(p_star, left, right)
        #check what kind of waves are in every direction (shock or rarefaction)
        wave_type_left = p_star > left.p ? Shock() : Rarefaction()
        wave_type_right = p_star > right.p ? Shock() : Rarefaction()
    end

    rho_star_left = rho_star_calc(p_star, left, wave_type_left)
    status_left_star = HydroStatusSimple(rho_star_left, u_star, p_star) #the profile in the left near the contact discontinuity
    head_speed_left = head_speed_calc(p_star, u_star, left, LeftSide(), wave_type_left)
    tail_speed_left = tail_speed_calc(p_star, u_star, left, LeftSide(), wave_type_left)

    rho_star_right = rho_star_calc(p_star, right, wave_type_right)
    status_right_star = HydroStatusSimple(rho_star_right, u_star, p_star) #the profile in the right near the contact discontinuity
    head_speed_right = head_speed_calc(p_star, u_star, right, RightSide(), wave_type_right)
    tail_speed_right = tail_speed_calc(p_star, u_star, right, RightSide(), wave_type_right)


    for i = 1:length(x)
        S = x_over_t[i] #this is like the S which Toro use in Section 4.5

        #see Figure 4.14 in Toro for the flow of the following lines
        side = S < u_star ? LeftSide() : RightSide()
        status = choose_side(left, right, side)
        head_speed = choose_side(head_speed_left, head_speed_right, side)
        tail_speed = choose_side(tail_speed_left, tail_speed_right, side)

        #this is used to flip the direction of the inequallity in the branching between the left and the right side
        #the xor which will be in the following lines uses this boolean like some "controlled not"
        #when right_condition is "true" the inequallity will be flipped
        #when right_condition is "false" the inequallity will be stay the same
        right_condition = isa(side, RightSide)

        if xor(S < head_speed, right_condition)
            profile[i] = choose_side(status_left, status_right, side)
        elseif xor(S < tail_speed, right_condition) #can only happen in Rarefaction, because in Shock  head_speed == tail_speed
            profile[i] = rarefaction_profile(S, status, side)
        else
            profile[i] = choose_side(status_left_star, status_right_star, side)
        end
    end

    return profile
end

"""
Sample the riemann problem in the situation where the gas in "side" is actually a vaccum
"""
function sample_riemann_side_vaccum(
    x,
    t::Float64,
    status::HydroStatus,
    side::T,
) where {T<:Side}

    profile = similar(x, HydroStatusSimple) #we return the values of rho, u and p in each point x in space
    vaccum = HydroStatusSimple(0.0, 0.0, 0.0) #vaccum
    simple = HydroStatusSimple(status)

    x_over_t = x ./ t

    tail_speed = tail_speed_calc(0.0, 0.0, status, side, Vaccum()) #the head of the rarefaciton wave
    head_speed = head_speed_calc(0.0, 0.0, status, side, Vaccum()) #the tail of the rarefaciton wave (=the boundary with the vaccum)

    #see "sample_riemann_regular" for explnation about this boolean
    right_condition = isa(side, RightSide)

    for i = 1:length(x)
        S = x_over_t[i]

        if xor(S < head_speed, right_condition)
            profile[i] = simple
        elseif xor(S < tail_speed, right_condition)
            profile[i] = rarefaction_profile(S, status, side)
        else
            profile[i] = vaccum
        end
    end

    return profile
end

# this is the density, velocity and pressure profile of a rarefaction wave in some specific x/t	
# see (4.56) and (4.63) from Toro
function rarefaction_profile(x_over_t, status::HydroStatus, side::T) where {T<:Side}
    gamma = status.gamma
    gamma_ratio = (gamma - 1.0) / (gamma + 1.0)
    gamma_plus = 2.0 / (gamma + 1.0)
    gamma_minus = 2.0 / (gamma - 1.0)

    rarefaction_factor =
        (
            gamma_plus - side_factor(side) * gamma_ratio / status.c * (status.u - x_over_t)
        )^gamma_minus

    rho = status.rho * rarefaction_factor
    u =
        gamma_plus *
        (-side_factor(side) * status.c + (1.0 / gamma_minus) * status.u + x_over_t)
    p = status.p * rarefaction_factor^gamma

    return HydroStatusSimple(rho, u, p)
end

# function (4.6) and (4.7) from Toro
f(p, status::HydroStatus, ::Shock) = (p - status.p) * sqrt(status.A / (p + status.B))
f(p, status::HydroStatus, ::Rarefaction) =
    2.0 * status.c / (status.gamma - 1.0) *
    ((p / status.p)^((status.gamma - 1.0) / 2.0 / status.gamma) - 1.0)
function f(p, status::HydroStatus)
    wave_type = p > status.p ? Shock() : Rarefaction()
    return f(p, status, wave_type)
end

# solver for p_star based on the first guess of the pressure(see PressureGuess.jl)
# uses secant method to find solve the equation
function p_star_calc(
    left::HydroStatus,
    right::HydroStatus,
    guess_scheme::T,
    TOL::Float64,
) where {T<:FirstGuessScheme}
    p0 = max(TOL, p_guess(left, right, guess_scheme))

    #we define p = exp(u) so inside the interation it won't get negative value
    u0 = log(p0) #p = exp(u) => u = log(p)
    f_star(u) = f(exp(u), left) + f(exp(u), right) + right.u - left.u #function (4.5) from Toro
    u = find_zero(f_star, u0, Order1(), rtol = TOL) #Secant method
    return exp(u) #p = exp(u) 
end

# equation (4.9) from Toro
u_star_calc(p_star, left::HydroStatus, right::HydroStatus) =
    0.5 * (left.u + right.u) + 0.5 * (f(p_star, right) - f(p_star, left))


# equations (4.23) and (4.32) from Toro
rho_star_calc(p_star, status::HydroStatus, ::Rarefaction) =
    status.rho * (p_star / status.p)^(1.0 / status.gamma)

# equation (4.25) and (4.34) from Toro
c_star_calc(p_star, status::HydroStatus, ::Rarefaction) =
    status.c * (p_star / status.p)^((status.gamma - 1.0) / 2.0 / status.gamma)

# equations (4.50) and (4.57) from Toro
function rho_star_calc(p_star, status::HydroStatus, ::Shock)
    gamma_ratio = (status.gamma - 1.0) / (status.gamma + 1.0)
    pressure_ratio = p_star / status.p
    return status.rho * (pressure_ratio + gamma_ratio) /
           (pressure_ratio * gamma_ratio + 1.0)
end

# equation (4.55) and (4.62) from Toro 
head_speed_calc(
    p_star,
    u_star,
    status::HydroStatus,
    side::T,
    ::Rarefaction,
) where {T<:Side} = status.u + side_factor(side) * status.c

# equation (4.55) and (4.62) from Toro 
function tail_speed_calc(
    p_star,
    u_star,
    status::HydroStatus,
    side::T,
    ::Rarefaction,
) where {T<:Side}
    c_star = c_star_calc(p_star, status, Rarefaction())
    return u_star + side_factor(side) * c_star
end

# equations (4.52) and (4.59) from Toro
function head_speed_calc(
    p_star,
    u_star,
    status::HydroStatus,
    side::T,
    ::Shock,
) where {T<:Side}
    gamma = status.gamma
    gamma_ratio1 = (gamma + 1.0) / 2.0 / gamma
    gamma_ratio2 = (gamma - 1.0) / 2.0 / gamma
    return status.u +
           side_factor(side) *
           status.c *
           sqrt(gamma_ratio1 * p_star / status.p + gamma_ratio2)
end

# for shock the head and the tail are the same
tail_speed_calc(p_star, u_star, status::HydroStatus, side::T, ::Shock) where {T<:Side} =
    head_speed_calc(p_star, u_star, status, side, Shock())


# for vaccum generation the density in the middle is zero (= vaccum)
rho_star_calc(p_star, status::HydroStatus, ::Vaccum) = 0.0

# those are the velocities of the head and the tail of the rarefaciton wave when a vaccum is generated
# see (4.76), (4.77), (4.79) and (4.80) from Toro
head_speed_calc(p_star, u_star, status::HydroStatus, side::T, ::Vaccum) where {T<:Side} =
    status.u + side_factor(side) * status.c
tail_speed_calc(p_star, u_star, status::HydroStatus, side::T, ::Vaccum) where {T<:Side} =
    status.u - side_factor(side) * 2 * status.c / (status.gamma - 1.0)
