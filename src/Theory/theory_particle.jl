"""
    sample_particle!(ptc::Particle1D, m, x, v, idx, tb)
    sample_particle!(ptc::Particle1D, KS, ctr, idx, t0 = 1.0)

Sample particles from local flow conditions

"""
function sample_particle!(ptc::Particle1D, m, x, v, idx, tb)
    ptc.m = m
    ptc.x = x
    ptc.v = v
    ptc.idx = idx
    ptc.tb = tb
end

function sample_particle!(ptc::Particle1D, KS, ctr, idx, t0 = 1.0)

    ptc.m = KS.gas.m
    ptc.v .= sample_velocity(ctr.prim)
    ptc.x = ctr.x + (rand() - 0.5) * ctr.dx
    ptc.idx = idx
    ptc.tb = boundary_time(ptc.x, ptc.v, ctr.x, ctr.dx, t0)

end


"""
    sample_velocity(prim::AbstractArray{<:Real,1})

Sample particle velocity from local flow conditions

"""
function sample_velocity(prim::T) where {T<:AbstractArray{<:Real,1}}
    rd1 = rand(3)
    rd2 = rand(3)

    return @. sqrt(-log(rd1) / prim[end]) * sin(2.0 * π * rd2) + [prim[2], 0., 0.]
end


"""
    boundary_time(x, v, x0, dx0, t0 = 1.0)

Calculate particle flight time to cell boundary

"""
function boundary_time(x, v, x0, dx0, t0 = 1.0)
    if v[1] < 0.
        tb = (x0 - 0.5 * dx0 - x) / v[1]
    elseif v[1] > 0.
        tb = (x0 + 0.5 * dx0 - x) / v[1]
    else
        tb = t0
    end

    return tb
end