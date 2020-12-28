"""
Solution algorithm

* 1D solver: `solve!(KS::SolverSet, ctr::AbstractArray{<:AbstractControlVolume1D,1},
    face::Array{<:AbstractInterface1D,1}, simTime::Float64)`

* @return: ending time

"""
function solve!(
    KS::X,
    ctr::Y,
    face::Z,
    simTime,
) where {
    X<:AbstractSolverSet,
    Y<:AbstractArray{<:AbstractControlVolume,1},
    Z<:AbstractArray{<:AbstractInterface1D,1},
}

    #--- initial checkpoint ---#
    write_jld(KS, ctr, simTime)

    #--- setup ---#
    iter = 0
    t = deepcopy(simTime)
    dt = timestep(KS, ctr, simTime)
    nt = Int(floor(KS.set.maxTime / dt)) + 1
    res = zeros(axes(KS.ib.wL))

    #--- main loop ---#
    #while true
    @showprogress for iter = 1:nt

        #dt = timestep(KS, ctr, simTime)
        reconstruct!(KS, ctr)
        evolve!(KS, ctr, face, dt)
        update!(KS, ctr, face, dt, res)

        #iter += 1
        t += dt

        if iter % 100 == 0
            println("iter: $(iter), time: $(simTime), dt: $(dt), res: $(res[1:end])")

            #if iter%1000 == 0
            #    write_jld(KS, ctr, iter)
            #end
        end

        if t > KS.set.maxTime || maximum(res) < 5.e-7
            break
        end

    end # loop

    write_jld(KS, ctr, simTime)
    return t

end # function


"""
Timestep calculator

* 1D solver: `timestep(KS::SolverSet, ctr::AbstractArray{<:AbstractControlVolume1D,1},
    simTime::Real)`

* @return: Δt

"""
function timestep(
    KS::X,
    ctr::Y,
    simTime,
) where {X<:AbstractSolverSet,Y<:AbstractArray{<:AbstractControlVolume,1}}

    tmax = 0.0

    if KS.set.nSpecies == 1

        @inbounds Threads.@threads for i = 1:KS.pSpace.nx
            prim = ctr[i].prim
            sos = sound_speed(prim, KS.gas.γ)
            vmax = max(KS.vSpace.u1, abs(prim[2])) + sos
            tmax = max(tmax, vmax / ctr[i].dx)
        end

    elseif KS.set.nSpecies == 2

        @inbounds Threads.@threads for i = 1:KS.pSpace.nx
            prim = ctr[i].prim
            sos = sound_speed(prim, KS.gas.γ)
            vmax = max(maximum(KS.vSpace.u1), maximum(abs.(prim[2, :]))) + sos

            if KS.set.space[3:4] in ["3f", "4f"]
                tmax = max(tmax, vmax / ctr[i].dx, KS.gas.sol / ctr[i].dx)
            else
                tmax = max(tmax, vmax / ctr[i].dx)
            end
        end

    end

    dt = KS.set.cfl / tmax
    dt = ifelse(dt < (KS.set.maxTime - simTime), dt, KS.set.maxTime - simTime)

    return dt

end
