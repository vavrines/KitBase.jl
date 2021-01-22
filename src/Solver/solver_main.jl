"""
    solve!(
        KS::X,
        ctr::Y,
        face::Z,
        simTime,
    ) where {
        X<:AbstractSolverSet,
        Y<:AbstractArray{<:AbstractControlVolume,1},
        Z<:AbstractArray{<:AbstractInterface1D,1},
    }

    solve!(
        KS::X,
        ctr::Y,
        a1face::Z,
        a2face::Z,
        simTime,
    ) where {
        X<:AbstractSolverSet,
        Y<:AbstractArray{<:AbstractControlVolume,2},
        Z<:AbstractArray{<:AbstractInterface2D,2},
    }

Solution algorithm
- pre-process
- timestep calculation
- reconstruction
- evolution
- update

@args: solver setup
@args: array of control volumes
@args: array of interfaces
@args & return: time instant

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
    res = zero(KS.ib.wL)

    #--- main loop ---#
    #while true
    @showprogress for iter = 1:nt

        #dt = timestep(KS, ctr, simTime)
        reconstruct!(KS, ctr)
        evolve!(KS, ctr, face, dt; mode = Symbol(KS.set.flux), bc = Symbol(KS.set.boundary))
        update!(KS, ctr, face, dt, res; coll = Symbol(KS.set.collision), bc = Symbol(KS.set.boundary))

        #iter += 1
        t += dt

        if iter % 100 == 0
            println("iter: $(iter), time: $(simTime), dt: $(dt), res: $(res)")

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

end

function solve!(
    KS::X,
    ctr::Y,
    a1face::Z,
    a2face::Z,
    simTime,
) where {
    X<:AbstractSolverSet,
    Y<:AbstractArray{<:AbstractControlVolume,2},
    Z<:AbstractArray{<:AbstractInterface2D,2},
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
    @showprogress for iter = 1:nt
        reconstruct!(KS, ctr)
        evolve!(KS, ctr, a1face, a2face, dt; mode = Symbol(KS.set.flux), bc = Symbol(KS.set.boundary))
        update!(KS, ctr, a1face, a2face, dt, res; coll = Symbol(KS.set.collision), bc = Symbol(KS.set.boundary))

        t += dt

        if iter % 100 == 0
            println("iter: $(iter), time: $(simTime), dt: $(dt), res: $(res[1:end])")
        end

        if t > KS.set.maxTime || maximum(res) < 5.e-7
            break
        end
    end

    write_jld(KS, ctr, simTime)
    return t

end


"""
    timestep(KS, ctr, simTime)

Calculate timestep

@return: Δt

"""
function timestep(
    KS::X,
    ctr::Y,
    simTime,
) where {X<:AbstractSolverSet,Y<:AbstractArray{<:AbstractControlVolume,1}}

    tmax = 0.0

    if ctr[1].w isa Number
        
        @inbounds Threads.@threads for i = 1:KS.pSpace.nx
            prim = ctr[i].prim
            vmax = abs(ctr[i].prim[2])
            tmax = max(tmax, vmax / ctr[i].dx)
        end

    elseif KS.set.nSpecies == 1

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

function timestep(
    KS::X,
    ctr::Y,
    simTime,
) where {X<:AbstractSolverSet,Y<:AbstractArray{<:AbstractControlVolume,2}}

    tmax = 0.0

    if KS.set.nSpecies == 1

        @inbounds Threads.@threads for j = 1:KS.pSpace.ny
            for i = 1:KS.pSpace.nx
                prim = ctr[i, j].prim
                sos = sound_speed(prim, KS.gas.γ)
                umax = max(KS.vSpace.u1, abs(prim[2])) + sos
                vmax = max(KS.vSpace.v1, abs(prim[3])) + sos
                tmax = max(tmax, umax / ctr[i, j].dx + vmax / ctr[i, j].dy)
            end
        end

    elseif KS.set.nSpecies == 2

        @inbounds Threads.@threads for j = 1:KS.pSpace.ny
            for i = 1:KS.pSpace.nx
                prim = ctr[i, j].prim
                sos = sound_speed(prim, KS.gas.γ)
                umax = max(maximum(KS.vSpace.u1), maximum(abs.(prim[2, :]))) + sos
                vmax = max(maximum(KS.vSpace.v1), maximum(abs.(prim[3, :]))) + sos

                if KS.set.space[3:4] in ["3f", "4f"]
                    tmax = max(
                        tmax,
                        umax / ctr[i, j].dx + vmax / ctr[i, j].dy,
                        KS.gas.sol / ctr[i, j].dx + KS.gas.sol / ctr[i, j].dy
                    )
                else
                    tmax = max(tmax, umax / ctr[i, j].dx + vmax / ctr[i, j].dy)
                end
            end
        end

    end

    dt = KS.set.cfl / tmax
    dt = ifelse(dt < (KS.set.maxTime - simTime), dt, KS.set.maxTime - simTime)

    return dt

end