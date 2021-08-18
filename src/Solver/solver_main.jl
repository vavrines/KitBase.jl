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

- @args: solver setup
- @args: array of control volumes
- @args: array of interfaces
- @args & return: time instant

"""
function solve!(
    KS::X,
    ctr::Y,
    face::Z,
    simTime,
) where {
    X<:AbstractSolverSet,
    Y<:AbstractArray{<:Union{ControlVolume1D,ControlVolume1D1F,ControlVolume1D2F},1},
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
        update!(
            KS,
            ctr,
            face,
            dt,
            res;
            coll = Symbol(KS.set.collision),
            bc = Symbol(KS.set.boundary),
        )

        #iter += 1
        t += dt

        if iter % 500 == 0
            println("iter: $(iter), time: $(t), dt: $(dt), res: $(res)")

            #if iter%1000 == 0
            #    write_jld(KS, ctr, iter)
            #end
        end

        if t > KS.set.maxTime || maximum(res) < 5.e-7
            break
        end

    end

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
        evolve!(
            KS,
            ctr,
            a1face,
            a2face,
            dt;
            mode = Symbol(KS.set.flux),
            bc = Symbol(KS.set.boundary),
        )
        update!(
            KS,
            ctr,
            a1face,
            a2face,
            dt,
            res;
            coll = Symbol(KS.set.collision),
            bc = Symbol(KS.set.boundary),
        )

        t += dt

        if iter % 500 == 0
            println("iter: $(iter), time: $(t), dt: $(dt), res: $(res[1:end])")
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

Calculate timestep based on current solutions

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
            tmax = max(tmax, vmax / KS.ps.dx[i])
        end

    elseif KS.set.nSpecies == 1

        @inbounds Threads.@threads for i = 1:KS.pSpace.nx
            prim = ctr[i].prim
            sos = sound_speed(prim, KS.gas.γ)
            vmax = begin
                if KS.vs isa Nothing
                    abs(prim[2]) + sos
                else
                    max(KS.vSpace.u1, abs(prim[2])) + sos
                end
            end
            tmax = max(tmax, vmax / KS.ps.dx[i])
        end

    elseif KS.set.nSpecies == 2

        @inbounds Threads.@threads for i = 1:KS.pSpace.nx
            prim = ctr[i].prim
            sos = sound_speed(prim, KS.gas.γ)
            vmax = begin
                if KS.vs isa Nothing
                    maximum(abs.(prim[2, :])) + sos
                else
                    max(maximum(KS.vSpace.u1), maximum(abs.(prim[2, :]))) + sos
                end
            end

            if KS.set.space[3:4] in ["3f", "4f"]
                tmax = max(tmax, vmax / KS.ps.dx[i], KS.gas.sol / KS.ps.dx[i])
            else
                tmax = max(tmax, vmax / KS.ps.dx[i])
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

    nx, ny, dx, dy = begin
        if KS.ps isa CSpace2D
            KS.ps.nr, KS.ps.nθ, KS.ps.dr, KS.ps.darc
        else
            KS.ps.nx, KS.ps.ny, KS.ps.dx, KS.ps.dy
        end
    end

    tmax = 0.0

    if KS.set.nSpecies == 1

        @inbounds Threads.@threads for j = 1:ny
            for i = 1:nx
                prim = ctr[i, j].prim
                sos = sound_speed(prim, KS.gas.γ)
                umax = ifelse(
                    KS.vs isa Nothing,
                    abs(prim[2]) + sos,
                    max(KS.vSpace.u1, abs(prim[2])) + sos,
                )
                vmax = ifelse(
                    KS.vs isa Nothing,
                    abs(prim[3]) + sos,
                    max(KS.vSpace.v1, abs(prim[3])) + sos,
                )
                tmax = max(tmax, umax / dx[i, j] + vmax / dy[i, j])
            end
        end

    elseif KS.set.nSpecies == 2

        @inbounds Threads.@threads for j = 1:ny
            for i = 1:nx
                prim = ctr[i, j].prim
                sos = sound_speed(prim, KS.gas.γ)
                umax = max(maximum(KS.vSpace.u1), maximum(abs.(prim[2, :]))) + sos
                vmax = max(maximum(KS.vSpace.v1), maximum(abs.(prim[3, :]))) + sos

                if KS.set.space[3:4] in ["3f", "4f"]
                    tmax = max(
                        tmax,
                        umax / dx[i, j] + vmax / dy[i, j],
                        KS.gas.sol / dx[i, j] + KS.gas.sol / dy[i, j],
                    )
                else
                    tmax = max(tmax, umax / dx[i, j] + vmax / dy[i, j])
                end
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
) where {X<:AbstractSolverSet,Y<:AbstractVector{<:AbstractUnstructControlVolume}}

    tmax = 0.0

    if KS.set.nSpecies == 1

        @inbounds Threads.@threads for i in eachindex(ctr)
            prim = ctr[i].prim
            sos = sound_speed(prim, KS.gas.γ)
            umax = KS.vSpace.u1 + sos
            vmax = KS.vSpace.v1 + sos
            pd =
                abs.([dot([umax, vmax], ctr[i].n[j]) for j in eachindex(ctr[i].n)]) ./ ctr[i].dx
            #tmax = max(tmax, maximum(pd))
            tmax = max(tmax, sqrt(umax^2 + vmax^2) / minimum(ctr[i].dx))
        end

    elseif KS.set.nSpecies == 2

    end

    dt = KS.set.cfl / tmax
    dt = ifelse(dt < (KS.set.maxTime - simTime), dt, KS.set.maxTime - simTime)

    return dt

end
