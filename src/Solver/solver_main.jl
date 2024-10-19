"""
$(SIGNATURES)

Solution algorithm for 1D structured and unstructured mesh

## Arguments
- `KS`: SolverSet
- `ctr`: vector of cell-centered solution
- `face`: vector of cell interface
- `simTime`: simulation time
"""
function solve!(
    KS::AbstractSolverSet,
    ctr::AV{<:AbstractControlVolume},
    face::AV{<:AbstractInterface},
    simTime;
    steady=false,
)

    #--- initial checkpoint ---#
    write_jld(KS, ctr, simTime)

    #--- setup ---#
    iter = 0
    t = deepcopy(simTime)
    dt = timestep(KS, ctr, simTime)
    nt = Int(floor(KS.set.maxTime / dt)) + 1
    res = zero(ctr[1].w)

    #--- main loop ---#
    #while true
    @showprogress for iter in 1:nt
        #dt = timestep(KS, ctr, simTime)
        reconstruct!(KS, ctr)
        evolve!(
            KS,
            ctr,
            face,
            dt;
            mode=symbolize(KS.set.flux),
            bc=symbolize(KS.set.boundary),
        )
        update!(
            KS,
            ctr,
            face,
            dt,
            res;
            coll=symbolize(KS.set.collision),
            bc=symbolize(KS.set.boundary),
        )

        #iter += 1
        t += dt

        if iter % 500 == 0
            println("iter: $(iter), time: $(t), dt: $(dt), res: $(res)")

            if iter % 1000 == 0
                write_jld(KS, ctr, iter)
            end
        end

        if t > KS.set.maxTime
            break
        end
        if steady == true
            if maximum(res) < 5.e-7
                break
            end
        end
    end

    write_jld(KS, ctr, simTime)
    return t
end

"""
$(SIGNATURES)

Solution algorithm for 2D structured mesh

## Arguments
- `KS`: SolverSet
- `ctr`: matrix of cell-centered solution
- `a1face`: matrix of cell interface perpendicular to `x` axis
- `a2face`: matrix of cell interface perpendicular to `y` axis
- `simTime`: simulation time
"""
function solve!(
    KS::AbstractSolverSet,
    ctr::AM{<:AbstractControlVolume},
    a1face::T,
    a2face::T,
    simTime;
    steady=false,
) where {T<:AA{<:AbstractInterface,2}}

    #--- initial checkpoint ---#
    write_jld(KS, ctr, simTime)

    #--- setup ---#
    iter = 0
    t = deepcopy(simTime)
    dt = timestep(KS, ctr, simTime)
    nt = Int(floor(KS.set.maxTime / dt)) + 1
    res = zero(ctr[1].w)

    #--- main loop ---#
    @showprogress for iter in 1:nt
        reconstruct!(KS, ctr)
        evolve!(
            KS,
            ctr,
            a1face,
            a2face,
            dt;
            mode=symbolize(KS.set.flux),
            bc=symbolize(KS.set.boundary),
        )
        update!(
            KS,
            ctr,
            a1face,
            a2face,
            dt,
            res;
            coll=symbolize(KS.set.collision),
            bc=symbolize(KS.set.boundary),
        )

        t += dt

        if iter % 500 == 0
            println("iter: $(iter), time: $(t), dt: $(dt), res: $(res[1:end])")
        end

        if t > KS.set.maxTime
            break
        end
        if steady == true
            if maximum(res) < 5.e-7
                break
            end
        end
    end

    write_jld(KS, ctr, simTime)
    return t
end

"""
$(SIGNATURES)

Calculate timestep based on the current solution

## Arguments
- `KS`: SolverSet
- `ctr`: array of cell-centered solution
- `simTime`: simulation time
"""
function timestep(KS::AbstractSolverSet, ctr::AV{<:AbstractControlVolume}, simTime=0.0)
    tmax = 0.0

    if ctr[1].w isa Number
        @inbounds @threads for i in 1:KS.ps.nx
            prim = ctr[i].prim
            vmax = abs(ctr[i].prim[2])
            tmax = max(tmax, vmax / KS.ps.dx[i])
        end

    elseif KS.set.nSpecies == 1
        @inbounds @threads for i in 1:KS.ps.nx
            prim = ctr[i].prim
            sos = sound_speed(prim, KS.gas.γ)
            vmax = begin
                if KS.vs isa Nothing
                    abs(prim[2]) + sos
                else
                    max(KS.vs.u1, abs(prim[2])) + sos
                end
            end
            tmax = max(tmax, vmax / KS.ps.dx[i])
        end

    elseif KS.set.nSpecies == 2
        @inbounds @threads for i in 1:KS.ps.nx
            prim = ctr[i].prim
            sos = sound_speed(prim, KS.gas.γ)
            vmax = begin
                if KS.vs isa Nothing
                    maximum(abs.(prim[2, :])) + sos
                else
                    max(maximum(KS.vs.u1), maximum(abs.(prim[2, :]))) + sos
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

function timestep(KS::AbstractSolverSet, ctr::AM{<:AbstractControlVolume}, simTime=0.0)
    nx, ny, dx, dy = begin
        if KS.ps isa CSpace2D
            KS.ps.nr, KS.ps.nθ, KS.ps.dr, KS.ps.darc
        else
            KS.ps.nx, KS.ps.ny, KS.ps.dx, KS.ps.dy
        end
    end

    tmax = 0.0

    if KS.set.nSpecies == 1
        @inbounds @threads for j in 1:ny
            for i in 1:nx
                prim = ctr[i, j].prim
                sos = sound_speed(prim, KS.gas.γ)
                umax, vmax = begin
                    if KS.vs isa Nothing
                        abs(prim[2]) + sos, abs(prim[3]) + sos
                    else
                        max(KS.vs.u1, abs(prim[2])) + sos, max(KS.vs.v1, abs(prim[3])) + sos
                    end
                end
                tmax = max(tmax, umax / dx[i, j] + vmax / dy[i, j])
            end
        end

    elseif KS.set.nSpecies == 2
        @inbounds @threads for j in 1:ny
            for i in 1:nx
                prim = ctr[i, j].prim
                sos = sound_speed(prim, KS.gas.γ)
                umax = max(maximum(KS.vs.u1), maximum(abs.(prim[2, :]))) + sos
                vmax = max(maximum(KS.vs.v1), maximum(abs.(prim[3, :]))) + sos

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
    KS::AbstractSolverSet,
    ctr::AV{<:AbstractUnstructControlVolume},
    simTime=0.0,
)
    tmax = 0.0

    if KS.set.nSpecies == 1
        @inbounds @threads for i in eachindex(ctr)
            prim = ctr[i].prim
            sos = sound_speed(prim, KS.gas.γ)
            umax = KS.vs.u1 + sos
            vmax = KS.vs.v1 + sos
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
