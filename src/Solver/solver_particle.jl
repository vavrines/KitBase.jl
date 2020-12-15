# ============================================================
# Update Algorithm for Particle Simulations
# ============================================================

function update!(
    KS::T,
    ctr::AbstractArray{ControlVolumeParticle1D,1},
    ptc::AbstractArray{Particle1D,1},
    ptc_temp::AbstractArray{Particle1D,1},
    face::AbstractArray{Interface1D,1},
    dt,
    res;
    coll = :bgk::Symbol,
    bc = :fix::Symbol,
) where {T<:AbstractSolverSet}

    res .= 0.

    update_transport!(KS, ctr, ptc, ptc_temp, dt)
    update_collision!(KS, ctr, ptc_temp, face, res, coll)
    update_boundary!(KS, ctr, ptc, ptc_temp, face, dt, coll, bc)

    ptc = deepcopy(ptc_temp)

    return nothing

end


"""
Update algorithm for particle transports

    update_transport!(KS::AbstractSolverSet, ctr::ControlVolumeParticle1D, ptc::Particle1D, ptc_temp::Particle1D, dt)

"""
function update_transport!(
    KS::T,
    ctr::AbstractArray{ControlVolumeParticle1D,1},
    ptc::AbstractArray{Particle1D,1},
    ptc_temp::AbstractArray{Particle1D,1},
    dt,
) where {T<:AbstractSolverSet}

    @inbounds Threads.@threads for i = 1:KS.pSpace.nx
        ctr[i].wg .= ctr[i].w
    end

    no_particle_temp = 0
    np = KS.gas.np
    @inbounds @simd for i in 1:np
        cell_no = ptc[i].idx
        tc = -ctr[cell_no].τ * log(rand())

        if tc > dt # class 1: free transport particles
            ptc[i].x += dt * ptc[i].v[1]

            ctr[cell_no].wg[1] -= ptc[i].m / ctr[cell_no].dx
            ctr[cell_no].wg[2] -= ptc[i].m * ptc[i].v[1] / ctr[cell_no].dx
            ctr[cell_no].wg[3] -= 0.5 * ptc[i].m * sum(ptc[i].v .^ 2) / ctr[cell_no].dx

            if ptc[i].x >= KS.pSpace.x0 && ptc[i].x <= KS.pSpace.x1
                no_particle_temp += 1

                cell_no_temp = Int(ceil((ptc[i].x - KS.pSpace.x0) / ctr[1].dx))
                #cell_no_temp = argmin(abs.(KS.pSpace.x[1:KS.pSpace.nx] .- ptc[i].x))
                #cell_no_temp = find_idx(KS.pSpace.x[1:KS.pSpace.nx], ptc[i].x, mode=:uniform)

                tb = boundary_time(ptc[i].x, ptc[i].v, ctr[cell_no_temp].x, ctr[cell_no_temp].dx)
                sample_particle!(ptc_temp[no_particle_temp], ptc[i].m, ptc[i].x, ptc[i].v, cell_no_temp, tb)
            end
       
        elseif tc > ptc[i].tb  # class 2: transfering collision particles

            ptc[i].x += dt * ptc[i].v[1]

            ctr[cell_no].wg[1] -= ptc[i].m / ctr[cell_no].dx
            ctr[cell_no].wg[2] -= ptc[i].m * ptc[i].v[1] / ctr[cell_no].dx
            ctr[cell_no].wg[3] -= 0.5 * ptc[i].m * sum(ptc[i].v .^ 2) / ctr[cell_no].dx

            if ptc[i].x >= KS.pSpace.x0 && ptc[i].x <= KS.pSpace.x1
                cell_no_temp = Int(ceil((ptc[i].x - KS.pSpace.x0) / ctr[1].dx)) 
                #cell_no_temp = find_idx(KS.pSpace.x[1:KS.pSpace.nx], ptc[i].x, mode=:uniform)

                ctr[cell_no_temp].wg[1] += ptc[i].m / ctr[cell_no_temp].dx
                ctr[cell_no_temp].wg[2] += ptc[i].m * ptc[i].v[1] / ctr[cell_no_temp].dx
                ctr[cell_no_temp].wg[3] +=
                    0.5 * ptc[i].m * sum(ptc[i].v .^ 2) / ctr[cell_no_temp].dx
            end
        end
    end

    @inbounds Threads.@threads for i = 1:KS.pSpace.nx
        ctr[i].w .= 0.0
    end
    @inbounds Threads.@threads for i = 1:no_particle_temp
        particle_cell_temp = ptc_temp[i].idx

        ctr[particle_cell_temp].w[1] += ptc_temp[i].m / ctr[particle_cell_temp].dx
        ctr[particle_cell_temp].w[2] +=
            ptc_temp[i].m * ptc_temp[i].v[1] / ctr[particle_cell_temp].dx
        ctr[particle_cell_temp].w[3] +=
            0.5 * ptc_temp[i].m * sum(ptc_temp[i].v .^ 2) / ctr[particle_cell_temp].dx
    end

    KS.gas.np = no_particle_temp
    return nothing

end


"""
Update algorithm for particle collisions

    update_collision!(KS::AbstractSolverSet, ctr::ControlVolumeParticle1D, ptc::Particle1D, ptc_temp::Particle1D, dt)

"""
function update_collision!(
    KS::T,
    ctr::AbstractArray{ControlVolumeParticle1D,1},
    ptc_temp::AbstractArray{Particle1D,1},
    face::AbstractArray{Interface1D,1},
    res::AbstractArray{<:AbstractFloat,1},
    coll = :bgk::Symbol,
) where {T<:AbstractSolverSet}

    no_particle_temp = KS.gas.np
    @inbounds Threads.@threads for i = 1:KS.pSpace.nx
        wold = deepcopy(ctr[i].wg)

        # fluid variables
        @. ctr[i].wg += (face[i].fw - face[i+1].fw) / ctr[i].dx
        if ctr[i].wg[1] < 0.0
            ctr[i].wg .= zero(ctr[i].wg) .+ 1e-8
        else
            primG = conserve_prim(ctr[i].wg, KS.gas.γ)
            if primG[end] < 0.0
                primG[end] = eltype(primG)(1e8)
            end
        end

        @. res += (ctr[i].wg - wold)^2

        @. ctr[i].w += ctr[i].wg
        if ctr[i].w[1] < 0.0
            ctr[i].w .= zero(ctr[i].w) .+ 1e-8
        else
            ctr[i].prim = conserve_prim(ctr[i].w, KS.gas.γ)
            if ctr[i].prim[end] < 0.0
                ctr[i].prim[end] = eltype(ctr[i].prim)(1e8)
                ctr[i].w .= prim_conserve(ctr[i].prim, KS.gas.γ)
            end
            ctr[i].τ = vhs_collision_time(ctr[i].prim, KS.gas.μᵣ, KS.gas.ω)
        end

        # particles
        no_particle_cell_temp = Int(round(ctr[i].wg[1] * ctr[i].dx / KS.gas.m))
        for j in 1:no_particle_cell_temp
            no_particle_temp += 1
            sample_particle!(ptc_temp[no_particle_temp], KS, ctr[i], i)
        end
    end

    KS.gas.np = no_particle_temp
    return nothing

end


function update_boundary!(
    KS::T,
    ctr::AbstractArray{ControlVolumeParticle1D,1},
    ptc::AbstractArray{Particle1D,1},
    ptc_temp::AbstractArray{Particle1D,1},
    face::AbstractArray{Interface1D,1},
    dt,
    coll = :bgk::Symbol,
    bc = :fix::Symbol,
) where {T<:AbstractSolverSet}

    np = KS.gas.np

    if bc == :fix

        ng = 1 - first(eachindex(KS.pSpace.x)) # ghost cell

        # left boundary
        @inbounds for i = 1:ng
            idx = 1 - i
            npc = ceil(ctr[idx].w[1] * ctr[idx].dx / KS.gas.m) |> Int
            for j = 1:npc
                np += 1
                sample_particle!(ptc_temp[np], KS, ctr[idx], idx)
            end
        end

        # right boundary
        @inbounds for i = 1:ng
            idx = KS.pSpace.nx + i
            npc = ceil(ctr[idx].w[1] * ctr[idx].dx / KS.gas.m) |> Int
            for j = 1:npc
                np += 1
                sample_particle!(ptc_temp[np], KS, ctr[idx], idx)
            end
        end

    end

    KS.gas.np = np
    return nothing

end
