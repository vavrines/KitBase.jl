# ============================================================
# Update Algorithm for Particle Simulations
# ============================================================

function update!(
    KS::T,
    ctr::ControlVolumeParticle1D,
    ptc::Particle1D,
    ptc_temp::Particle1D,
    face::Interface1D,
    dt,
    residual;
    coll = :bgk::Symbol,
    bc = :fix::Symbol,
) where {T<:AbstractSolverSet}

    update_transport!(KS, ctr, ptc, ptc_temp, dt)
    update_collision!(KS, ctr, ptc, ptc_temp, face, dt, residual; coll = coll)
    update_boundary!(KS, ctr, ptc, ptc_temp, face, dt; coll = coll, bc = bc)

    no_particle = no_particle_temp
    particle = particle_temp

    return nothing

end


"""
Update algorithm for particle transports

    update_transport!(KS::AbstractSolverSet, ctr::ControlVolumeParticle1D, ptc::Particle1D, ptc_temp::Particle1D, dt)

"""
function update_transport!(
    KS::T,
    ctr::ControlVolumeParticle1D,
    ptc::Particle1D,
    ptc_temp::Particle1D,
    dt,
) where {T<:AbstractSolverSet}

    for i = 1:KS.pSpace.nx
        ctr[i].wg .= ctr[i].w
    end

    np_tmp = 0


    no_particle_temp = 0
    np = KS.gas.np
    for i = 1:np
        cell_no = ptc[i].idx

        tc = -ctr[cell_no].τ * log(rand())

        if tc > dt
            ptc[i].x += dt * ptc[i].v[1]

            cell_no_temp = ptc[i].idx

            ctr[cell_no].wg[1] -= ptc[i].m / ctr[cell_no].dx
            ctr[cell_no].wg[2] -= ptc[i].m * ptc[i].v[1] / ctr[cell_no].dx
            ctr[cell_no].wg[3] -= 0.5 * ptc[i].m * sum(ptc[i].v .^ 2) / ctr[cell_no].dx

            if ptc[i].x > KS.pSpace.x0 && ptc[i].x < KS.pSpace.x1
                no_particle_temp += 1

                ptc_temp[no_particle_temp].m = ptc[i].m
                ptc_temp[no_particle_temp].x = ptc[i].x
                ptc_temp[no_particle_temp].v .= ptc[i].v

                cell_no_temp = Int(ceil((ptc[i].x - KS.pSpace.x0) / ctr[cell_no_temp].dx))
                #cell_no_temp = argmin(abs.(KS.pSpace.x[1:KS.pSpace.nx] .- ptc[i].x))

                ptc_temp[no_particle_temp].idx = cell_no_temp

                if ptc_temp[no_particle_temp].v[1] < 0.0
                    ptc_temp[no_particle_temp].tb =
                        (
                            ctr[cell_no_temp].x - 0.5 * ctr[cell_no_temp].dx -
                            ptc_temp[no_particle_temp].x
                        ) / ptc_temp[no_particle_temp].v[1]
                elseif ptc_temp[no_particle_temp].v[1] > 0.0
                    ptc_temp[no_particle_temp].tb =
                        (
                            ctr[cell_no_temp].x + 0.5 * ctr[cell_no_temp].dx -
                            ptc_temp[no_particle_temp].x
                        ) / ptc_temp[no_particle_temp].v[1]
                else
                    ptc_temp[no_particle_temp].tb = 1.0
                end
            end

        elseif tc > ptc[i].tb

            ptc[i].x += dt * ptc[i].v[1]

            cell_no_temp = ptc[i].idx

            ctr[cell_no].wg[1] -= ptc[i].m / dx
            ctr[cell_no].wg[2] -= ptc[i].m * ptc[i].v[1] / ctr[cell_no].dx
            ctr[cell_no].wg[3] -= 0.5 * ptc[i].m * sum(ptc[i].v .^ 2) / ctr[cell_no].dx

            if ptc[i].x > KS.pSpace.x0 && ptc[i].x < KS.pSpace.x1
                cell_no_temp = Int(ceil((ptc[i].x - KS.pSpace.x0) / ctr[cell_no_temp].dx))
                ctr[cell_no_temp].wg[1] += ptc[i].m / dx
                ctr[cell_no_temp].wg[2] += ptc[i].m * ptc[i].v[1] / ctr[cell_no_temp].dx
                ctr[cell_no_temp].wg[3] +=
                    0.5 * ptc[i].m * sum(ptc[i].v .^ 2) / ctr[cell_no_temp].dx
            end
        end
    end

    for i = 1:KS.pSpace.nx
        ctr[i].w .= 0.0
    end

    for i = 1:no_particle_temp
        particle_cell_temp = ptc_temp[i].idx

        ctr[particle_cell_temp].w[1] += ptc_temp[i].m / ctr[particle_cell_temp].dx
        ctr[particle_cell_temp].w[2] +=
            ptc_temp[i].m * ptc_temp[i].v[1] / ctr[particle_cell_temp].dx
        ctr[particle_cell_temp].w[3] +=
            0.5 * ptc_temp[i].m * sum(ptc_temp[i].v .^ 2) / ctr[particle_cell_temp].dx
    end

    return nothing

end


"""
Update algorithm for particle collisions

    update_collision!(KS::AbstractSolverSet, ctr::ControlVolumeParticle1D, ptc::Particle1D, ptc_temp::Particle1D, dt)

"""
function update_collision!(
    KS::T,
    ctr::ControlVolumeParticle1D,
    ptc::Particle1D,
    ptc_temp::Particle1D,
    face::Interface1D,
    dt,
    residual;
    coll = :bgk::Symbol,
) where {T<:AbstractSolverSet}

    for i = 1:KS.pSpace.nx
        @. ctr[i].wg += (face[i].fw - face[i+1].fw) / ctr[i].dx

        if ctr[i].wg[1] < 0.0
            ctr[i].wg .= zero(ctr[i].wg)
        else
            ctr[i].prim = conserve_prim(ctr[i].wg, KS.gas.γ)
            if ctr[i].prim[end] < 0.0
                ctr[i].prim[end] = eltype(ctr[i].prim)(1e8)
            end
        end

        @. ctr[i].w += ctr[i].wg

        prim = conserve_prim(ctr[i].w, KS.gas.γ)
        if ctr[i].w[1] < 0.0
            ctr[i].w .= zero(ctr[i].w)
        else
            prim = conserve_prim(ctr[i].w, KS.gas.γ)
            if ctr[i].prim[end] < 0.0
                ctr[i].prim[end] = eltype(ctr[i].prim)(1e8)
                ctr[i].w .= prim_conserve(prim, KS.gas.γ)
            end
        end
        ctr[i].τ = vhs_collision_time(prim, KS.gas.μᵣ, KS.gas.ω)
    end


    for i = 1:KS.pSpace.nx
        no_particle_cell_temp = Int(round(ctr[i].wg[1] * ctr[i].dx / KS.gas.m))
        for j = 1:no_particle_cell_temp
            no_particle_temp = no_particle_temp + 1

            rd1 = rand(3)
            rd2 = rand(3)
            rd = rand()

            ptc_temp[no_particle_temp].mass = KS.gas.m
            @. ptc_temp[no_particle_temp].v =
                sqrt(-log(rd1) / ctr[i].prim[end]) * sin(2.0 * π * rd2)
            ptc_temp[no_particle_temp].v[1] += ctr[i].prim[2]
            ptc_temp[no_particle_temp].x = ctr[i].x + (rd - 0.5) * ctr[i].dx
            ptc_temp[no_particle_temp].idx = i

            if ptc_temp[no_particle_temp].v[1] < 0
                ptc_temp[no_particle_temp].tb =
                    (ctr[i].x - 0.5 * ctr[i].dx - ptc_temp[no_particle_temp].x) /
                    ptc_temp[no_particle_temp].v[1]
            elseif ptc_temp[no_particle_temp].v[1] > 0
                ptc_temp[no_particle_temp].tb =
                    (ctr[i].x + 0.5 * ctr[i].dx - ptc_temp[no_particle_temp].x) /
                    ptc_temp[no_particle_temp].v[1]
            else
                ptc_temp[no_particle_temp].tb = 1.0
            end
        end
    end

    return nothing

end

function update_boundary!(
    KS::T,
    ctr::ControlVolumeParticle1D,
    ptc::Particle1D,
    ptc_temp::Particle1D,
    face::Interface1D,
    dt,
    coll = :bgk::Symbol,
    bc = :fix::Symbol,
) where {T<:AbstractSolverSet}

    ng = 1 - first(eachindex(KS.pSpace.x))

    if bc == :fix

        for i = 1:ng
            idx = 1 - i

            no_particle_cell_temp = round(ctr[idx].w[idx] * ctr[idx].dx / KS.gas.m) |> Int
            for j = 1:no_particle_cell_temp
                no_particle_temp = no_particle_temp + 1

                rd1 = rand(3)
                rd2 = rand(3)
                rd = rand()

                ptc_temp[no_particle_temp].mass = KS.gas.m
                @. ptc_temp[no_particle_temp].v =
                    sqrt(-log(rd1) / ctr[idx].prim[end]) * sin(2.0 * π * rd2)
                ptc_temp[no_particle_temp].v[1] += ctr[idx].prim[2]
                ptc_temp[no_particle_temp].x = ctr[idx].x + (rd - 0.5) * ctr[idx].dx
                ptc_temp[no_particle_temp].idx = 0

                if ptc_temp[no_particle_temp].v[1] < 0
                    ptc_temp[no_particle_temp].tb =
                        (ctr[idx].x - 0.5 * ctr[idx].dx - ptc_temp[no_particle_temp].x) /
                        ptc_temp[no_particle_temp].v[1]
                elseif ptc_temp[no_particle_temp].v[1] > 0
                    ptc_temp[no_particle_temp].tb =
                        (ctr[idx].x + 0.5 * ctr[idx].dx - ptc_temp[no_particle_temp].x) /
                        ptc_temp[no_particle_temp].v[1]
                else
                    ptc_temp[no_particle_temp].tb = 1.0
                end
            end

        end

        for i = 1:ng
            idx = KS.pSpace.nx + i

            no_particle_cell_temp = round(ctr[idx].w[idx] * ctr[idx].dx / KS.gas.m) |> Int
            for j = 1:no_particle_cell_temp
                no_particle_temp = no_particle_temp + 1

                rd1 = rand(3)
                rd2 = rand(3)
                rd = rand()

                ptc_temp[no_particle_temp].mass = KS.gas.m
                @. ptc_temp[no_particle_temp].v =
                    sqrt(-log(rd1) / ctr[idx].prim[end]) * sin(2.0 * π * rd2)
                ptc_temp[no_particle_temp].v[1] += ctr[idx].prim[2]
                ptc_temp[no_particle_temp].x = ctr[idx].x + (rd - 0.5) * ctr[idx].dx
                ptc_temp[no_particle_temp].idx = 0

                if ptc_temp[no_particle_temp].v[1] < 0
                    ptc_temp[no_particle_temp].tb =
                        (ctr[idx].x - 0.5 * ctr[idx].dx - ptc_temp[no_particle_temp].x) /
                        ptc_temp[no_particle_temp].v[1]
                elseif ptc_temp[no_particle_temp].v[1] > 0
                    ptc_temp[no_particle_temp].tb =
                        (ctr[idx].x + 0.5 * ctr[idx].dx - ptc_temp[no_particle_temp].x) /
                        ptc_temp[no_particle_temp].v[1]
                else
                    ptc_temp[no_particle_temp].tb = 1.0
                end
            end

        end

    end

    return nothing

end
