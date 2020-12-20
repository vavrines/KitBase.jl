# ============================================================
# Update Algorithm for Particle Simulations
# ============================================================

"""
    sample_particle!(ptc::Particle1D, m, x, v, idx, tb)
    sample_particle!(ptc::Particle1D, KS, ctr, idx, t0 = 1.0)

Sample particles from local flow conditions

"""
function sample_particle!(ptc::Particle1D, m, x, v, e, idx, flag, tc)
    ptc.m = m
    ptc.x = x
    ptc.v = v
    ptc.e = e
    ptc.idx = idx
    ptc.flag = flag
    ptc.tc = tc

    return nothing
end

function sample_particle!(ptc::Particle1D, m, prim::T, x, dx, idx, μᵣ, ω, flag = 0) where {T<:AbstractArray{<:Real,1}}
    ptc.m = m
    ptc.x = x + (rand() - 0.5) * dx
    ptc.v = sample_maxwell(prim)
    ptc.e = 0.5 / prim[end]
    ptc.idx = idx
    ptc.flag = flag
    τ = vhs_collision_time(prim, μᵣ, ω)
    ptc.tc = next_collision_time(τ)

    return nothing
end

function sample_particle!(ptc::Particle1D, m, prim::T, umin, umax, x, dx, idx, μᵣ, ω, flag = 0) where {T<:AbstractArray{<:Real,1}}
    ptc.m = m
    ptc.x = x + (rand() - 0.5) * dx
    ptc.v = sample_maxwell(prim, umin, umax)
    ptc.e = 0.5 / prim[end]
    ptc.idx = idx
    ptc.flag = flag
    τ = vhs_collision_time(prim, μᵣ, ω)
    ptc.tc = next_collision_time(τ)

    return nothing
end

function sample_particle!(ptc::Particle1D, KS, ctr, idx)
    ptc.m = KS.gas.m
    ptc.x = ctr.x + (rand() - 0.5) * ctr.dx
    ptc.v = sample_maxwell(ctr.prim)
    ptc.e = 0.5 / ctr.prim[end]
    ptc.idx = idx
    ptc.flag = 0
    τ = vhs_collision_time(ctr.prim, KS.gas.μᵣ, KS.gas.ω)
    ptc.tc = next_collision_time(τ)

    return nothing
end


"""
    free_transport!(KS, x, v, flag, dt, np=length(x))

Flight particles along trajectories

"""
function transport!(KS, x, v, flag, dt, np = length(x))

    @inbounds for i in 1:np
        x[i] += v[i,1]*dt

        flag[i] = 0
        if x[i] < KS.pSpace.x0
            flag[i] = 1
        elseif x[i] > KS.pSpace.x1
            flag[i] = 2
        end
        
        if flag[i] != 0
            x0 = x[i] - v[i, 1]*dt
            vi = @view v[i, :]
            x[i] = boundary!(x[i], vi, (KS.ib, [KS.pSpace.x0, KS.pSpace.x1], x0, flag[i], dt))
        end
    end

    return nothing

end


function boundary!(x, v, p; bc = :maxwell::Symbol)

    ib, xw, x0, flag, dt = p

    primB = [ib.primL, ib.primR]
    bound = [[0., Inf], [-Inf, 0.]]
    vw = [ib.vL, ib.vR]

    if bc == :maxwell
        x = maxwell_reflection!(x, v, (primB[flag], bound[flag], xw[flag], vw[flag], x0, dt))
    else
        throw("unavailable boundary condition")
    end

    return x

end


"""
    particle_boundary_maxwell!(x, v, p)

Maxwell diffusive boundary for particle reflection

"""
function maxwell_reflection!(x, v, p)

    primw, bound, xw, vw, x0, dt = p

    #v[1] = sqrt(-log(1.0-rand()))
    v[1] = sample_maxwell(primw[end], bound[1], bound[2])
    v[2] = sample_maxwell(primw[end], vw[2])
    v[3] = sample_maxwell(primw[end], vw[3])
    
    dtr = dt * (x - xw) / (x - x0)
    x = xw + v[1] * dtr

    return x

end


function sort!(KS, ctr::T, x, idx, ref, np=length(idx); mode=:uniform) where {T<:AbstractArray{<:AbstractControlVolume1D,1}}

    # calculate cell indices of particles
    @inbounds for i in 1:np
        idx[i] = find_idx(KS.pSpace.x[1:end], x[i], mode=mode)
    end

    # count the number of particles in each cell
    @inbounds for i in eachindex(ctr)
        ctr[i].np = 0
    end
    for i in 1:np
        ctr[idx[i]].np += 1
    end

    # build index list as cumulative sum of the number of particles in each cell
    m = 1
    @inbounds for jcell in eachindex(ctr)
        ctr[jcell].ip = m
        m += ctr[jcell].np
    end

    # build cross-reference list
    temp = zeros(Int, axes(ctr))
    @inbounds for i=1:np
        jcell = idx[i]
        k = ctr[jcell].ip + temp[jcell]
        ref[k] = i
        temp[jcell] += 1
    end

    return nothing

end


function dsmc!(KS, ctr, ref, v, dt, ne=1)

    vcm = zeros(3)
    vrel = zeros(3)
    col = 0

    for jcell in eachindex(ctr)

        number = ctr[jcell].np
        if number > 1

            #select = coeff*number**2*vrmax(jcell) + remainder(jcell)
            select = ne * 0.5 / √2 * dt / ctr[jcell].dx / KS.gas.Kn * number^2*ctr[jcell].vrmax + ctr[jcell].remainder

            nsel = Int(floor(select))
            ctr[jcell].remainder = select-nsel 
            crm = ctr[jcell].vrmax

            for isel=1:nsel

                # pick two particles at random out of this cell
                k = floor( rand()*number ) |> Int
                kk = mod( Int(floor(k+rand()*(number-1))+1), number )
                ip1 = ref[ k+ctr[jcell].ip ]  # first particle
                ip2 = ref[ kk+ctr[jcell].ip ] # second particle

                # calculate pair's relative speed
                cr = sqrt( (v[ip1,1]-v[ip2,1])^2 + (v[ip1,2]-v[ip2,2])^2 + (v[ip1,3]-v[ip2,3])^2 )
                if cr > crm    # If relative speed larger than crm,
                    crm = cr                # then reset crm to larger value
                end

                # accept or reject candidate pair according to relative speed
                if cr/ctr[jcell].vrmax > rand()
                    # If pair accepted, select post-collision velocities
                    col += 1                     # Collision counter
                    for k=1:3
                        vcm[k] = 0.5*(v[ip1,k] + v[ip2,k])       # Center of mass velocity
                    end
                    cos_th = 1.0 - 2.0*rand()       # Cosine and sine of
                    sin_th = sqrt(1.0 - cos_th^2)      # collision angle theta
                    phi = 2.0*pi*rand()             # Collision angle phi
                    vrel[1] = cr*cos_th                 # Compute post-collision
                    vrel[2] = cr*sin_th*cos(phi)        # relative velocity
                    vrel[3] = cr*sin_th*sin(phi)
                    for  k=1:3
                        v[ip1,k] = vcm[k] + 0.5*vrel[k]   # Update post-collision
                        v[ip2,k] = vcm[k] - 0.5*vrel[k]   # velocities
                    end

                end
                ctr[jcell].vrmax = crm     # Update max relative speed

            end
        end
    end

    return col

end


function stat!(KS, ctr, ptc)
    @inbounds Threads.@threads for i in 1:KS.pSpace.nx
        ctr[i].w .= 0.
    end

    @inbounds Threads.@threads for i in eachindex(ptc.x)
        ctr[ptc.idx[i]].w[1] += ptc.m[i] / ctr[ptc.idx[i]].dx
        ctr[ptc.idx[i]].w[2] += ptc.m[i] * ptc.v[i, 1] / ctr[ptc.idx[i]].dx
        ctr[ptc.idx[i]].w[3] += ptc.m[i] * ptc.v[i, 2] / ctr[ptc.idx[i]].dx
        ctr[ptc.idx[i]].w[4] += 0.5 * ptc.m[i] * sum(ptc.v[i, :] .^ 2) / ctr[ptc.idx[i]].dx
    end

    @inbounds Threads.@threads for i in 1:KS.pSpace.nx
        ctr[i].prim .= KitBase.conserve_prim(ctr[i].w, KS.gas.γ)
    end
end


function update!(
    KS::AbstractSolverSet,
    ctr::AbstractArray{ControlVolumeParticle1D,1},
    ptc::AbstractArray{Particle1D,1},
    ptc_temp::AbstractArray{Particle1D,1},
    face::AbstractArray{Interface1D,1},
    dt,
    res;
    coll = :bgk::Symbol,
    bc = :fix::Symbol,
)

    res .= 0.

    particle_transport!(KS, ctr, ptc, ptc_temp, dt)
    particle_collision!(KS, ctr, ptc_temp, face, res, coll)
    particle_boundary!(KS, ctr, ptc, ptc_temp, face, dt, coll, bc)
    duplicate!(ptc, ptc_temp, KS.gas.np)

    return nothing

end


"""
Update algorithm for particle transports

    update_transport!(KS::AbstractSolverSet, ctr::ControlVolumeParticle1D, ptc::Particle1D, ptc_temp::Particle1D, dt)

"""
function particle_transport!(
    KS::T1,
    ctr::T2,
    ptc::T3,
    ptc_tmp::T3,
    dt,
) where {
    T1<:AbstractSolverSet,
    T2<:AbstractArray{ControlVolumeParticle1D,1},
    T3<:AbstractArray{Particle1D,1},
}

    @inbounds Threads.@threads for i = 1:KS.pSpace.nx
        ctr[i].wf .= ctr[i].w
    end

    npt = 0
    @inbounds for i in 1:KS.gas.np
        cellid = ptc[i].idx
        ptc[i].tc = -vhs_collision_time(ctr[cellid].prim, KS.gas.μᵣ, KS.gas.ω) * log(rand())

        if ptc[i].tc > dt # class 1: free transport particles
            ptc[i].x += dt * ptc[i].v[1]

            ctr[cellid].wf[1] -= ptc[i].m / ctr[cellid].dx
            ctr[cellid].wf[2] -= ptc[i].m * ptc[i].v[1] / ctr[cellid].dx
            ctr[cellid].wf[3] -= 0.5 * ptc[i].m * sum(ptc[i].v .^ 2) / ctr[cellid].dx

            if ptc[i].x >= KS.pSpace.x0 && ptc[i].x <= KS.pSpace.x1
                npt += 1

                cellid_tmp = Int(ceil((ptc[i].x - KS.pSpace.x0) / ctr[1].dx))
                #cellid_tmp = find_idx(KS.pSpace.x[1:KS.pSpace.nx], ptc[i].x, mode=:uniform)

                sample_particle!(ptc_tmp[npt], ptc[i].m, ptc[i].x, ptc[i].v, 0.5 / ctr[cellid_tmp].prim[end], cellid_tmp, ptc[i].tc)
            end
        elseif ptc[i].tc > boundary_time(ptc[i].x, ptc[i].v, ctr[cellid].x, ctr[cellid].dx)   # class 2: jumping collision particles
            #ptc[i].x += ptc[i].tc * ptc[i].v[1]
            ptc[i].x += dt * ptc[i].v[1]

            if ptc[i].x > ctr[cellid].x + 0.5 * ctr[cellid].dx && ptc[i].x < ctr[cellid].x - 0.5 * ctr[cellid].dx

                ctr[cellid].wf[1] -= ptc[i].m / ctr[cellid].dx
                ctr[cellid].wf[2] -= ptc[i].m * ptc[i].v[1] / ctr[cellid].dx
                ctr[cellid].wf[3] -= 0.5 * ptc[i].m * sum(ptc[i].v .^ 2) / ctr[cellid].dx

                if ptc[i].x >= KS.pSpace.x0 && ptc[i].x <= KS.pSpace.x1
                    cellid_tmp = Int(ceil((ptc[i].x - KS.pSpace.x0) / ctr[1].dx)) 
                    #cellid_tmp = find_idx(KS.pSpace.x[1:KS.pSpace.nx], ptc[i].x, mode=:uniform)

                    ctr[cellid_tmp].wf[1] += ptc[i].m / ctr[cellid_tmp].dx
                    ctr[cellid_tmp].wf[2] += ptc[i].m * ptc[i].v[1] / ctr[cellid_tmp].dx
                    ctr[cellid_tmp].wf[3] +=
                        0.5 * ptc[i].m * sum(ptc[i].v .^ 2) / ctr[cellid_tmp].dx
                end
            end
        end
    end

    @inbounds Threads.@threads for i = 1:KS.pSpace.nx
        ctr[i].wp .= 0.0
    end
    @inbounds for i = 1:npt
        cellid = ptc_tmp[i].idx

        ctr[cellid].wp[1] += ptc_tmp[i].m / ctr[cellid].dx
        ctr[cellid].wp[2] += ptc_tmp[i].m * ptc_tmp[i].v[1] / ctr[cellid].dx
        ctr[cellid].wp[3] += 0.5 * ptc_tmp[i].m * sum(ptc_tmp[i].v .^ 2) / ctr[cellid].dx
    end

    KS.gas.np = npt
    return nothing

end


"""
Update algorithm for particle collisions

    update_collision!(KS::AbstractSolverSet, ctr::ControlVolumeParticle1D, ptc::Particle1D, ptc_temp::Particle1D, dt)

"""
function particle_collision!(
    KS::T,
    ctr::AbstractArray{ControlVolumeParticle1D,1},
    ptc_temp::AbstractArray{Particle1D,1},
    face::AbstractArray{Interface1D,1},
    res::AbstractArray{<:AbstractFloat,1},
    coll = :bgk::Symbol,
) where {T<:AbstractSolverSet}

    sum_res = zeros(3)
    sum_avg = zeros(3)

    npt = KS.gas.np
    @inbounds for i = 1:KS.pSpace.nx
        # fluid variables
        wold = deepcopy(ctr[i].w)
        @. ctr[i].wf += (face[i].fw - face[i+1].fw) / ctr[i].dx

        if ctr[i].wf[1] < 0.0
            ctr[i].wf .= zero(ctr[i].wf)
        #end
        else
            primG = conserve_prim(ctr[i].wf, KS.gas.γ)
            if primG[end] < 0.0
                primG[end] = 1e8
                ctr[i].wf .= prim_conserve(primG, KS.gas.γ)
            end
        end

        @. ctr[i].w = ctr[i].wf + ctr[i].wp
        if ctr[i].w[1] < 0.0
            ctr[i].w .= zero(ctr[i].w) .+ 1e-8
        else
            ctr[i].prim = conserve_prim(ctr[i].w, KS.gas.γ)
            if ctr[i].prim[end] < 0.0
                ctr[i].prim[end] = eltype(ctr[i].prim)(1e8)
                ctr[i].w .= prim_conserve(ctr[i].prim, KS.gas.γ)
            end
        end
        @. sum_res += (ctr[i].w - wold)^2
        @. sum_avg += abs.(ctr[i].w)

        ctr[i].τ = vhs_collision_time(ctr[i].prim, KS.gas.μᵣ, KS.gas.ω)
        
        # particles
        no_cellid = Int(round(ctr[i].wf[1] * ctr[i].dx / KS.gas.m))
        for j in 1:no_cellid
            npt += 1
            sample_particle!(ptc_temp[npt], KS.gas.m, primG, ctr[i].x, ctr[i].dx, i, KS.gas.μᵣ, KS.gas.ω)
            #sample_particle!(ptc_temp[npt], KS.gas.m, primG, KS.vSpace.u0, KS.vSpace.u1, ctr[i].x, ctr[i].dx, i, KS.gas.μᵣ, KS.gas.ω)
        end
    end

    @. res = sqrt(KS.pSpace.nx*sum_res)/(sum_avg+1e-8)
    KS.gas.np = npt
    return nothing


end


function particle_boundary!(
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
        # ghost cell
        ng = 1 - first(eachindex(KS.pSpace.x))

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


function duplicate!(ptc, ptc_new, n=length(ptc))
    @inbounds Threads.@threads for i in 1:n
        ptc[i].m = ptc_new[i].m
        ptc[i].x = ptc_new[i].x
        ptc[i].v .= ptc_new[i].v
        ptc[i].e = ptc_new[i].e
        ptc[i].idx = ptc_new[i].idx
        ptc[i].tc = ptc_new[i].tc
    end

    return nothing
end
