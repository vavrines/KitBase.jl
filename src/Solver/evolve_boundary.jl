function evolve_boundary!(
    KS::SolverSet,
    ctr::T1,
    a1face::T2,
    a2face::T2,
    dt,
    mode,
    bc,
) where {T1<:AbstractArray{ControlVolume2D2F,2},T2<:AbstractArray{Interface2D2F,2}}

    nx, ny, dx, dy = begin
        if KS.ps isa CSpace2D
            KS.ps.nr, KS.ps.nθ, KS.ps.dr, KS.ps.darc
        else
            KS.ps.nx, KS.ps.ny, KS.ps.dx, KS.ps.dy
        end
    end

    bcs = ifelse(bc isa Symbol, [bc, bc, bc, bc], bc)

    if bcs[1] == :maxwell
        @inbounds Threads.@threads for j = 1:ny
            vn = KS.vSpace.u .* a1face[1, j].n[1] .+ KS.vSpace.v .* a1face[1, j].n[2]
            vt = KS.vSpace.v .* a1face[1, j].n[1] .- KS.vSpace.u .* a1face[1, j].n[2]
            xc = (KS.ps.vertices[1, j, 1, 1] + KS.ps.vertices[1, j, 4, 1]) / 2
            yc = (KS.ps.vertices[1, j, 1, 2] + KS.ps.vertices[1, j, 4, 2]) / 2
            bcL = local_frame(KS.ib.bc(xc, yc), a1face[1, j].n[1], a1face[1, j].n[2])
            flux_boundary_maxwell!(
                a1face[1, j].fw,
                a1face[1, j].fh,
                a1face[1, j].fb,
                bcL, # left
                ctr[1, j].h,
                ctr[1, j].b,
                vn,
                vt,
                KS.vSpace.weights,
                KS.gas.K,
                dt,
                dy[1, j],
                1,
            )
            a1face[1, j].fw .=
                global_frame(a1face[1, j].fw, a1face[1, j].n[1], a1face[1, j].n[2])
        end
    end

    if bcs[2] == :maxwell
        @inbounds Threads.@threads for j = 1:ny
            vn = KS.vSpace.u .* a1face[nx+1, j].n[1] .+ KS.vSpace.v .* a1face[nx+1, j].n[2]
            vt = KS.vSpace.v .* a1face[nx+1, j].n[1] .- KS.vSpace.u .* a1face[nx+1, j].n[2]
            xc = (KS.ps.vertices[nx, j, 2, 1] + KS.ps.vertices[nx, j, 3, 1]) / 2
            yc = (KS.ps.vertices[nx, j, 2, 2] + KS.ps.vertices[nx, j, 3, 2]) / 2
            bcR = local_frame(KS.ib.bc(xc, yc), a1face[nx+1, j].n[1], a1face[nx+1, j].n[2])
            flux_boundary_maxwell!(
                a1face[nx+1, j].fw,
                a1face[nx+1, j].fh,
                a1face[nx+1, j].fb,
                bcR, # right
                ctr[nx, j].h,
                ctr[nx, j].b,
                vn,
                vt,
                KS.vSpace.weights,
                KS.gas.K,
                dt,
                dy[nx, j],
                -1,
            )
            a1face[nx+1, j].fw .=
                global_frame(a1face[nx+1, j].fw, a1face[nx+1, j].n[1], a1face[nx+1, j].n[2])
        end
    end

    if bcs[3] == :maxwell
        @inbounds Threads.@threads for i = 1:nx
            vn = KS.vSpace.u .* a2face[i, 1].n[1] .+ KS.vSpace.v .* a2face[i, 1].n[2]
            vt = KS.vSpace.v .* a2face[i, 1].n[1] .- KS.vSpace.u .* a2face[i, 1].n[2]
            xc = (KS.ps.vertices[i, 1, 1, 1] + KS.ps.vertices[i, 1, 2, 1]) / 2
            yc = (KS.ps.vertices[i, 1, 1, 2] + KS.ps.vertices[i, 1, 2, 2]) / 2
            bcD = local_frame(KS.ib.bc(xc, yc), a2face[i, 1].n[1], a2face[i, 1].n[2])

            flux_boundary_maxwell!(
                a2face[i, 1].fw,
                a2face[i, 1].fh,
                a2face[i, 1].fb,
                bcD, # left
                ctr[i, 1].h,
                ctr[i, 1].b,
                vn,
                vt,
                KS.vSpace.weights,
                KS.gas.K,
                dt,
                dx[i, 1],
                1,
            )
            a2face[i, 1].fw .=
                global_frame(a2face[i, 1].fw, a2face[i, 1].n[1], a2face[i, 1].n[2])
        end
    end

    if bcs[4] == :maxwell
        @inbounds Threads.@threads for i = 1:nx
            vn = KS.vSpace.u .* a2face[i, ny+1].n[1] .+ KS.vSpace.v .* a2face[i, ny+1].n[2]
            vt = KS.vSpace.v .* a2face[i, ny+1].n[1] .- KS.vSpace.u .* a2face[i, ny+1].n[2]
            xc = (KS.ps.vertices[i, ny, 3, 1] + KS.ps.vertices[i, ny, 4, 1]) / 2
            yc = (KS.ps.vertices[i, ny, 3, 2] + KS.ps.vertices[i, ny, 4, 2]) / 2
            bcU = local_frame(KS.ib.bc(xc, yc), a2face[i, ny+1].n[1], a2face[i, ny+1].n[2])

            flux_boundary_maxwell!(
                a2face[i, ny+1].fw,
                a2face[i, ny+1].fh,
                a2face[i, ny+1].fb,
                bcU, # right
                ctr[i, ny].h,
                ctr[i, ny].b,
                vn,
                vt,
                KS.vSpace.weights,
                KS.gas.K,
                dt,
                dx[i, ny],
                -1,
            )
            a2face[i, ny+1].fw .=
                global_frame(a2face[i, ny+1].fw, a2face[i, ny+1].n[1], a2face[i, ny+1].n[2])
        end
    end

end
