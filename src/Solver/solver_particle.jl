# ============================================================
# Update Algorithm for Particle Simulations
# ============================================================

function update!(
    KS::T, 
    ctr::ControlVolumeParticle1D, 
    ptc::Particle1D, 
    face::Interface1D, 
    dt, 
    residual; 
    coll = :bgk::Symbol,
    bc = :fix::Symbol,
) where {T<:AbstractSolverSet}

    np = KS.gas.np

    KS.pSpace.nx

#=


    allocate(rand_time(no_particle))

    call RANDOM_NUMBER(rand_time)

    forall(i=ixmin:ixmax)
        ctr(i)%wg=ctr(i)%w
    endforall

    no_particle_temp=0
 
    for i=1:np
        
    
    
    
    
    
    
    
    
    
    
        cell_no=particle(i)%cell_no
        time_c=-ctr(cell_no)%tau*log(rand_time(i))
        
        
        
        
        
        
        
        
        
        
        
        if(time_c>dt)then
            particle(i)%x=particle(i)%x+dt*particle(i)%v(1)
            cell_no_temp=particle(i)%cell_no
            ctr(cell_no)%wg(1)=ctr(cell_no)%wg(1)-particle(i)%mass/dx
            ctr(cell_no)%wg(2)=ctr(cell_no)%wg(2)-particle(i)%mass*particle(i)%v(1)/dx
            ctr(cell_no)%wg(3)=ctr(cell_no)%wg(3)-0.5*particle(i)%mass*sum(particle(i)%v(:)**2)/dx
            if(particle(i)%x>0 .and. particle(i)%x<50)then
                no_particle_temp=no_particle_temp+1
                particle_temp(no_particle_temp)%mass=particle(i)%mass
                particle_temp(no_particle_temp)%x=particle(i)%x
                particle_temp(no_particle_temp)%v=particle(i)%v
                cell_no_temp=ceiling(particle(i)%x/dx)
                particle_temp(no_particle_temp)%cell_no=cell_no_temp
                if(particle_temp(no_particle_temp)%v(1)<0)then
                    particle_temp(no_particle_temp)%time_b=(ctr(cell_no_temp)%x-&
                        0.5*dx-particle_temp(no_particle_temp)%x)/particle_temp(no_particle_temp)%v(1)
                elseif(particle_temp(no_particle_temp)%v(1)>0)then
                    particle_temp(no_particle_temp)%time_b=(ctr(cell_no_temp)%x+&
                        0.5*dx-particle_temp(no_particle_temp)%x)/particle_temp(no_particle_temp)%v(1)
                else
                    particle_temp(no_particle_temp)%time_b=1
                endif
            endif
        elseif(time_c>particle(i)%time_b)then
            particle(i)%x=particle(i)%x+dt*particle(i)%v(1)
            cell_no_temp=particle(i)%cell_no
            ctr(cell_no)%wg(1)=ctr(cell_no)%wg(1)-particle(i)%mass/dx
            ctr(cell_no)%wg(2)=ctr(cell_no)%wg(2)-particle(i)%mass*particle(i)%v(1)/dx
            ctr(cell_no)%wg(3)=ctr(cell_no)%wg(3)-0.5*particle(i)%mass*sum(particle(i)%v(:)**2)/dx
            if(particle(i)%x>0 .and. particle(i)%x<50)then
                cell_no_temp=ceiling(particle(i)%x/dx)
                ctr(cell_no_temp)%wg(1)=ctr(cell_no_temp)%wg(1)+particle(i)%mass/dx
                ctr(cell_no_temp)%wg(2)=ctr(cell_no_temp)%wg(2)+particle(i)%mass*particle(i)%v(1)/dx
                ctr(cell_no_temp)%wg(3)=ctr(cell_no_temp)%wg(3)+0.5*particle(i)%mass*sum(particle(i)%v(:)**2)/dx
            endif
        endif
    enddo

    forall(i=ixmin:ixmax)
        ctr(i)%w=0
    endforall
    do i=1,no_particle_temp
        particle_cell_temp=particle_temp(i)%cell_no
        ctr(particle_cell_temp)%w(1)=ctr(particle_cell_temp)%w(1)+particle_temp(i)%mass/dx
        ctr(particle_cell_temp)%w(2)=ctr(particle_cell_temp)%w(2)+particle_temp(i)%mass*particle_temp(i)%v(1)/dx
        ctr(particle_cell_temp)%w(3)=ctr(particle_cell_temp)%w(3)+0.5*particle_temp(i)%mass*sum(particle_temp(i)%v(:)**2)/dx
    enddo
=#

end