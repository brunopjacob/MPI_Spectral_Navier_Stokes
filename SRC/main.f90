program main

use mpi
use declare_variables
use read_input
use vtk
use constants
use incompressible_solver

implicit none

type(VTK_file_handle)	:: fd !VTK file handle


CALL MPI_INIT( ierr )
CALL MPI_COMM_RANK( MPI_COMM_WORLD, id, ierr )
CALL MPI_COMM_SIZE( MPI_COMM_WORLD, npu, ierr )
CALL MPI_GET_PROCESSOR_NAME(procname, namelength,ierr)


if (id==0)then
!Program initialization
write(*,*)'#################################SPECTRAL 3D##################################'
write(*,*)'#                                                                            #'
write(*,*)'#                     |                                                      #'
write(*,*)'#    SSSSS            | SPECTRAL 3D:   Incompressible Three-Dimensional      #'
write(*,*)'#    S       SPECTRAL |                Navier-Stokes Solver Using Fourier    #'
write(*,*)'#    SSSSS       3    |                Pseudospectral method                 #'
write(*,*)'#        S       D    | Version: 1.0                                         #'
write(*,*)'#    SSSSS            |                                                      #'
write(*,*)'#                     | Developed by: Bruno J., UCSB                         #'
write(*,*)'#                                                                            #'
write(*,*)'##############################################################################'
endif


!==========================================================
!CPU Elapsed Time (start counter)                         =
!==========================================================
call cpu_time(t_elapsed_start)


!==========================================================
!Read input parameters                                    =
!==========================================================
call i_o_variables()

!==========================================================
!Domain decomposition                                     =
!==========================================================
call mpi_decomposition()


!==========================================================
!Allocate fields					  =
!==========================================================
call allocate_fields()


!==========================================================
!Read restart file (if any)	      		          =
!==========================================================
call restart_read()


!==========================================================
!Define constants					  =
!==========================================================
call define_constants()


!==========================================================
!Evaluate dimensions					  =
!==========================================================
call evaluate_dimensions()


!==========================================================
!Initialize variables					  =
!==========================================================
call initialize_variables()

 
!==========================================================
!create mesh                                              =
!==========================================================
call create_mesh()


!==========================================================
!evaluate wavenumbers                                     =
!==========================================================
call evaluate_wavenumbers()


!==========================================================
!evaluate filter                                          =
!==========================================================
call evaluate_filter()


!==========================================================
!Evaluate time step dt, based on CFL condition            =
!==========================================================
call evaluate_timestep()


!==========================================================
!Evaluate initial conditions                              =
!==========================================================
if (restart_flag .eqv. .FALSE.) then
	
	call evaluate_initial_conditions_incompressible(ic_set,domain_set,x,yp,zp,Lx,Ly,&
	Lz,nx,ny,nz,nyp,nzp,icommy,icommz,npuy,npuz,mu,amplitude_perturb,u,v,w,max_u0,max_v0,max_w0)
!endif
	
endif


call MPI_BARRIER (MPI_COMM_WORLD, ierr)

!==========================================================
!Time integration 					  =
!==========================================================

!********************************************************************************
!start time loop:		

!===========================================================
!Incompressible solver                                     =
!===========================================================
			
	do tc=tc_restart,nt
	
		t=(tc-1)*dt
	
		call temporal_integration_incompressible(nx,ny,nz,nyp,nzp,x,yp,zp,tc,kx,kyp,kzp,k_2p,icommy,icommz,npuy,npuz,imag,mu,dt,alphark,betark,crk,&
		temp_integration_set,u,v,w,pressure,vort_1,vort_2,vort_3,vort_mag,qcrit,ab_u,ab_v,ab_w,dx,dy,dz)

		unew      = u
		vnew      = v
		wnew      = w

		!===========================================================
		!Filter solution fields                                	   =
		!===========================================================
		call filter_field(nx,ny,nz,nyp,nzp,filter,unew)
		call filter_field(nx,ny,nz,nyp,nzp,filter,vnew)
		call filter_field(nx,ny,nz,nyp,nzp,filter,wnew)
		call filter_field(nx,ny,nz,nyp,nzp,filter,pressure)
		
		!===========================================================
		!Filter vort_1, vort_2, vort_3 and Q                       =
		!===========================================================
		aux = vort_1
		call PZFFT3DV (aux, vort_1, nx, ny, nz, icommy, icommz, npuy, npuz, -1)
		aux = vort_2
		call PZFFT3DV (aux, vort_2, nx, ny, nz, icommy, icommz, npuy, npuz, -1)
		aux = vort_3
		call PZFFT3DV (aux, vort_3, nx, ny, nz, icommy, icommz, npuy, npuz, -1)
		aux = qcrit
		call PZFFT3DV (aux, qcrit, nx, ny, nz, icommy, icommz, npuy, npuz, -1)
		
		call filter_field(nx,ny,nz,nyp,nzp,filter,vort_1)
		call filter_field(nx,ny,nz,nyp,nzp,filter,vort_2)
		call filter_field(nx,ny,nz,nyp,nzp,filter,vort_3)
		call filter_field(nx,ny,nz,nyp,nzp,filter,qcrit)
		
		aux = vort_1
		call PZFFT3DV (aux, vort_1, nx, ny, nz, icommy, icommz, npuy, npuz, +1)
		aux = vort_2
		call PZFFT3DV (aux, vort_2, nx, ny, nz, icommy, icommz, npuy, npuz, +1)
		aux = vort_3
		call PZFFT3DV (aux, vort_3, nx, ny, nz, icommy, icommz, npuy, npuz, +1)
		aux = qcrit
		call PZFFT3DV (aux, qcrit, nx, ny, nz, icommy, icommz, npuy, npuz, +1)


		!===========================================================
		!evaluate results in Euclidean space (R)  		   =
		!===========================================================
		aux = u
		call PZFFT3DV (aux, unew, nx, ny, nz, icommy, icommz, npuy, npuz, +1)
		aux = v
		call PZFFT3DV (aux, vnew, nx, ny, nz, icommy, icommz, npuy, npuz, +1)
		aux = w
		call PZFFT3DV (aux, wnew, nx, ny, nz, icommy, icommz, npuy, npuz, +1)
		aux = pressure
		call PZFFT3DV (aux, pressure, nx, ny, nz, icommy, icommz, npuy, npuz, +1)
		


		!===========================================================
		!syncronize MPI data                     		   =
		!===========================================================
		call MPI_BARRIER (MPI_COMM_WORLD, ierr)
		
		
		
		
		!===========================================================
		!print results						   =
		!===========================================================
		if (mod(tc,freqsave)==0)then
			!generate fields to print
			call grid_to_print(unew,usave,nx,nyp,nzp,node)
			call grid_to_print(vnew,vsave,nx,nyp,nzp,node)
			call grid_to_print(wnew,wsave,nx,nyp,nzp,node)
			call grid_to_print(pressure,psave,nx,nyp,nzp,node)
			call grid_to_print(vort_1,vort_1_save,nx,nyp,nzp,node)
			call grid_to_print(vort_2,vort_2_save,nx,nyp,nzp,node)
			call grid_to_print(vort_3,vort_3_save,nx,nyp,nzp,node)
			call grid_to_print(qcrit,qcrit_save,nx,nyp,nzp,node)
		
			call VTK_open_file(PREFIX=TRIM(str2), proc_rank=id, num_procs=npu, FD=fd, VTK_COUNTER=vtk_counter)

			!call VTK_open_file(PREFIX=TRIM(str2), FD=fd, VTK_COUNTER=vtk_counter)
			call VTK_write_mesh(FD=fd, X=xsave, Y=ypsave, Z=zpsave)
			call VTK_write_var(FD=fd, NAME="Velocity", VX=usave, VY=vsave, VZ=wsave)
			call VTK_write_var(FD=fd, NAME="Pressure", FIELD=psave)
			call VTK_write_var(FD=fd, NAME="Vorticity", VX=vort_1_save, VY=vort_2_save, VZ=vort_3_save)
			call VTK_write_var(FD=fd, NAME="Q", FIELD=qcrit_save)
			call VTK_close_file(FD=fd,VTK_COUNTER=vtk_counter)
			if (id == 0) then
			write(*,*)'#############################################'
			write(*,*)'DATA SAVED FOR t = ', t+dt
			write(*,*)'#############################################'
			endif
		endif
		!===========================================================
		!print restart info  					   =
		!===========================================================
		if (mod(tc,freq_restart)==0)then
			tc_restart = tc
			open(unit = npu+3, status='replace',file = TRIM(str1)//'restart_info.bin',form='unformatted',access='stream')  ! create a new file, or overwrite an existing one
			write(npu+3) t,tc_restart,vtk_counter ! write the data in array x to the file
			close(npu+3) ! close the file
			
		        write(rank,'(i8)') id
			str3=TRIM(str1)//'restart_'//trim(adjustl(rank))//".bin"
			open(unit = id+2, status='replace',file = TRIM(str3),form='unformatted',access='stream')  ! create a new file, or overwrite an existing one
			write(id+2) u,v,w  ! write the data in array x to the file
			close(id+2) ! close the file
		endif
		
                !===========================================================
		!print statistics   					   =
		!===========================================================
		if (mod(tc,freq_spectrum)==0)then
			call evaluate_spectrum(nx,nyp,nzp,n3,k0,dk,kmax,k_2p,freq_spectrum,tc,u,v,w)
		endif
		
		!===========================================================
		!print screen info  					   =
		!===========================================================
		if (id == 0) then
		if (mod(tc,freqsave_screen)==0)then
		
			max_u    = maxval(dabs(real(unew)))
			max_v    = maxval(dabs(real(vnew)))
			max_w    = maxval(dabs(real(wnew)))
			
			call screen_info(max_u,max_v,max_w,tc,t,dt,solver_set)
		
		endif
		endif
	enddo

!end time loop
!********************************************************************************
!call VTK_collect_file(FD=fd)

call MPI_COMM_FREE(icommy,ierr)
call MPI_COMM_FREE(icommz,ierr)
call MPI_FINALIZE(ierr) 


!==========================================================
!CPU Elapsed Time (pause counter)                         =
!==========================================================
call cpu_time(t_elapsed_end)

write(*,*)'*********************************************'
write(*,*)'SIMULATION COMPLETED'
write(*,*)'CPU ELAPSED TIME = ',t_elapsed_end-t_elapsed_start
write(*,*)'*********************************************'

end program
!===============================================================================================



