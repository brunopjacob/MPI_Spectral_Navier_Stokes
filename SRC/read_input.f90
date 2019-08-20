module read_input

use mpi
use declare_variables

public:: i_o_variables
private:: read_file
      

contains


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
subroutine i_o_variables()
   
   
	   open (unit=1, file='input_parameters.dat', status='old', action='read' )
		
  	   call read_file (line)
	   read (line, *) turbulence_modeling_set

	   call read_file (line)
	   read (line, *) ic_set
	   
!	   call read_file (line)
!	   read (line, *) domain_set
	   
	   call read_file (line)
	   read (line, *) nx
!	   if (domain_set==2) nx = 1
	
	   call read_file (line)
	   read (line, *) ny
	   
	   call read_file (line)
	   read (line, *) nz
	   
	   call read_file (line)
	   read (line, *) Lx
	
	   call read_file (line)
	   read (line, *) Ly
	   
	   call read_file (line)
	   read (line, *) Lz
		
	   call read_file (line)
	   read (line, *) temp_integration_set
	
	   call read_file (line)
	   read (line, *) temp_input_set
	   
	   call read_file (line)
	   read (line, *) filter_set
	      
	   call read_file (line)
	   read (line, *) mu
	   	   
	   call read_file (line)
	   read (line, *) amplitude_perturb
	   	   	   
	   call read_file (line)
	   read (line, *) dt
	   
	   call read_file (line)
	   read (line, *) CFL
	   
	   call read_file (line)
	   read (line, *) nt
	   
	   call read_file (line)
	   read (line, *) freqsave
	   
	   call read_file (line)
	   read (line, *) freqsave_screen
	   
	   call read_file (line)
	   read (line, *) freq_restart
	   
	   call read_file (line)
	   read (line, *) freq_spectrum
	     
	   call read_file (line)
	   read (line, *) str2
	   
	   close (1)
	   
end subroutine i_o_variables
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
subroutine read_file (line)


  implicit none
  character*80:: line
  integer::var_read
  
10 continue
   read (1, *) line
   if (line(1:1) .eq. '#') goto 10
   return
   
end subroutine read_file
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
subroutine restart_read()

	!Input file (or restart of the simulation)
	str1 = '../RESULTS/FLUID_FIELD/'
	inquire (file = TRIM(str1)//'restart_info.bin', exist = restart_flag)
	if (restart_flag .eqv. .TRUE.) then
		if (id == 0) then
			write(*,*)''
			write(*,*)'RESTART DATA FOUND ... LOADING DATA ...'
		endif
		open(unit = npu+3, status='old',file = TRIM(str1)//'restart_info.bin',form='unformatted',access='stream')  ! open an existing file
		read(npu+3) t,tc_restart,vtk_counter ! read the data into array x, of the appropriate data type
		close(npu+3) ! close the file
		tc_restart=tc_restart+1 !in order to do not evaluate twice the same time step
		
		
	        write(rank,'(i8)') id
	        str3=TRIM(str1)//'restart_'//trim(adjustl(rank))//".bin"
		open(unit = id+2, status='old',file = TRIM(str3),form='unformatted',access='stream')  ! open an existing file
		read(id+2) u,v,w ! read the data into array x, of the appropriate data type
		close(id+2) ! close the file
		
		
		 if (id == 0) then
		 write(*,*)'LOAD COMPLETED'
		 write(*,*)''
		 write(*,*)'RESTARTING FROM TIME STEP = ',tc_restart
		 write(*,*)'#############################################'
		 write(*,*)''
		 endif
	endif

end subroutine restart_read
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


end module read_input
