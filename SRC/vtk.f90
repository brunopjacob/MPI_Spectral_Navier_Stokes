!##########################################################
      MODULE vtk
!##########################################################
!	This module is responsibe to binary VTK results 
!
!##########################################################
	
	use declare_variables
	use mpi

	IMPLICIT NONE

	private

	character(len=1), parameter :: newline=achar(10)
	integer                     :: iproc=0, nb_procs=1
	type, public :: VTK_file_handle
	private
	character(len=80) :: prefix
	integer           :: unit
	integer           :: ni, nj, nk
	integer           :: counter=0
	integer           :: restart=0
	logical           :: first=.true.
	end type VTK_file_handle

	private :: handle_error, handle_warning, handle_info

	private :: VTK_write_mesh_2d,   VTK_write_mesh_3d,   &
		   VTK_write_vector_2d, VTK_write_vector_3d, &
		   VTK_write_scalar_2d, VTK_write_scalar_3d

	public :: VTK_open_file,  &
		  VTK_write_mesh, &
		  VTK_write_var,  &
		  VTK_close_file, &
		  VTK_collect_file

	interface VTK_write_mesh
	  module procedure VTK_write_mesh_2d, VTK_write_mesh_3d
	end interface

	interface VTK_write_var
	  module procedure VTK_write_scalar_2d, VTK_write_scalar_3d, &
			   VTK_write_vector_2d, VTK_write_vector_3d
	end interface
 
	contains

	subroutine handle_error(name, message)
    
	IMPLICIT NONE

	character(len=*), intent(in) :: name, message
	integer :: code, errcode=92

	print '(/,"   *** Error *** ", A,": ", A,/)',name, message
	stop
	end subroutine handle_error

	subroutine handle_warning(name, message)
    
	IMPLICIT NONE
    
	character(len=*), intent(in) :: name, message

	print '(/,"   *** Warning *** ",A,": ", A,/)',name, message
	end subroutine handle_warning

	subroutine handle_info(name, message)
	
	IMPLICIT NONE
	
	character(len=*), intent(in) :: name, message

	print '(/,"   *** Info *** ",A,": ", A,/)',name, message
	end subroutine handle_info

	subroutine VTK_open_file(prefix, proc_rank, num_procs, restart, fd, vtk_counter)
    
	IMPLICIT NONE
	
	character(len=*), intent(in)         :: prefix
	integer, optional, intent(in)        :: proc_rank, num_procs, restart
	type(VTK_file_handle), intent(inout) :: fd
	integer,intent(in):: vtk_counter  !!!
	character(len=10) :: rank, snapshot
	character(len=100) :: fluid_results_path='../RESULTS/FLUID_FIELD/'
	character(len=80) :: f
	character(len=256):: MAIN_header
	integer           :: err
	logical           :: file_opened

	
	fd%counter = vtk_counter
	
	!... Looking for a none connected logical file unit.
	fd%prefix=trim(prefix)
	fd%unit=99
	inquire(unit=fd%unit, opened=file_opened)
	do while (file_opened .and. fd%unit /= 0)

	  fd%unit = fd%unit - 1
	  inquire(unit=fd%unit, opened=file_opened)

	end do
    
	if (fd%unit == 0 .and. file_opened) then

	  call handle_error("VTK_open_file","All file units from 0 to 99 are already connected.")
	  stop

	else

	if ( present(proc_rank) .and. present(num_procs) ) then

	  iproc    = proc_rank
	  nb_procs = num_procs

	else if ( present(proc_rank) ) then

	  call handle_error("VTK_open_file","Both PROC_RANK and NUM_PROCS arguments must be present.")

	else if ( present(num_procs) ) then

	  call handle_error("VTK_open_file","Both PROC_RANK and NUM_PROCS arguments must be present.")

	end if
	
	!if (restart_flag == 1) then
	!	fd%counter=fd%counter
	!	!restart_flag == 0
	!endif
      
	if ((fd%first) .and. (present(restart))) then
	  fd%restart=restart
	  fd%counter=restart
	  fd%first=.false.
	end if
      
	fd%counter=fd%counter+1
	WRITE(snapshot,'(i8)') fd%counter
      
	if ( present(proc_rank) ) then
	  WRITE(rank,'(i8)') iproc
	  f=trim(fd%prefix)//"_"//trim(adjustl(rank))//"_"//trim(adjustl(snapshot))//".vtk"
	else
	  f=trim(fd%prefix)//"_"//trim(adjustl(snapshot))//".vtk"
	end if
    !  open(unit=fd%unit, file=trim(adjustl(f)), form="UNFORMATTED", access="STREAM", status="replace", &
    !       action="WRITE", iostat=err) 
	
    !open(unit=fd%unit,ACCESS='SEQUENTIAL',FORM='BINARY',file=trim(fluid_results_path)//trim(adjustl(f)), status="REPLACE",convert='big_endian')
	open(unit=fd%unit,ACCESS='STREAM',FORM='UNFORMATTED',file=trim(fluid_results_path)//trim(adjustl(f)), status="REPLACE",convert='big_endian')
	if(err /= 0) print '("Problem creating file ",a,".")', trim(f)
	end if
    
	MAIN_header="# vtk DataFile Version 3.0"//newline//"specNSF"//newline//"BINARY"//newline
	WRITE(unit=fd%unit) trim(MAIN_header)
  
	end subroutine VTK_open_file

	subroutine VTK_write_mesh_2d(fd, x, y)
    
	IMPLICIT NONE

	type(VTK_file_handle), intent(inout)   :: fd
	real(kind=rp), intent(in), dimension(:) :: x, y
	character(len=30)  :: buf1, buf2
	character(len=256) :: GRID_header
	integer, parameter :: nk=1

	fd%ni=size(x) ; fd%nj=size(y)
	WRITE(buf1,'(i8," ",i8," ",i8)') fd%ni,fd%nj,fd%nk
	GRID_header="DATASET RECTILINEAR_GRID"//newline//"DIMENSIONS "//trim(adjustl(buf1))//newline
	WRITE(unit=fd%unit) trim(GRID_header)
	WRITE(buf2,'(i8)') fd%ni
	GRID_header="X_COORDINATES "//trim(adjustl(buf2))//" float"//newline
	WRITE(unit=fd%unit) trim(GRID_header),real(x(1:fd%ni),kind=rp),newline
	WRITE(buf2,'(i8)') fd%nj
	GRID_header="Y_COORDINATES "//trim(adjustl(buf2))//" float"//newline
	WRITE(unit=fd%unit) trim(GRID_header),real(y(1:fd%nj),kind=rp),newline
	WRITE(buf2,'(i8)') fd%nk
	GRID_header="Z_COORDINATES "//trim(adjustl(buf2))//" float"//newline
	WRITE(unit=fd%unit) trim(GRID_header),real(0.0_rp,kind=rp),newline
	WRITE(buf2,'(i8)') fd%ni*fd%nj
	GRID_header="POINT_DATA "//trim(adjustl(buf2))//newline
	WRITE(unit=fd%unit) trim(GRID_header)

	end subroutine VTK_write_mesh_2d

	subroutine VTK_write_mesh_3d(fd, x, y, z)
    
	IMPLICIT NONE

	type(VTK_file_handle), intent(inout)   :: fd
	real(kind=rp), intent(in), dimension(:) :: x, y, z
	character(len=30)  :: buf1, buf2
	character(len=256) :: GRID_header
	

	fd%ni=size(x) ; fd%nj=size(y) ; fd%nk=size(z)
	WRITE(buf1,'(i8," ",i8," ",i8)') fd%ni,fd%nj,fd%nk
	GRID_header="DATASET RECTILINEAR_GRID"//newline//"DIMENSIONS "//trim(adjustl(buf1))//newline
	WRITE(unit=fd%unit) trim(GRID_header)
	WRITE(buf2,'(i8)') fd%ni
	GRID_header="X_COORDINATES "//trim(adjustl(buf2))//" float"//newline
	WRITE(unit=fd%unit) trim(GRID_header),real(x(1:fd%ni),kind=sp),newline	
	WRITE(buf2,'(i8)') fd%nj
	GRID_header="Y_COORDINATES "//trim(adjustl(buf2))//" float"//newline
	WRITE(unit=fd%unit) trim(GRID_header),real(y(1:fd%nj),kind=sp),newline
	WRITE(buf2,'(i8)') fd%nk
	GRID_header="Z_COORDINATES "//trim(adjustl(buf2))//" float"//newline
	WRITE(unit=fd%unit) trim(GRID_header),real(z(1:fd%nk),kind=sp),newline
	WRITE(buf2,'(i8)') fd%ni*fd%nj*fd%nk
	GRID_header="POINT_DATA "//trim(adjustl(buf2))//newline
	WRITE(unit=fd%unit) trim(GRID_header)
    
	end subroutine VTK_write_mesh_3d

	subroutine VTK_write_vector_2d(fd, name, vx, vy)
    
	IMPLICIT NONE

	type(VTK_file_handle), intent(in)        :: fd
	character(len=*), intent(in)             :: name
	real(kind=rp), intent(in), dimension(:,:) :: vx, vy
	real(kind=rp), allocatable, dimension(:,:,:) :: velocity
	integer                                     :: i, j, code=0
	character(len=256)                          :: uname, vname, VAR_header

	if ((size(vx,dim=1) /= fd%ni) .or. &
	    (size(vx,dim=2) /= fd%nj)) call handle_warning("VTK_WRITE_var","Incompatible X component and mesh sizes.")

	if ((size(vy,dim=1) /= fd%ni) .or. &
	    (size(vy,dim=2) /= fd%nj)) call handle_warning("VTK_WRITE_var","Incompatible Y component and mesh sizes.")

	if (.not.allocated(velocity)) then

	  allocate(velocity(3,fd%ni,fd%nj),STAT=code)
	  if ( code /= 0 ) &
	  call handle_error("VTK_WRITE_var","Not enough memory to allocate VELOCITY array")

	end if

	do j=1, fd%nj
	  do i=1, fd%ni

	      velocity(1,i,j) = vx(i,j)
	      velocity(2,i,j) = vy(i,j)
	      velocity(3,i,j) = 0.0_rp

	  end do
	end do

	VAR_header="VECTORS "//trim(adjustl(name))//" float "//newline
	WRITE(unit=fd%unit) trim(VAR_header),real(velocity(:,1:fd%ni,1:fd%nj),kind=sp),newline

	uname="X_"//name
	VAR_header="SCALARS "//trim(adjustl(uname))//" float "//newline//"LOOKUP_TABLE default"//newline
	  WRITE(unit=fd%unit) trim(VAR_header),real(velocity(1,1:fd%ni,1:fd%nj),kind=sp),newline

	vname="Y_"//name
	VAR_header="SCALARS "//trim(adjustl(vname))//" float "//newline//"LOOKUP_TABLE default"//newline
	  WRITE(unit=fd%unit) trim(VAR_header),real(velocity(2,1:fd%ni,1:fd%nj),kind=sp),newline

	if (allocated(velocity)) deallocate(velocity)

	end subroutine VTK_write_vector_2d

	subroutine VTK_write_vector_3d(fd, name, vx, vy, vz)
	IMPLICIT NONE

	type(VTK_file_handle), intent(in)          :: fd
	character(len=*), intent(in)               :: name
	real(kind=rp), intent(in), dimension(:,:,:) :: vx, vy, vz
	real(kind=rp), allocatable, dimension(:,:,:,:) :: velocity
	integer                                       :: i, j, k, code=0
	character(len=256)                            :: uname, vname, wname, VAR_header

	if ((size(vx,dim=1) /= fd%ni) .or. &
	    (size(vx,dim=2) /= fd%nj) .or. &
	    (size(vx,dim=3) /= fd%nk)) call handle_warning("VTK_WRITE_var","Incompatible X component and mesh sizes.")

	if ((size(vy,dim=1) /= fd%ni) .or. &
	    (size(vy,dim=2) /= fd%nj) .or. &
	    (size(vy,dim=3) /= fd%nk)) call handle_warning("VTK_WRITE_var","Incompatible Y component and mesh sizes.")

	if ((size(vz,dim=1) /= fd%ni) .or. &
	    (size(vz,dim=2) /= fd%nj) .or. &
	    (size(vz,dim=3) /= fd%nk)) call handle_warning("VTK_WRITE_var","Incompatible Z component and mesh sizes.")

	if (.not.allocated(velocity)) then
	  allocate(velocity(3,fd%ni,fd%nj,fd%nk),STAT=code)
	  if ( code /= 0 ) &
	  call handle_error("VTK_WRITE_var","Not enough memory to allocate VELOCITY array")
	end if

	do k=1, fd%nk
	  do j=1, fd%nj
	    do i=1, fd%ni

	      velocity(1,i,j,k) = vx(i,j,k)
	      velocity(2,i,j,k) = vy(i,j,k)
	      velocity(3,i,j,k) = vz(i,j,k)

	    end do
	  end do
	end do
      
	VAR_header="VECTORS "//trim(adjustl(name))//" float "//newline
	WRITE(unit=fd%unit) trim(VAR_header),real(velocity(:,1:fd%ni,1:fd%nj,1:fd%nk),kind=sp),newline

	uname="X_"//name
	VAR_header="SCALARS "//trim(adjustl(uname))//" float "//newline//"LOOKUP_TABLE default"//newline
	  WRITE(unit=fd%unit) trim(VAR_header),real(velocity(1,1:fd%ni,1:fd%nj,1:fd%nk),kind=sp),newline

	vname="Y_"//name
	VAR_header="SCALARS "//trim(adjustl(vname))//" float "//newline//"LOOKUP_TABLE default"//newline
	  WRITE(unit=fd%unit) trim(VAR_header),real(velocity(2,1:fd%ni,1:fd%nj,1:fd%nk),kind=sp),newline

	wname="Z_"//name
	VAR_header="SCALARS "//trim(adjustl(wname))//" float "//newline//"LOOKUP_TABLE default"//newline
	  WRITE(unit=fd%unit) trim(VAR_header),real(velocity(3,1:fd%ni,1:fd%nj,1:fd%nk),kind=sp),newline

	if (allocated(velocity)) deallocate(velocity)

	end subroutine VTK_write_vector_3d
	

	subroutine VTK_write_scalar_2d(fd, name, field)
    
	IMPLICIT NONE

	type(VTK_file_handle), intent(in)        :: fd
	character(len=*), intent(in)             :: name
	real(kind=rp), intent(in), dimension(:,:) :: field
	character(len=256) :: VAR_header

	if ((size(field,dim=1) /= fd%ni) .or. (size(field,dim=2) /= fd%nj)) &
	  call handle_warning("VTK_WRITE_var","Incompatible FIELD and MESH sizes.")

	VAR_header="SCALARS "//trim(adjustl(name))//" float "//newline//"LOOKUP_TABLE default"//newline
	WRITE(unit=fd%unit) trim(VAR_header),real(field(1:fd%ni,1:fd%nj),kind=sp),newline

	end subroutine VTK_write_scalar_2d

	subroutine VTK_write_scalar_3d(fd, name, field)
    
	IMPLICIT NONE

	type(VTK_file_handle), intent(in)          :: fd
	character(len=*), intent(in)               :: name
	real(kind=rp), intent(in), dimension(:,:,:) :: field
	character(len=256) :: VAR_header

	if ((size(field,dim=1) /= fd%ni) .or. &
	    (size(field,dim=2) /= fd%nj) .or. &
	    (size(field,dim=3) /= fd%nk)) call handle_warning("VTK_WRITE_var","Incompatible FIELD and MESH sizes.")

	VAR_header="SCALARS "//trim(adjustl(name))//" float "//newline//"LOOKUP_TABLE default"//newline
	WRITE(unit=fd%unit) trim(VAR_header),real(field(1:fd%ni,1:fd%nj,1:fd%nk),kind=sp),newline

	end subroutine VTK_write_scalar_3d

	subroutine VTK_close_file(fd,vtk_counter)
   
	IMPLICIT NONE

	type(VTK_file_handle), intent(in) :: fd
	integer :: vtk_counter
	logical :: file_opened

	inquire(unit=fd%unit, opened=file_opened)
	if (file_opened) then
	  close(unit=fd%unit)
	else
	  call handle_warning("VTK_close_file","No such file to close. Please, check file descriptor.")
	end if

	vtk_counter = fd%counter
	end subroutine VTK_close_file

	subroutine VTK_collect_file(fd)

	IMPLICIT NONE
	
	type(VTK_file_handle), intent(inout) :: fd
	character(len=10) :: rank, snapshot
	character(len=80) :: f, vtrfile
	character(len=100) :: fluid_results_path='../RESULTS/FLUID_FIELD/'
	integer           :: shot, code, err, nt, np, k
	logical           :: file_opened

	!... Looking for a none connected logical file unit.
	if (iproc == 0) then
	  fd%unit=99
	  inquire(unit=fd%unit, opened=file_opened)
	  do while (file_opened .and. fd%unit /= 0)
	    fd%unit = fd%unit - 1
	    inquire(unit=fd%unit, opened=file_opened)
	  end do
	  if (fd%unit == 0 .and. file_opened) then
	    call handle_warning("VTK_open_file","warning, all file units from 0 to 99 are already connected.")
	  else
	    f=trim(adjustl(fd%prefix))//".pvd"
	    open(unit=fd%unit, file=trim(fluid_results_path)//trim(adjustl(f)), form="FORMATTED", status="replace", &
	      action="WRITE", iostat=err)
	    if(err /= 0) print '("VTK_collect_file: Error, problem creating file ",a,".")', trim(f)
	    WRITE(unit=fd%unit,fmt='(100A)')   '<?xml version="1.0"?>'
	    WRITE(unit=fd%unit,fmt='(100A)')   '<VTKFile type="Collection" version="0.1" format="ascii">'
	    WRITE(unit=fd%unit,fmt='(100A)')   '  <Collection>'
	    nt=len_trim(fd%prefix)
	    np=scan(STRING=fd%prefix, SET="/", BACK=.true.)
	    vtrfile=fd%prefix(np+1:nt)
	    if ( nb_procs == 1 ) then

	      do shot = 1, fd%counter
		WRITE(snapshot, '(i6)') shot
		WRITE(unit=fd%unit,fmt='(100A)') '    <DataSet timestep="'//trim(adjustl(snapshot))//&
					  &'" part="0'//'" file="'//trim(adjustl(vtrfile))//&
					  &"_"//trim(adjustl(snapshot))//'.vtr"/>'
	      end do

	    else

	      do k = 0, nb_procs-1
		WRITE(rank, '(i6)') k
		do shot = 1, fd%counter
		  WRITE(snapshot, '(i6)') shot
		  WRITE(unit=fd%unit,fmt='(100A)') '    <DataSet timestep="'//trim(adjustl(snapshot))//&
					    &'" part="'//trim(adjustl(rank))//'" file="'//&
					    &trim(adjustl(vtrfile))//"_"//trim(adjustl(rank))//&
					    &"_"//trim(adjustl(snapshot))//'.vtr"/>'
		end do
	      end do

	    end if
	    WRITE(unit=fd%unit,fmt='(100A)')    '  </Collection>'
	    WRITE(unit=fd%unit,fmt='(100A)')    '</VTKFile>'
	    close(unit=fd%unit)
	  end if
	end if
	fd%counter=0 ; fd%restart=0 ; fd%first=.true. ; iproc=0 ; nb_procs=1
	end subroutine VTK_collect_file

      END MODULE vtk
