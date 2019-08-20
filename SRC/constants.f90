module constants

use mpi
use declare_variables


public::mpi_decomposition,allocate_fields,define_constants,evaluate_dimensions,&
initialize_variables,create_mesh,evaluate_wavenumbers,evaluate_timestep


contains


	!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	subroutine mpi_decomposition()

	! MPI domain decomposition
		call FACTOR (npu,lnpu) ! factors npu into powers of 2, 3 AND 5 AND stores the results into
		! the components of lnpu - from ffte's fft235.f
		
		! choosing adequate values for the number of subdivisions on y- AND z- directions
		npuy = (2**(lnpu(1)/2))*(3**(lnpu(2)/2))*(5**(lnpu(3)/2))
		npuz = npu/npuy
		
		if ( MOD(ny,npuy) > 0 ) then !! impossible to parallelize
			call MPI_BARRIER (MPI_COMM_WORLD, ierr)
			if (id == 0) print *,'bad decomposition (ny=',ny, ', npuy=',npuy,')! choose another number of processes or change the number of points on y!'
			call MPI_FINALIZE (ierr)
			stop
		endif
		
		if ( MOD(nz,npuz) > 0 ) then !! impossible to parallelize
			call MPI_BARRIER (MPI_COMM_WORLD, ierr)
			if (id == 0) print *,'bad decomposition (nz=',nz, ', npuz=',npuz,')! choose another number of processes or change the number of points on z!'
			call MPI_FINALIZE (ierr)
			stop
		endif
		
		if (nx<npuy .or. nx<npuz) then 
			call MPI_BARRIER (MPI_COMM_WORLD, ierr)
			if (id == 0) print *,'bad decomposition (nx=',nx, ', npuy=', npuy, ', npuz=',npuz,')! nx must be greater than both npuy and npuz!'
			call MPI_FINALIZE (ierr)
			stop
		endif
	
	nyp = ny/npuy
	nzp = nz/npuz
	
	idy = mod(id,npuy)
	idz = id/npuy
	
	if (id == 0) then
		write(*,*) 'npuy=', npuy, 'nyp =', nyp
		write(*,*) 'npuz=', npuz, 'nzp =', nzp 
	endif
	
	call MPI_COMM_SPLIT (MPI_COMM_WORLD, idz, 0, icommy, ierr)
	call MPI_COMM_SPLIT (MPI_COMM_WORLD, idy, 0, icommz, ierr)	
	
	print *,'process ', id, ' of ', npu, ' is alive in node: ', trim(procname)
	
	call MPI_BARRIER (MPI_COMM_WORLD, ierr)


	! MPI mapping
	! north
	node%n = mod (idz + 1, npuz)*npuy + idy
	! northeast
	node%ne = mod (idz + 1, npuz)*npuy + mod (idy + 1, npuy)
	! east
	node%e = idz*npuy + mod (idy + 1, npuy)
	! southeast
	node%se = mod (idz + npuz - 1, npuz)*npuy + mod (idy + 1, npuy)
	! south
	node%s = mod (idz + npuz - 1, npuz)*npuy + idy
	! southwest
	node%sw = mod (idz + npuz - 1, npuz)*npuy + mod (idy + npuy - 1, npuy)
	! west
	node%w = idz*npuy + mod (idy + npuy - 1, npuy)
	! northwest
	node%nw = mod (idz + 1, npuz)*npuy + mod (idy + npuy - 1, npuy)

	
	end subroutine mpi_decomposition
	!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 



	!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	subroutine allocate_fields()
	
	
	allocate(kx(nx),ky(ny),kz(nz))
	allocate(x(nx),y(ny),z(nz))
	allocate (yp (nyp), kyp (nyp), zp(nzp), kzp(nzp))
	allocate(k_2p(nx,nyp,nzp),u(nx,nyp,nzp),ua(nx,nyp,nzp),v(nx,nyp,nzp),&
	va(nx,nyp,nzp),w(nx,nyp,nzp),wa(nx,nyp,nzp),unew(nx,nyp,nzp),vnew(nx,nyp,nzp),&
	wnew(nx,nyp,nzp),pressure(nx,nyp,nzp),vort_1(nx,nyp,nzp),vort_2(nx,nyp,nzp),&
	vort_3(nx,nyp,nzp),vort_mag(nx,nyp,nzp),qcrit(nx,nyp,nzp),ab_u(nx,nyp,nzp,3),&
	ab_v(nx,nyp,nzp,3),ab_w(nx,nyp,nzp,3),filter(nx,nyp,nzp),aux(nx,nyp,nzp),&
    nueff(nx,nyp,nzp),fc(nx,nyp,nzp))
	allocate(xsave(0:nx),ypsave(0:nyp),zpsave(0:nzp))
	allocate(usave(0:nx,0:nyp,0:nzp),vsave(0:nx,0:nyp,0:nzp),wsave(0:nx,0:nyp,0:nzp),&
	psave(0:nx,0:nyp,0:nzp),vort_1_save(0:nx,0:nyp,0:nzp),vort_2_save(0:nx,0:nyp,0:nzp),&
	vort_3_save(0:nx,0:nyp,0:nzp),qcrit_save(0:nx,0:nyp,0:nzp))

	end subroutine allocate_fields
	!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  



	!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
	subroutine define_constants()


		!Imaginary unit
		imag=(0.0_rp,1.0_rp)
	
	
		!Parameters of HALE-RK7 scheme
		alphark(1) = 0.0_rp
		alphark(2) = -0.6479007459340_rp
		alphark(3) = -2.7047608632040_rp
		alphark(4) = -0.4600805501180_rp
		alphark(5) = -0.5005817877850_rp
		alphark(6) = -1.9065322559130_rp
		alphark(7) = -1.450_rp
		betark(1)  = 0.1173221468690_rp
		betark(2)  = 0.5032702621270_rp
		betark(3)  = 0.2336632816580_rp
		betark(4)  = 0.2834196346250_rp
		betark(5)  = 0.5403674140230_rp
		betark(6)  = 0.3714994146200_rp
		betark(7)  = 0.1366700993850_rp
		crk(1)     = 0.0_rp
		crk(2)     = 0.1173221468690_rp
		crk(3)     = 0.2945232307580_rp
		crk(4)     = 0.3056586221310_rp
		crk(5)     = 0.5828641484030_rp
		crk(6)     = 0.8586642735990_rp
		crk(7)     = 0.8686642735990_rp
	
	end subroutine define_constants
	!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  



	!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
	subroutine evaluate_dimensions()
		 
	 pi=2.0_rp*acos(0.0_rp)

	 dy=Ly/(ny)
	 dz=Lz/(nz)
	
	 !if (domain_set == 3) then !three-dimensional simulation
	 dx=Lx/(nx)
	
	 !else !two-dimensional simulation
	 !	dx = dy
	 !	Lx = dx
	 !endif
	 
	Re = 1_rp/mu       !Reynolds number
	
	if (id == 0) write(*,*)'Reynolds number = ',Re

		
	end subroutine evaluate_dimensions
	!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  



	!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
	subroutine initialize_variables()
	

	!Initialize fields
	
	!Sensitive data (should not be initialized if restart = true)
	if (restart_flag .eqv. .FALSE.) then
		u          = 0.0_rp
		v          = 0.0_rp
		w          = 0.0_rp
		t          = 0.0_rp
		tc_restart = 1
	endif
	
	!General data (should always be initialized)
	x      = 0.0_rp
	y      = 0.0_rp
	z      = 0.0_rp
	

	!Initialize FFTE
	call PZFFT3DV (u, unew, nx, ny, nz, icommy, icommz, npuy, npuz, 0)
	call PZFFT3DV (v, vnew, nx, ny, nz, icommy, icommz, npuy, npuz, 0)
	call PZFFT3DV (w, wnew, nx, ny, nz, icommy, icommz, npuy, npuz, 0)
	call PZFFT3DV (pressure, pressure, nx, ny, nz, icommy, icommz, npuy, npuz, 0)

	end subroutine initialize_variables
	!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
	
	
	
	!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
	subroutine create_mesh()

	!x-direction
	do i=1,nx
	   x(i) = i*dx
   	   xsave(i) = x(i)
	enddo
	xsave(0) = x(1) - dx

	
	!y-direction
	do j=1,ny
	   y(j)=j*dy
	enddo
	
	!z-direction
	do k=1,nz
	   z(k)=k*dz
	enddo
	
	end subroutine create_mesh
	!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  



	!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
	subroutine evaluate_wavenumbers()

	! wavenumber vectors
	do i = 1,nx/2+1
		kx(i) = 2.0_rp*pi*(i-1)/Lx
	enddo
	do i = nx/2+2,nx
		kx(i) = 2.0_rp*pi*(i-1-nx)/Lx
	enddo
	

	do j = 1, nyp
		yp (j) = (idy*nyp+j)*dy
		ypsave (j) = yp (j)
		if (idy*nyp+j <= ny/2+1) then
			kyp(j) = 2.0_rp*pi*(idy*nyp+j-1)/Ly
		else
			kyp(j) = 2.0_rp*pi*(idy*nyp+j-1-ny)/Ly
		endif
	enddo
	ypsave(0) = yp(1)-dy
	
	
	do k = 1, nzp
		zp (k) = (idz*nzp+k)*dz
		zpsave (k) = zp(k)
		if (idz*nzp+k <= nz/2+1) then
			kzp(k) = 2.0_rp*pi*(idz*nzp+k-1)/Lz
		else
			kzp(k) = 2.0_rp*pi*(idz*nzp+k-1-nz)/Lz
		endif
	enddo
	zpsave(0) = zp(1)-dz


	!evaluate k**2
	small=1E-30
	do k = 1, nzp
		do j = 1, nyp
			do i = 1, nx
				k_2p(i,j,k) = kx(i)**2 + kyp(j)**2 + kzp(k)**2 + small
			enddo
		enddo
	enddo


	k0 = 2.0_rp*pi*min(1.0_rp/Lx,1.0_rp/Ly,1.0_rp/Lz)
	dk = k0
	kmax = 2.0_rp*pi*max(nx/Lx, ny/Ly, nz/Lz) !with or without 2.0_rp?????????
	n3 = nx*ny*nz
	n3p = nx*nyp*nzp


	end subroutine evaluate_wavenumbers
	!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
	
	
	
	!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
	subroutine evaluate_timestep()

	if (temp_input_set == 2) then
		if (domain_set == 2) then
	 		dt1 = min((dy/max_v0),(dz/max_w0))
	 		dt2 = (1.0_rp/(2.0_rp*(mu/max_rho0)))*((1.0_rp/dy**2)+&
	 		(1.0_rp/dz**2))**(-1)
	 		dt = CFL*(min(dt1,dt2))
	 	else
	 		dt1 = min((dx/max_u0), (dy/max_v0), (dz/max_w0))
			dt2 = (1.0_rp/(2.0_rp*(mu/max_rho0)))*((1.0_rp/dx**2)+&
			(1.0_rp/dy**2)+(1.0_rp/dz**2))**(-1)
	 		dt = CFL*(min(dt1,dt2))
	 	endif
	 endif
	 
	 if (id==0) then
	 write(*,*)'dt = ',dt 
	 write(*,*)''
	 write(*,*)'INITIALIZING CALCULATIONS'
	 write(*,*)'#############################################'
	 write(*,*)''
	 write(*,*)''
	 write(*,*)'INITIALIZING INCOMPRESSIBLE SOLVER'
	 write(*,*)'#############################################'
	 write(*,*)''
	 endif
	 
	end subroutine evaluate_timestep
	!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  


	!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
	subroutine evaluate_filter()
	
	if (filter_set == 1) then !no filter
		filter = 1.0_rp

		else if (filter_set == 2) then !Sharpened Raised Cosine (SRC)
			 do k=1,nzp
				do j=1,nyp
					do i=1,nx
						filter(i,j,k) = (1.0_rp/8.0_rp)*(1.0_rp + cos(Lx*kx(i)/nx))*&
						(1.0_rp + cos(Ly*kyp(j)/ny))*(1.0_rp + cos(Lz*kzp(k)/nz))
					enddo
				enddo
			 enddo
			
			 filter = (filter**4)*(35.0_rp -84.0_rp*filter + 70.0_rp*&
			 (filter**2) - 20.0_rp*(filter**3))
	
	
		else if (filter_set == 3) then !Raised Cosine (RC)
			 do k=1,nzp
				do j=1,nyp
					do i=1,nx
						filter(i,j,k) = (1.0_rp/8.0_rp)*(1.0_rp + cos(Lx*kx(i)/nx))*&
						(1.0_rp + cos(Ly*kyp(j)/ny))*(1.0_rp + cos(Lz*kzp(k)/nz))
					enddo
				end do
			end do
			

		else if (filter_set == 4) then !Lanczos' filter
			 do k=1,nzp
				do j=1,nyp
					do i=1,nx
						filter(i,j,k) = (sin(Lx*kx(i)/nx)/(Lx*kx(i)/nx))*&
						(sin(Ly*kyp(j)/ny)/(Ly*kyp(j)/ny))*(sin(Lz*kz(k)/nz)/(Lz*kzp(k)/nz))
					enddo
				end do
			end do
			filter(1:nx,1,1)=1.0_rp
			filter(1,1:nyp,1)=1.0_rp
			filter(1,1,1:nzp)=1.0_rp
			
			
		else if (filter_set == 5) then !Exponential filter
			 do k=1,nzp
				do j=1,nyp
					do i=1,nx
						filter(i,j,k) = (exp(-1.0_rp*(Lx*kx(i)/nx)**2))*&
						(exp(-1.0_rp*(Ly*kyp(j)/ny)**2))
					enddo
				end do
			end do
			filter(1:nx,1,1)=filter(1:nx,nyp,nzp)
			filter(1,1:npy,1)=filter(nx,1:nyp,nzp)
			filter(1,1,1:npz)=filter(nx,npy,1:nzp)
	
	endif



	end subroutine evaluate_filter
	!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  


end module constants
