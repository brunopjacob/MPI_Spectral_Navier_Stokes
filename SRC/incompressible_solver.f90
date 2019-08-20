module incompressible_solver

use mpi

use declare_variables

public::evaluate_initial_conditions_incompressible,&
evaluate_rhs_incompressible,filter_field,projection,random_generator,&
screen_info,temporal_integration_incompressible,grid_to_print,&
evaluate_spectrum

contains


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!Subroutine: evaluate initial conditions in Fourier space (R -> F) (incompressible)
 subroutine evaluate_initial_conditions_incompressible(ic_set,domain_set,x,yp,zp,Lx,Ly,&
			Lz,nx,ny,nz,nyp,nzp,icommy,icommz,npuy,npuz,mu,amplitude_perturb,u,v,w,max_u0,max_v0,max_w0)
 
 implicit none
 integer,parameter:: rp=selected_real_kind(8,200),sp=selected_int_kind(8)

 !intent in
 integer(sp)::nx,ny,nz,nyp,nzp,ic_set,domain_set,icommy,icommz,npuy,npuz,ierr
 real(rp),dimension(nx)::x
 real(rp),dimension(nyp)::yp
 real(rp),dimension(nzp)::zp
 real(rp)::Lx,Ly,Lz,t,amplitude_perturb,mu
 
 !intent in/out
 integer(sp)::i,j,k
 complex(rp),dimension(nx,nyp,nzp)::ue,ve,we,uf,vf,wf,aux
 real(rp),dimension(nx,nyp,nzp)::randomsetu,randomsetv,randomsetw
 real(rp)::d_cp,r_cp,s_cp,rc
 
 !intent out
 complex(rp),dimension(nx,nyp,nzp)::u,v,w
 real(rp)::max_u0,max_v0,max_w0
 
 
 !Initialize FFTE
 call PZFFT3DV (ue, uf, nx, ny, nz, icommy, icommz, npuy, npuz, 0)
 call PZFFT3DV (ve, vf, nx, ny, nz, icommy, icommz, npuy, npuz, 0)
 call PZFFT3DV (we, wf, nx, ny, nz, icommy, icommz, npuy, npuz, 0)

 t=0.0_rp

 !Evaluate initial u,v,w,rho fields

 if (ic_set==1) then   !temporal jet
	 do k=1,nzp
		do j=1,nyp
			do i=1,nx
	    		u(i,j,k)    =  0.0_rp 
	    		v(i,j,k)    =  0.0_rp 
	    		w(i,j,k)    =  1.0_rp - 0.692_rp*exp(-0.69315_rp*( (yp(j)-Ly/2.0_rp )**2))
	    	enddo
		enddo
	 enddo
	
	elseif (ic_set==2) then !temporal wake
	 do k=1,nzp
		do j=1,nyp
			do i=1,nx
	    		u(i,j,k)    = 0.0_rp 
	    		v(i,j,k)    = 0.0_rp 
	    		w(i,j,k)    = -( (1.0_rp)/2.0_rp + (1.0_rp)/2.0_rp*tanh( (-abs( yp(j) - Ly/2.0_rp ) +&
	    		              0.5_rp*(Lz/16.0_rp) ) / (2.0_rp * (Lz/16.0_rp)/20.0_rp) ) ) + (1.0_rp)     	
			enddo
	    enddo
	 enddo
	 
	 elseif (ic_set==3) then !circular jet
	 d_cp     = 2.0_rp       !diameter of circular profile
	 r_cp     = d_cp/2.0_rp  !radius of circular profile
	 s_cp     = d_cp/Lx      !slope of circular profile 
	 u = 0.0_rp
	 v = 0.0_rp
 	 do k=1,nzp
		do j=1,nyp
			do i=1,nx
    	        rc          = sqrt((x(i)-Lx/2.0_rp)**2+(yp(j)-Ly/2.0_rp)**2) + 1E-16 !circular jet and wake r-coordinate system
    	        
	    		if (rc <= (r_cp - s_cp)) then  
					w(i,j,k) = 1.0_rp
	    		elseif (rc >= (r_cp - s_cp)) then
					w(i,j,k) = 0.0_rp
	    		else  
					w(i,j,k) = 0.5_rp*(1.0_rp - tanh((rc-r_cp)/(2.0_rp*s_cp)))
				endif
	    		
			enddo
	    enddo
	 enddo
	 
endif

 
 !Impose random perturbation
 call random_generator(nx,ny,nz,nyp,nzp,randomsetu,randomsetv,randomsetw)
 
 if (domain_set==3) then
	 u      = real(u)   + randomsetu   *amplitude_perturb
 endif
 
 v      = real(v)   + randomsetv   *amplitude_perturb
 w      = real(w)   + randomsetw   *amplitude_perturb



 !Defining variables in Euclidean space
 ue      = u  	  !ue    = u in Euclidean space
 ve      = v  	  !ve    = v in Euclidean space
 we      = w  	  !we    = w in Euclidean space


 !Defining variables in Fourier space
 !uf      = ue       !uf    = u in Fourier space
 !vf      = ve       !vf    = v in Fourier space
 !wf      = we       !wf    = w in Fourier space


 !Inverse Fourier transform (R-> F) 
 aux = ue
 call PZFFT3DV (aux, uf, nx, ny, nz, icommy, icommz, npuy, npuz, -1)
 aux = ve
 call PZFFT3DV (aux, vf, nx, ny, nz, icommy, icommz, npuy, npuz, -1)
 aux = we
 call PZFFT3DV (aux, wf, nx, ny, nz, icommy, icommz, npuy, npuz, -1)


 u      = uf
 v      = vf
 w      = wf
  
  
 call MPI_ALLREDUCE(maxval(dabs(real(ue))), max_u0, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr )
 call MPI_ALLREDUCE(maxval(dabs(real(ve))), max_v0, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr )
 call MPI_ALLREDUCE(maxval(dabs(real(we))), max_w0, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr )
  
 
 call MPI_BARRIER (MPI_COMM_WORLD, ierr)
 
 end subroutine evaluate_initial_conditions_incompressible
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++




!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 !subroutine: evaluate rhs in Fourier space (R -> F) (Incompressible solver)
 subroutine evaluate_rhs_incompressible(nx,ny,nz,nyp,nzp,kx,kyp,kzp,k_2p,icommy,icommz,npuy,npuz,mu,u,v,w,imag,&
 rhsx,rhsy,rhsz,ns_x_1,ns_y_1,ns_z_1,vort_1,vort_2,vort_3,vort_mag,qcrit,vort_flag,tc)
 
 
 implicit none
 integer,parameter:: rp=selected_real_kind(8,200),sp=selected_int_kind(8)

 !intent in
 integer(sp)::nx,ny,nz,nyp,nzp,tc,vort_flag,icommy,icommz,npuy,npuz
 complex(rp),dimension(nx)::kx
 complex(rp),dimension(nyp)::kyp
 complex(rp),dimension(nzp)::kzp
 complex(rp),dimension(nx,nyp,nzp)::k_2p,u,v,w
 complex(rp)::imag
 real(rp)::mu

 !intent in/out
 integer(sp)::i,j,k
 complex(rp),dimension(nx,nyp,nzp)::uf,vf,wf,ue,ve,we,&
 dudxf,dudyf,dudzf,dvdxf,dvdyf,dvdzf,dwdxf,dwdyf,dwdzf,&
 dudxe,dudye,dudze,dvdxe,dvdye,dvdze,dwdxe,dwdye,dwdze,&
 d2udxdxf,d2udydyf,d2udzdzf,d2udxdyf,d2udxdzf,d2udydzf,&
 d2vdxdxf,d2vdydyf,d2vdzdzf,d2vdxdyf,d2vdxdzf,d2vdydzf,d2wdxdxf,d2wdydyf,d2wdzdzf,&
 d2wdxdyf,d2wdxdzf,d2wdydzf,d2udxdxe,d2udydye,d2udzdze,&
 d2udxdye,d2udxdze,d2udydze,d2vdxdxe,d2vdydye,d2vdzdze,d2vdxdye,d2vdxdze,d2vdydze,&
 d2wdxdxe,d2wdydye,d2wdzdze,d2wdxdye,d2wdxdze,d2wdydze,&
 ns_x_1,ns_x_2,ns_y_1,ns_y_2,ns_z_1,ns_z_2,&
 vort_1,vort_2,vort_3,vort_mag,&
 s_11,s_12,s_13,s_21,s_22,s_23,s_31,s_32,s_33,w_11,w_12,w_13,w_21,w_22,w_23,w_31,w_32,w_33,&
 qs,qw,qcrit,rcrit,rs,sigma_stret,aux,&
 nueff_s11,nueff_s12,nueff_s13,nueff_s21,nueff_s22,nueff_s23,&
 nueff_s31,nueff_s32,nueff_s33,ddx_nueff_s11,ddy_nueff_s12,ddz_nueff_s13, ddx_nueff_s21,&
 ddy_nueff_s22,ddz_nueff_s23,ddx_nueff_s31,ddy_nueff_s32,ddz_nueff_s33
 
 real(rp)::kappa,beta
 real(rp)::a1,b1,kmax1,l,cs,l1,l2
 complex(rp),dimension(nx,nyp,nzp)::ns !,nueff


!dynamic Smagorinsky model
 complex(rp),allocatable, dimension(:,:,:):: ufte,vfte,wfte,&
 uuf,uvf,uwf,vvf,vwf,wwf,&
 uufte,vvfte,wwfte,uvfte,uwfte,vwfte,uuftf,vvftf,wwftf,uvftf,uwftf,vwftf,&
 uftf,vftf,wftf,&
 aux1_M_11,aux1_M_12,aux1_M_13,aux1_M_21,aux1_M_22,aux1_M_23,aux1_M_31,aux1_M_32,aux1_M_33,aux2_M_11,aux2_M_12,aux2_M_13,aux2_M_21,aux2_M_22,aux2_M_23,aux2_M_31,aux2_M_32,aux2_M_33,&
 L_11,L_12,L_13,L_21,L_22,L_23,L_31,L_32,L_33,M_11,M_12,M_13,M_21,M_22,M_23,M_31,M_32,M_33,&
 s_11ft,s_12ft,s_13ft,s_21ft,s_22ft,s_23ft,s_31ft,s_32ft,s_33ft,nsft,MM,ML,csd

 
 !intent out
 complex(rp),dimension(nx,nyp,nzp)::rhsx,rhsy,rhsz  
  
 !initialize FFTE
 if (tc==0)then
 
 	call PZFFT3DV (uf, ue, nx, ny, nz, icommy, icommz, npuy, npuz, 0)
	call PZFFT3DV (vf, ve, nx, ny, nz, icommy, icommz, npuy, npuz, 0)
	call PZFFT3DV (wf, we, nx, ny, nz, icommy, icommz, npuy, npuz, 0)
	call PZFFT3DV (dudxf, dudxe, nx, ny, nz, icommy, icommz, npuy, npuz, 0)
	call PZFFT3DV (dudyf, dudye, nx, ny, nz, icommy, icommz, npuy, npuz, 0)
	call PZFFT3DV (dudzf, dudze, nx, ny, nz, icommy, icommz, npuy, npuz, 0)
	call PZFFT3DV (dvdxf, dvdxe, nx, ny, nz, icommy, icommz, npuy, npuz, 0)
	call PZFFT3DV (dvdyf, dvdye, nx, ny, nz, icommy, icommz, npuy, npuz, 0)
	call PZFFT3DV (dvdzf, dvdze, nx, ny, nz, icommy, icommz, npuy, npuz, 0)
	call PZFFT3DV (dwdxf, dwdxe, nx, ny, nz, icommy, icommz, npuy, npuz, 0)
	call PZFFT3DV (dwdyf, dwdye, nx, ny, nz, icommy, icommz, npuy, npuz, 0)
	call PZFFT3DV (dwdzf, dwdze, nx, ny, nz, icommy, icommz, npuy, npuz, 0)
	call PZFFT3DV (d2udxdxf, d2udxdxe, nx, ny, nz, icommy, icommz, npuy, npuz, 0)
	call PZFFT3DV (d2udydyf, d2udydye, nx, ny, nz, icommy, icommz, npuy, npuz, 0)
	call PZFFT3DV (d2udzdzf, d2udzdze, nx, ny, nz, icommy, icommz, npuy, npuz, 0)
	call PZFFT3DV (d2udxdyf, d2udxdye, nx, ny, nz, icommy, icommz, npuy, npuz, 0)
	call PZFFT3DV (d2udxdzf, d2udxdze, nx, ny, nz, icommy, icommz, npuy, npuz, 0)
	call PZFFT3DV (d2udydzf, d2udydze, nx, ny, nz, icommy, icommz, npuy, npuz, 0)
	call PZFFT3DV (d2vdxdxf, d2vdxdxe, nx, ny, nz, icommy, icommz, npuy, npuz, 0)
	call PZFFT3DV (d2vdydyf, d2vdydye, nx, ny, nz, icommy, icommz, npuy, npuz, 0)
	call PZFFT3DV (d2vdzdzf, d2vdzdze, nx, ny, nz, icommy, icommz, npuy, npuz, 0)
	call PZFFT3DV (d2vdxdyf, d2vdxdye, nx, ny, nz, icommy, icommz, npuy, npuz, 0)
	call PZFFT3DV (d2vdxdzf, d2vdxdze, nx, ny, nz, icommy, icommz, npuy, npuz, 0)
	call PZFFT3DV (d2vdydzf, d2vdydze, nx, ny, nz, icommy, icommz, npuy, npuz, 0)
	call PZFFT3DV (d2wdxdxf, d2wdxdxe, nx, ny, nz, icommy, icommz, npuy, npuz, 0)
	call PZFFT3DV (d2wdydyf, d2wdydye, nx, ny, nz, icommy, icommz, npuy, npuz, 0)
	call PZFFT3DV (d2wdzdzf, d2wdzdze, nx, ny, nz, icommy, icommz, npuy, npuz, 0)
	call PZFFT3DV (d2wdxdyf, d2wdxdye, nx, ny, nz, icommy, icommz, npuy, npuz, 0)
	call PZFFT3DV (d2wdxdzf, d2wdxdze, nx, ny, nz, icommy, icommz, npuy, npuz, 0)
	call PZFFT3DV (d2wdydzf, d2wdydze, nx, ny, nz, icommy, icommz, npuy, npuz, 0)
	call PZFFT3DV (ns_x_1, ns_x_1, nx, ny, nz, icommy, icommz, npuy, npuz, 0)
	call PZFFT3DV (ns_x_2, ns_x_2, nx, ny, nz, icommy, icommz, npuy, npuz, 0)
	call PZFFT3DV (ns_y_1, ns_y_1, nx, ny, nz, icommy, icommz, npuy, npuz, 0)
	call PZFFT3DV (ns_y_2, ns_y_2, nx, ny, nz, icommy, icommz, npuy, npuz, 0)
	call PZFFT3DV (ns_z_1, ns_z_1, nx, ny, nz, icommy, icommz, npuy, npuz, 0)
	call PZFFT3DV (ns_z_2, ns_z_2, nx, ny, nz, icommy, icommz, npuy, npuz, 0)
    call PZFFT3DV (nueff_s11, nueff_s11, nx, ny, nz, icommy, icommz, npuy, npuz, 0)
	call PZFFT3DV (nueff_s12, nueff_s12, nx, ny, nz, icommy, icommz, npuy, npuz, 0)
	call PZFFT3DV (nueff_s13, nueff_s13, nx, ny, nz, icommy, icommz, npuy, npuz, 0)
	call PZFFT3DV (nueff_s21, nueff_s21, nx, ny, nz, icommy, icommz, npuy, npuz, 0)
	call PZFFT3DV (nueff_s22, nueff_s22, nx, ny, nz, icommy, icommz, npuy, npuz, 0)
	call PZFFT3DV (nueff_s23, nueff_s23, nx, ny, nz, icommy, icommz, npuy, npuz, 0)
	call PZFFT3DV (nueff_s31, nueff_s31, nx, ny, nz, icommy, icommz, npuy, npuz, 0)
	call PZFFT3DV (nueff_s32, nueff_s32, nx, ny, nz, icommy, icommz, npuy, npuz, 0)
	call PZFFT3DV (nueff_s33, nueff_s33, nx, ny, nz, icommy, icommz, npuy, npuz, 0)
	call PZFFT3DV (ddx_nueff_s11, ddx_nueff_s11, nx, ny, nz, icommy, icommz, npuy, npuz, 0)
	call PZFFT3DV (ddy_nueff_s12, ddy_nueff_s12, nx, ny, nz, icommy, icommz, npuy, npuz, 0)
	call PZFFT3DV (ddz_nueff_s13, ddz_nueff_s13, nx, ny, nz, icommy, icommz, npuy, npuz, 0)
	call PZFFT3DV (ddx_nueff_s21, ddx_nueff_s21, nx, ny, nz, icommy, icommz, npuy, npuz, 0)
	call PZFFT3DV (ddy_nueff_s22, ddy_nueff_s22, nx, ny, nz, icommy, icommz, npuy, npuz, 0)
	call PZFFT3DV (ddz_nueff_s23, ddz_nueff_s23, nx, ny, nz, icommy, icommz, npuy, npuz, 0)
	call PZFFT3DV (ddx_nueff_s31, ddx_nueff_s31, nx, ny, nz, icommy, icommz, npuy, npuz, 0)
	call PZFFT3DV (ddy_nueff_s32, ddy_nueff_s32, nx, ny, nz, icommy, icommz, npuy, npuz, 0)
	call PZFFT3DV (ddz_nueff_s33, ddz_nueff_s33, nx, ny, nz, icommy, icommz, npuy, npuz, 0)
	
 endif

!Initialize variables
 vort_1   = 0.0_rp
 vort_2   = 0.0_rp
 vort_3   = 0.0_rp
 vort_mag = 0.0_rp
 qcrit    = 0.0_rp


!Defining variables in Fourier space
 uf       = u       !uf    = u in Fourier space
 vf       = v       !vf    = v in Fourier space
 wf       = w       !wf    = w in Fourier space

!Defining variables in Euclidean space
! ue       = u       !ue     = u in Euclidean space
! ve       = v       !ve     = v in Euclidean space
! we       = w       !we     = w in Euclidean space

 aux = uf
 call PZFFT3DV (aux, ue, nx, ny, nz, icommy, icommz, npuy, npuz, +1)
 aux = vf
 call PZFFT3DV (aux, ve, nx, ny, nz, icommy, icommz, npuy, npuz, +1)
 aux = wf
 call PZFFT3DV (aux, we, nx, ny, nz, icommy, icommz, npuy, npuz, +1)

 ue      = real(ue)
 ve      = real(ve)
 we      = real(we)
 

		
! 	!evaluate 1st order derivatives  ===================================
	
	do k=1,nzp
		do j=1,nyp
			do i=1,nx
				dudxf(i,j,k) = imag*kx(i)*uf(i,j,k)
				dudyf(i,j,k) = imag*kyp(j)*uf(i,j,k)
				dudzf(i,j,k) = imag*kzp(k)*uf(i,j,k)
				dvdxf(i,j,k) = imag*kx(i)*vf(i,j,k)
				dvdyf(i,j,k) = imag*kyp(j)*vf(i,j,k)
				dvdzf(i,j,k) = imag*kzp(k)*vf(i,j,k)
				dwdxf(i,j,k) = imag*kx(i)*wf(i,j,k)
				dwdyf(i,j,k) = imag*kyp(j)*wf(i,j,k)
				dwdzf(i,j,k) = imag*kzp(k)*wf(i,j,k)
			enddo
		enddo
	enddo


	aux = dudxf
	call PZFFT3DV (aux, dudxe, nx, ny, nz, icommy, icommz, npuy, npuz, +1)
	aux = dudyf
	call PZFFT3DV (aux, dudye, nx, ny, nz, icommy, icommz, npuy, npuz, +1)
	aux = dudzf
	call PZFFT3DV (aux, dudze, nx, ny, nz, icommy, icommz, npuy, npuz, +1)
	aux = dvdxf
	call PZFFT3DV (aux, dvdxe, nx, ny, nz, icommy, icommz, npuy, npuz, +1)
	aux = dvdyf
	call PZFFT3DV (aux, dvdye, nx, ny, nz, icommy, icommz, npuy, npuz, +1)
	aux = dvdzf
	call PZFFT3DV (aux, dvdze, nx, ny, nz, icommy, icommz, npuy, npuz, +1)
	aux = dwdxf
	call PZFFT3DV (aux, dwdxe, nx, ny, nz, icommy, icommz, npuy, npuz, +1)
	aux = dwdyf
	call PZFFT3DV (aux, dwdye, nx, ny, nz, icommy, icommz, npuy, npuz, +1)
	aux = dwdzf
	call PZFFT3DV (aux, dwdze, nx, ny, nz, icommy, icommz, npuy, npuz, +1)


	dudxe = real(dudxe)
	dudye = real(dudye)
	dudze = real(dudze)
	dvdxe = real(dvdxe)
	dvdye = real(dvdye)
	dvdze = real(dvdze)
	dwdxe = real(dwdxe)
	dwdye = real(dwdye)
	dwdze = real(dwdze)
	

	

    
    		
	!===================================================================
 	!Turbulence modeling
 	!===================================================================
 
	!DNS ===============================================================
	if (turbulence_modeling_set == 1) then
	nueff = mu
	fc = 1.0_rp
 	
 	
 	!Smagorinsky model =================================================
 	else if (turbulence_modeling_set == 2) then
		kmax1=kmax/2.0_rp
		
	 	l=(dx*dy*dz)**(1.0_rp/3.0_rp)
		cs=0.12_rp
		
		ns=sqrt(real(2.0_rp*(real(s_11)**2+real(s_22)**2+real(s_33)**2+2.0_rp*(real(s_12)**2+real(s_13)**2+real(s_23)**2))))
		nueff=mu+((cs*l)**2)*ns

 
 	!Dynamic Smagorinsky model =========================================
	else if (turbulence_modeling_set == 3) then		
		
		allocate(uuf(nx,nyp,nzp),uvf(nx,nyp,nzp),uwf(nx,nyp,nzp),vvf(nx,nyp,nzp),vwf(nx,nyp,nzp),wwf(nx,nyp,nzp))
		allocate(uufte(nx,nyp,nzp),vvfte(nx,nyp,nzp),wwfte(nx,nyp,nzp),uvfte(nx,nyp,nzp),uwfte(nx,nyp,nzp),vwfte(nx,nyp,nzp))
		allocate(ufte(nx,nyp,nzp),vfte(nx,nyp,nzp),wfte(nx,nyp,nzp))
		allocate(aux1_M_11(nx,nyp,nzp),aux1_M_12(nx,nyp,nzp),aux1_M_13(nx,nyp,nzp),aux1_M_21(nx,nyp,nzp),aux1_M_22(nx,nyp,nzp),aux1_M_23(nx,nyp,nzp),aux1_M_31(nx,nyp,nzp),aux1_M_32(nx,nyp,nzp),aux1_M_33(nx,nyp,nzp))
		allocate(aux2_M_11(nx,nyp,nzp),aux2_M_12(nx,nyp,nzp),aux2_M_13(nx,nyp,nzp),aux2_M_21(nx,nyp,nzp),aux2_M_22(nx,nyp,nzp),aux2_M_23(nx,nyp,nzp),aux2_M_31(nx,nyp,nzp),aux2_M_32(nx,nyp,nzp),aux2_M_33(nx,nyp,nzp))
		allocate(L_11(nx,nyp,nzp),L_12(nx,nyp,nzp),L_13(nx,nyp,nzp),L_21(nx,nyp,nzp),L_22(nx,nyp,nzp),L_23(nx,nyp,nzp),L_31(nx,nyp,nzp),L_32(nx,nyp,nzp),L_33(nx,nyp,nzp))
		allocate(M_11(nx,nyp,nzp),M_12(nx,nyp,nzp),M_13(nx,nyp,nzp),M_21(nx,nyp,nzp),M_22(nx,nyp,nzp),M_23(nx,nyp,nzp),M_31(nx,nyp,nzp),M_32(nx,nyp,nzp),M_33(nx,nyp,nzp))
		allocate(s_11ft(nx,nyp,nzp),s_12ft(nx,nyp,nzp),s_13ft(nx,nyp,nzp),s_21ft(nx,nyp,nzp),s_22ft(nx,nyp,nzp),s_23ft(nx,nyp,nzp),s_31ft(nx,nyp,nzp),s_32ft(nx,nyp,nzp),s_33ft(nx,nyp,nzp))
		allocate(nsft(nx,nyp,nzp),MM(nx,nyp,nzp),ML(nx,nyp,nzp))
		allocate(csd(nx,nyp,nzp))
		
		call PZFFT3DV(uuf, uuf, nx, ny, nz, icommy, icommz, npuy, npuz, 0)
		call PZFFT3DV(vvf, vvf, nx, ny, nz, icommy, icommz, npuy, npuz, 0)
		call PZFFT3DV(wwf, wwf, nx, ny, nz, icommy, icommz, npuy, npuz, 0)
		call PZFFT3DV(uvf, uvf, nx, ny, nz, icommy, icommz, npuy, npuz, 0)
		call PZFFT3DV(uwf, uwf, nx, ny, nz, icommy, icommz, npuy, npuz, 0)
		call PZFFT3DV(vwf, vwf, nx, ny, nz, icommy, icommz, npuy, npuz, 0)
		call PZFFT3DV(uufte, uufte, nx, ny, nz, icommy, icommz, npuy, npuz, 0)
		call PZFFT3DV(vvfte, vvfte, nx, ny, nz, icommy, icommz, npuy, npuz, 0)
		call PZFFT3DV(wwfte, wwfte, nx, ny, nz, icommy, icommz, npuy, npuz, 0)
		call PZFFT3DV(uvfte, uvfte, nx, ny, nz, icommy, icommz, npuy, npuz, 0)
		call PZFFT3DV(uwfte, uwfte, nx, ny, nz, icommy, icommz, npuy, npuz, 0)
		call PZFFT3DV(vwfte, vwfte, nx, ny, nz, icommy, icommz, npuy, npuz, 0)
		call PZFFT3DV(ufte, ufte, nx, ny, nz, icommy, icommz, npuy, npuz, 0)
		call PZFFT3DV(vfte, vfte, nx, ny, nz, icommy, icommz, npuy, npuz, 0)
		call PZFFT3DV(wfte, wfte, nx, ny, nz, icommy, icommz, npuy, npuz, 0)
		call PZFFT3DV(aux2_M_11, aux2_M_11, nx, ny, nz, icommy, icommz, npuy, npuz, 0)
		call PZFFT3DV(aux2_M_12, aux2_M_12, nx, ny, nz, icommy, icommz, npuy, npuz, 0)
		call PZFFT3DV(aux2_M_13, aux2_M_13, nx, ny, nz, icommy, icommz, npuy, npuz, 0)
		call PZFFT3DV(aux2_M_21, aux2_M_21, nx, ny, nz, icommy, icommz, npuy, npuz, 0)
		call PZFFT3DV(aux2_M_22, aux2_M_22, nx, ny, nz, icommy, icommz, npuy, npuz, 0)
		call PZFFT3DV(aux2_M_23, aux2_M_23, nx, ny, nz, icommy, icommz, npuy, npuz, 0)
		call PZFFT3DV(aux2_M_31, aux2_M_31, nx, ny, nz, icommy, icommz, npuy, npuz, 0)
		call PZFFT3DV(aux2_M_32, aux2_M_32, nx, ny, nz, icommy, icommz, npuy, npuz, 0)
		call PZFFT3DV(aux2_M_33, aux2_M_33, nx, ny, nz, icommy, icommz, npuy, npuz, 0)


		!Compute filter matrix
		kmax1=kmax/2.0_rp
		a1=5.
		b1=4.
		do k=1,nzp
			do j =1,nyp
				do i=1,nx
					if (sqrt(real(k_2p(i,j,k)))<=kmax1) then
						fc(i,j,k)=1.0_rp
					else
						fc(i,j,k)=0.0_rp!exp(-a1*(sqrt(real(k_2(i,j,k)))-kmax1)**b1)
					endif
				enddo
			enddo
		enddo
		
        
		!Compute products of non-filtered velocities
		uuf = ue*ue
		uvf = ue*ve
		uwf = ue*we
		vvf = ve*ve
		vwf = ve*we
		wwf = we*we
		

   	    call PZFFT3DV (uuf, uuf, nx, ny, nz, icommy, icommz, npuy, npuz, -1)
   	    call PZFFT3DV (vvf, vvf, nx, ny, nz, icommy, icommz, npuy, npuz, -1)
   	    call PZFFT3DV (wwf, wwf, nx, ny, nz, icommy, icommz, npuy, npuz, -1)
   	    call PZFFT3DV (uvf, uvf, nx, ny, nz, icommy, icommz, npuy, npuz, -1)
   	    call PZFFT3DV (uwf, uwf, nx, ny, nz, icommy, icommz, npuy, npuz, -1)
   	    call PZFFT3DV (vwf, vwf, nx, ny, nz, icommy, icommz, npuy, npuz, -1)
		
		uufte = uuf*fc
		uvfte = uvf*fc
		uwfte = uwf*fc
		vvfte = vvf*fc	
		vwfte = vwf*fc	
		wwfte = wwf*fc	

		deallocate(uuf,uvf,uwf,vvf,vwf,wwf)

   	    call PZFFT3DV (uufte, uufte, nx, ny, nz, icommy, icommz, npuy, npuz, +1)
   	    call PZFFT3DV (vvfte, vvfte, nx, ny, nz, icommy, icommz, npuy, npuz, +1)
   	    call PZFFT3DV (wwfte, wwfte, nx, ny, nz, icommy, icommz, npuy, npuz, +1)
   	    call PZFFT3DV (uvfte, uvfte, nx, ny, nz, icommy, icommz, npuy, npuz, +1)
   	    call PZFFT3DV (uwfte, uwfte, nx, ny, nz, icommy, icommz, npuy, npuz, +1)
   	    call PZFFT3DV (vwfte, vwfte, nx, ny, nz, icommy, icommz, npuy, npuz, +1)

		
		!Compute filtered velocities
		ufte = uf*fc
		vfte = vf*fc
		wfte = wf*fc
			        
        call PZFFT3DV (ufte, ufte, nx, ny, nz, icommy, icommz, npuy, npuz, +1)
   	    call PZFFT3DV (vfte, vfte, nx, ny, nz, icommy, icommz, npuy, npuz, +1)
   	    call PZFFT3DV (wfte, wfte, nx, ny, nz, icommy, icommz, npuy, npuz, +1)
		
				
		!Leonard tensor (components filtered)
		L_11 = uufte - ufte*ufte
		L_12 = uvfte - ufte*vfte
		L_13 = uwfte - ufte*wfte
		L_21 = L_12
		L_22 = vvfte - vfte*vfte
		L_23 = vwfte - vfte*wfte
		L_31 = L_13
		L_32 = L_23
		L_33 = wwfte - wfte*wfte


		deallocate(ufte,vfte,wfte)
		deallocate(uufte,vvfte,wwfte,uvfte,uwfte,vwfte)


		!Shear stress tensor of filtered velocities
		s_11ft = s_11
		s_12ft = s_12
		s_13ft = s_13
		s_21ft = s_21
		s_22ft = s_22
		s_23ft = s_23
		s_31ft = s_31
		s_32ft = s_32
		s_33ft = s_33
		
        
        call PZFFT3DV (s_11ft, s_11ft, nx, ny, nz, icommy, icommz, npuy, npuz, -1)
   	    call PZFFT3DV (s_12ft, s_12ft, nx, ny, nz, icommy, icommz, npuy, npuz, -1)
   	    call PZFFT3DV (s_13ft, s_13ft, nx, ny, nz, icommy, icommz, npuy, npuz, -1)
   	    call PZFFT3DV (s_21ft, s_21ft, nx, ny, nz, icommy, icommz, npuy, npuz, -1)
   	    call PZFFT3DV (s_22ft, s_22ft, nx, ny, nz, icommy, icommz, npuy, npuz, -1)
   	    call PZFFT3DV (s_23ft, s_23ft, nx, ny, nz, icommy, icommz, npuy, npuz, -1)
   	    call PZFFT3DV (s_31ft, s_31ft, nx, ny, nz, icommy, icommz, npuy, npuz, -1)
   	    call PZFFT3DV (s_32ft, s_32ft, nx, ny, nz, icommy, icommz, npuy, npuz, -1)
   	    call PZFFT3DV (s_33ft, s_33ft, nx, ny, nz, icommy, icommz, npuy, npuz, -1)


		s_11ft = s_11ft*fc
		s_12ft = s_12ft*fc 
		s_13ft = s_13ft*fc
		s_21ft = s_21ft*fc
		s_22ft = s_22ft*fc
		s_23ft = s_23ft*fc
		s_31ft = s_31ft*fc
		s_32ft = s_32ft*fc
		s_33ft = s_33ft*fc


        call PZFFT3DV (s_11ft, s_11ft, nx, ny, nz, icommy, icommz, npuy, npuz, +1)
   	    call PZFFT3DV (s_12ft, s_12ft, nx, ny, nz, icommy, icommz, npuy, npuz, +1)
   	    call PZFFT3DV (s_13ft, s_13ft, nx, ny, nz, icommy, icommz, npuy, npuz, +1)
   	    call PZFFT3DV (s_21ft, s_21ft, nx, ny, nz, icommy, icommz, npuy, npuz, +1)
   	    call PZFFT3DV (s_22ft, s_22ft, nx, ny, nz, icommy, icommz, npuy, npuz, +1)
   	    call PZFFT3DV (s_23ft, s_23ft, nx, ny, nz, icommy, icommz, npuy, npuz, +1)
   	    call PZFFT3DV (s_31ft, s_31ft, nx, ny, nz, icommy, icommz, npuy, npuz, +1)
   	    call PZFFT3DV (s_32ft, s_32ft, nx, ny, nz, icommy, icommz, npuy, npuz, +1)
   	    call PZFFT3DV (s_33ft, s_33ft, nx, ny, nz, icommy, icommz, npuy, npuz, +1)


		!Define characteristic scales of filters 1 and 2
		l2=(dx*dy*dz)**(1.0_rp/3.0_rp)		
		l1=2.0_rp*l2
		
		!Compute norms of S
		ns   = sqrt(real(2.0_rp*(real(s_11)**2+real(s_22)**2+real(s_33)**2+2.0_rp*(real(s_12)**2+real(s_13)**2+real(s_23)**2)))) !norm of s (non-filtered)		
		nsft = sqrt(real(2.0_rp*(real(s_11ft)**2+real(s_22ft)**2+real(s_33ft)**2+2.0_rp*(real(s_12ft)**2+real(s_13ft)**2+real(s_23ft)**2)))) !norm of s (filtered)		


		!Compute product of filtered variables, l2**2 * nsft * Sftij = aux1
		aux1_M_11 = (l1**2)*nsft*s_11ft
		aux1_M_12 = (l1**2)*nsft*s_12ft
		aux1_M_13 = (l1**2)*nsft*s_13ft		
		aux1_M_21 = (l1**2)*nsft*s_21ft	
		aux1_M_22 = (l1**2)*nsft*s_22ft	
		aux1_M_23 = (l1**2)*nsft*s_23ft	
		aux1_M_31 = (l1**2)*nsft*s_31ft	
		aux1_M_32 = (l1**2)*nsft*s_32ft	
		aux1_M_33 = (l1**2)*nsft*s_33ft
	
		deallocate(s_11ft,s_12ft,s_13ft,s_21ft,s_22ft,s_23ft,s_31ft,s_32ft,s_33ft)

		
		!Filter products l1**2 * ns * Sij = aux1
		aux2_M_11 = (l2**2)*ns*s_11  
		aux2_M_12 = (l2**2)*ns*s_12
		aux2_M_13 = (l2**2)*ns*s_13
		aux2_M_21 = (l2**2)*ns*s_21
		aux2_M_22 = (l2**2)*ns*s_22
		aux2_M_23 = (l2**2)*ns*s_23
		aux2_M_31 = (l2**2)*ns*s_31
		aux2_M_32 = (l2**2)*ns*s_32
		aux2_M_33 = (l2**2)*ns*s_33



        call PZFFT3DV (aux2_M_11, aux2_M_11, nx, ny, nz, icommy, icommz, npuy, npuz, -1)
   	    call PZFFT3DV (aux2_M_12, aux2_M_12, nx, ny, nz, icommy, icommz, npuy, npuz, -1)
   	    call PZFFT3DV (aux2_M_13, aux2_M_13, nx, ny, nz, icommy, icommz, npuy, npuz, -1)
   	    call PZFFT3DV (aux2_M_21, aux2_M_21, nx, ny, nz, icommy, icommz, npuy, npuz, -1)
   	    call PZFFT3DV (aux2_M_22, aux2_M_22, nx, ny, nz, icommy, icommz, npuy, npuz, -1)
   	    call PZFFT3DV (aux2_M_23, aux2_M_23, nx, ny, nz, icommy, icommz, npuy, npuz, -1)
   	    call PZFFT3DV (aux2_M_31, aux2_M_31, nx, ny, nz, icommy, icommz, npuy, npuz, -1)
   	    call PZFFT3DV (aux2_M_32, aux2_M_32, nx, ny, nz, icommy, icommz, npuy, npuz, -1)
   	    call PZFFT3DV (aux2_M_33, aux2_M_33, nx, ny, nz, icommy, icommz, npuy, npuz, -1)


		aux2_M_11 = aux2_M_11*fc
		aux2_M_12 = aux2_M_12*fc
		aux2_M_13 = aux2_M_13*fc
		aux2_M_21 = aux2_M_21*fc
		aux2_M_22 = aux2_M_22*fc
		aux2_M_23 = aux2_M_23*fc
		aux2_M_31 = aux2_M_31*fc
		aux2_M_32 = aux2_M_32*fc
		aux2_M_33 = aux2_M_33*fc

		
        call PZFFT3DV (aux2_M_11, aux2_M_11, nx, ny, nz, icommy, icommz, npuy, npuz, +1)
   	    call PZFFT3DV (aux2_M_12, aux2_M_12, nx, ny, nz, icommy, icommz, npuy, npuz, +1)
   	    call PZFFT3DV (aux2_M_13, aux2_M_13, nx, ny, nz, icommy, icommz, npuy, npuz, +1)
   	    call PZFFT3DV (aux2_M_21, aux2_M_21, nx, ny, nz, icommy, icommz, npuy, npuz, +1)
   	    call PZFFT3DV (aux2_M_22, aux2_M_22, nx, ny, nz, icommy, icommz, npuy, npuz, +1)
   	    call PZFFT3DV (aux2_M_23, aux2_M_23, nx, ny, nz, icommy, icommz, npuy, npuz, +1)
   	    call PZFFT3DV (aux2_M_31, aux2_M_31, nx, ny, nz, icommy, icommz, npuy, npuz, +1)
   	    call PZFFT3DV (aux2_M_32, aux2_M_32, nx, ny, nz, icommy, icommz, npuy, npuz, +1)
   	    call PZFFT3DV (aux2_M_33, aux2_M_33, nx, ny, nz, icommy, icommz, npuy, npuz, +1)


		!Compute Mij
		M_11 = aux1_M_11 - aux2_M_11
		M_12 = aux1_M_12 - aux2_M_12
		M_13 = aux1_M_13 - aux2_M_13
		M_21 = aux1_M_21 - aux2_M_21
		M_22 = aux1_M_22 - aux2_M_22
		M_23 = aux1_M_23 - aux2_M_23
		M_31 = aux1_M_31 - aux2_M_31
		M_32 = aux1_M_32 - aux2_M_32
		M_33 = aux1_M_33 - aux2_M_33

		deallocate (aux1_M_11,aux1_M_12,aux1_M_13,aux1_M_21,aux1_M_22,aux1_M_23,aux1_M_31,aux1_M_32,aux1_M_33)
		deallocate (aux2_M_11,aux2_M_12,aux2_M_13,aux2_M_21,aux2_M_22,aux2_M_23,aux2_M_31,aux2_M_32,aux2_M_33)

		MM=M_11**2+M_22**2+M_33**2+2.0_rp*(M_12**2+M_13**2+M_23**2)+1E-20
		ML=M_11*L_11+M_22*L_22+M_33*L_33+2.0_rp*(M_12*L_12+M_13*L_13+M_23*L_23)

		deallocate (L_11,L_12,L_13,L_21,L_22,L_23,L_31,L_32,L_33)
		deallocate (M_11,M_12,M_13,M_21,M_22,M_23,M_31,M_32,M_33)

		csd = 0.50_rp*ML/MM 
		nueff=mu+real(((csd*l2)**2)*ns)
	
		do k=1,nzp
			do j =1,nyp
				do i=1,nx
					if (real(nueff(i,j,k))<0.0_rp) then
						nueff(i,j,k)=0.0_rp
					endif
				enddo
			enddo
		enddo
	
  	endif

	
	nueff_s11 = nueff*s_11 
	nueff_s12 = nueff*s_12
	nueff_s13 = nueff*s_13
	nueff_s21 = nueff*s_21
	nueff_s22 = nueff*s_22
	nueff_s23 = nueff*s_23
	nueff_s31 = nueff*s_31
	nueff_s32 = nueff*s_32
	nueff_s33 = nueff*s_33


    call PZFFT3DV (nueff_s11, nueff_s11, nx, ny, nz, icommy, icommz, npuy, npuz, -1)
   	call PZFFT3DV (nueff_s12, nueff_s12, nx, ny, nz, icommy, icommz, npuy, npuz, -1)
   	call PZFFT3DV (nueff_s13, nueff_s13, nx, ny, nz, icommy, icommz, npuy, npuz, -1)
   	call PZFFT3DV (nueff_s21, nueff_s21, nx, ny, nz, icommy, icommz, npuy, npuz, -1)
   	call PZFFT3DV (nueff_s22, nueff_s22, nx, ny, nz, icommy, icommz, npuy, npuz, -1)
   	call PZFFT3DV (nueff_s23, nueff_s23, nx, ny, nz, icommy, icommz, npuy, npuz, -1)
   	call PZFFT3DV (nueff_s31, nueff_s31, nx, ny, nz, icommy, icommz, npuy, npuz, -1)
   	call PZFFT3DV (nueff_s32, nueff_s32, nx, ny, nz, icommy, icommz, npuy, npuz, -1)
   	call PZFFT3DV (nueff_s33, nueff_s33, nx, ny, nz, icommy, icommz, npuy, npuz, -1)


	do k=1,nzp
		do j=1,nyp
			do i=1,nx
				ddx_nueff_s11(i,j,k) = imag*kx(i)*2.0_rp*nueff_s11(i,j,k)
				ddy_nueff_s12(i,j,k) = imag*kyp(j)*2.0_rp*nueff_s12(i,j,k)
				ddz_nueff_s13(i,j,k) = imag*kzp(k)*2.0_rp*nueff_s13(i,j,k)
				ddx_nueff_s21(i,j,k) = imag*kx(i)*2.0_rp*nueff_s21(i,j,k)
				ddy_nueff_s22(i,j,k) = imag*kyp(j)*2.0_rp*nueff_s22(i,j,k)
				ddz_nueff_s23(i,j,k) = imag*kzp(k)*2.0_rp*nueff_s23(i,j,k)
				ddx_nueff_s31(i,j,k) = imag*kx(i)*2.0_rp*nueff_s31(i,j,k)
				ddy_nueff_s32(i,j,k) = imag*kyp(j)*2.0_rp*nueff_s32(i,j,k)
				ddz_nueff_s33(i,j,k) = imag*kzp(k)*2.0_rp*nueff_s33(i,j,k)
			enddo
		enddo
	enddo
	
    
    call PZFFT3DV (ddx_nueff_s11, ddx_nueff_s11, nx, ny, nz, icommy, icommz, npuy, npuz, +1)
   	call PZFFT3DV (ddy_nueff_s12, ddy_nueff_s12, nx, ny, nz, icommy, icommz, npuy, npuz, +1)
   	call PZFFT3DV (ddz_nueff_s13, ddz_nueff_s13, nx, ny, nz, icommy, icommz, npuy, npuz, +1)
   	call PZFFT3DV (ddx_nueff_s21, ddx_nueff_s21, nx, ny, nz, icommy, icommz, npuy, npuz, +1)
   	call PZFFT3DV (ddy_nueff_s22, ddy_nueff_s22, nx, ny, nz, icommy, icommz, npuy, npuz, +1)
   	call PZFFT3DV (ddz_nueff_s23, ddz_nueff_s23, nx, ny, nz, icommy, icommz, npuy, npuz, +1)
   	call PZFFT3DV (ddx_nueff_s31, ddx_nueff_s31, nx, ny, nz, icommy, icommz, npuy, npuz, +1)
   	call PZFFT3DV (ddy_nueff_s32, ddy_nueff_s32, nx, ny, nz, icommy, icommz, npuy, npuz, +1)
   	call PZFFT3DV (ddz_nueff_s33, ddz_nueff_s33, nx, ny, nz, icommy, icommz, npuy, npuz, +1)




!	!compute products ==================================================
	
	!NS-x
	ns_x_1 = ue*dudxe + ve*dudye + we*dudze
	ns_x_2 = ddx_nueff_s11 + ddy_nueff_s12 + ddz_nueff_s13
			
	!NS-y
	ns_y_1 = ue*dvdxe + ve*dvdye + we*dvdze
	ns_y_2 = ddx_nueff_s21 + ddy_nueff_s22 + ddz_nueff_s23
			
	!NS-z
	ns_z_1 = ue*dwdxe + ve*dwdye + we*dwdze
	ns_z_2 = ddx_nueff_s31 + ddy_nueff_s32 + ddz_nueff_s33
			

		
!	!evaluate product in Fourier space (R -> F)
	aux = ns_x_1
	call PZFFT3DV (aux, ns_x_1, nx, ny, nz, icommy, icommz, npuy, npuz, -1)
	aux = ns_x_2
	call PZFFT3DV (aux, ns_x_2, nx, ny, nz, icommy, icommz, npuy, npuz, -1)
	aux = ns_y_1
	call PZFFT3DV (aux, ns_y_1, nx, ny, nz, icommy, icommz, npuy, npuz, -1)
	aux = ns_y_2
	call PZFFT3DV (aux, ns_y_2, nx, ny, nz, icommy, icommz, npuy, npuz, -1)
	aux = ns_z_1
	call PZFFT3DV (aux, ns_z_1, nx, ny, nz, icommy, icommz, npuy, npuz, -1)
	aux = ns_z_2
	call PZFFT3DV (aux, ns_z_2, nx, ny, nz, icommy, icommz, npuy, npuz, -1)


!	!evaluate rhsn in Fourier space (F)
	rhsx      = -ns_x_1 + ns_x_2 
	rhsy      = -ns_y_1 + ns_y_2
	rhsz      = -ns_z_1 + ns_z_2


! 	!evaluate projections (F)
 	call projection(nx,ny,nz,nyp,nzp,kx,kyp,kzp,k_2p,rhsx,rhsy,rhsz)


	!compute invariants ================================================
	if (vort_flag == 1) then
	
		!Components of vorticity vector:
		vort_1 = dwdye - dvdze
		vort_2 = dudze - dwdxe
		vort_3 = dvdxe - dudye
		
		!Components of tensors S and W: 
		s_11 = dudxe
		s_12 = 0.50_rp*(dvdxe + dudye)
		s_13 = 0.50_rp*(dwdxe + dudze)
		s_21 = s_12
		s_22 = dvdye
		s_23 = 0.50_rp*(dwdye + dvdze)
		s_31 = s_13 
		s_32 = s_23
		s_33 = dwdze
		
		w_11 = 0.0_rp
		w_12 = 0.50_rp*(dvdxe - dudye)
		w_13 = 0.50_rp*(dwdxe - dudze)
		w_21 = -w_12
		w_22 = 0.0_rp
		w_23 = 0.50_rp*(dwdye - dvdze)
		w_31 = -w_13
		w_32 = -w_23
		w_33 = 0.0_rp
		
		!Evaluate vorticity magnitude:
		vort_mag = sqrt(vort_1**2 + vort_2**2 + vort_3**2)
		
		
		!Evaluate 2nd invariant (Q) of the gradient of velocity tensor
		qs    = -0.50_rp*(s_11**2 + s_22**2 + s_33**2 + s_12*s_21 + s_13*s_31 + &
								 s_21*s_12 + s_23*s_32 + s_31*s_13 + s_32*s_23	)
		qw    = -0.50_rp*( w_12*w_21 + w_13*w_31 + w_21*w_12 + w_23*w_32 + &
								 w_31*w_13 + w_32*w_23 )			
		qcrit = qs + qw	
		
		!Evaluate 3rd invariant (R) of the gradient of velocity tensor
		rs 		    = (-1.0_rp/3.0_rp)*(s_11**3 + s_11*s_12*s_21 + s_11*s_13*s_31 + &
							    s_12*s_21*s_11 + s_12*s_22*s_21 + s_12*s_23*s_31 + &
							    s_13*s_31*s_11 + s_13*s_32*s_21 + s_13*s_33*s_31 + &
							    s_21*s_11*s_12 + s_21*s_12*s_22 + s_21*s_13*s_32 + &
							    s_22*s_21*s_12 + s_22**3 + s_22*s_23*s_32 + &
							    s_23*s_31*s_12 + s_23*s_32*s_22 + s_23*s_33*s_32 + &
							    s_31*s_11*s_13 + s_31*s_12*s_23 + s_31*s_13*s_33 + &
							    s_32*s_21*s_13 + s_32*s_22*s_23 + s_32*s_23*s_33 + &
							    s_33*s_31*s_13 + s_33*s_32*s_23 + s_33**3)
		
		sigma_stret  = -(w_11*w_11*s_11 + w_11*w_12*s_21 + w_11*w_13*s_31 + &
								w_12*w_21*s_11 + w_12*w_22*s_21 + w_12*w_23*s_31 + &
								w_13*w_31*s_11 + w_13*w_32*s_21 + w_13*w_33*s_31 + &
								w_21*w_11*s_12 + w_21*w_12*s_22 + w_21*w_13*s_32 + &
								w_22*w_21*s_12 + w_22*w_22*s_22 + w_22*w_23*s_32 + &
								w_23*w_31*s_12 + w_23*w_32*s_22 + w_23*w_33*s_32 + &
								w_31*w_11*s_13 + w_31*w_12*s_23 + w_31*w_13*s_33 + &
								w_32*w_21*s_13 + w_32*w_22*s_23 + w_32*w_23*s_33 + &
								w_33*w_31*s_13 + w_33*w_32*s_23 + w_33*w_33*s_33 )
		
		rcrit = rs + sigma_stret	
				
	endif

 end subroutine evaluate_rhs_incompressible
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++




!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   
subroutine temporal_integration_incompressible(nx,ny,nz,nyp,nzp,x,yp,zp,tc,kx,kyp,kzp,k_2p,icommy,icommz,npuy,npuz,imag,mu,dt,alphark,betark,crk,&
temp_integration_set,u,v,w,pressure,vort_1,vort_2,vort_3,vort_mag,qcrit,ab_u,ab_v,ab_w,dx,dy,dz)

implicit none
integer,parameter:: rp=selected_real_kind(8,200),sp=selected_int_kind(r=8)


 !intent in
 integer(sp)::nx,ny,nz,nyp,nzp,tc,icommy,icommz,npuy,npuz
 real(rp),dimension(nx),intent(in)::x,yp,zp
 complex(rp),dimension(nx),intent(in)::kx
 complex(rp),dimension(nyp),intent(in)::kyp
 complex(rp),dimension(nzp),intent(in)::kzp
 complex(rp),dimension(nx,nyp,nzp),intent(in)::k_2p
 complex(rp)::imag
 real(rp)::mu,dt,dx,dy,dz
 real(rp),dimension(7):: alphark,betark,crk


 !intent in/out
 integer(sp)::i,j,k,vort_flag,flag_AB,temp_integration_set
 complex(rp),dimension(nx,nyp,nzp)::urk,vrk,wrk,w_u,w_v,w_w,&
 rhsx,rhsy,rhsz,big_u,big_v,big_w,k1_u,k1_v,k1_w,&
 k2_u,k2_v,k2_w,k3_u,k3_v,k3_w,k4_u,k4_v,k4_w, ns_x_1,ns_y_1,ns_z_1
 real(rp)::trk
 complex(rp),dimension(nx,nyp,nzp)::u,v,w,pressure,vort_1,vort_2,vort_3,vort_mag,qcrit
 complex(rp),dimension(nx,nyp,nzp,4)::ab_u,ab_v,ab_w


 urk      = 0.0_rp
 vrk      = 0.0_rp
 wrk      = 0.0_rp
 rhsx     = 0.0_rp
 rhsy     = 0.0_rp
 rhsz     = 0.0_rp


!===============================================================
!evaluate fields in Fourier space (F) using AB4				   =
!===============================================================

if (temp_integration_set==1) then

flag_AB=1
	
	if (tc<=3) then
	temp_integration_set=2

	
	else if (tc>3) then
					            
	    !Evaluate fields
		u      = u + dt*(23.0_rp*ab_u(:,:,:,3) - 16.0_rp*ab_u(:,:,:,2) + 5.0_rp*ab_u(:,:,:,1) )/12.0_rp
		v      = v + dt*(23.0_rp*ab_v(:,:,:,3) - 16.0_rp*ab_v(:,:,:,2) + 5.0_rp*ab_v(:,:,:,1) )/12.0_rp
		w      = w + dt*(23.0_rp*ab_w(:,:,:,3) - 16.0_rp*ab_w(:,:,:,2) + 5.0_rp*ab_w(:,:,:,1) )/12.0_rp		 
		 
		 		
		!Evaluate pressure field
		do k=1,nzp
			do j=1,nyp
				do i=1,nx
					pressure(i,j,k)=imag*1.0_rp/k_2p(i,j,k)*( (ns_x_1(i,j,k))*&
					kx(i) + (ns_y_1(i,j,k))*kyp(j) + (ns_z_1(i,j,k))*kzp(k) ) 
				enddo
			enddo
		enddo
		 
		!Update Adams-Bashforth coefficients
		
			!Evaluate Adams-Bashforth coefficients for present time step
			call evaluate_rhs_incompressible(nx,ny,nz,nyp,nzp,kx,kyp,kzp,k_2p,icommy,icommz,npuy,npuz,mu,u,v,w,imag,&
				rhsx,rhsy,rhsz,ns_x_1,ns_y_1,ns_z_1,vort_1,vort_2,vort_3,vort_mag,qcrit,1,tc)
			
			!Update Adams-Bashforth coefficients
			!Coefficients 1,2 (1 <- 2, 2 <- 3)
			do i=1,2
				ab_u(:,:,:,i)      = ab_u(:,:,:,i+1) 
				ab_v(:,:,:,i)      = ab_v(:,:,:,i+1) 
				ab_w(:,:,:,i)      = ab_w(:,:,:,i+1) 
			enddo
			
			!Coefficient 3: (3 <- rhs from fields just calculated)			              
			ab_u(:,:,:,3)      = rhsx(:,:,:)
			ab_v(:,:,:,3)      = rhsy(:,:,:)
			ab_w(:,:,:,3)      = rhsz(:,:,:)
			
	endif

endif

!===============================================================
!evaluate fields in Fourier space (F) using RK46 Scheme 	   =
!===============================================================

if (temp_integration_set==2) then
	urk    = u
	vrk    = v
	wrk    = w
	
	w_u      = 0.0_rp
	w_v      = 0.0_rp
	w_w      = 0.0_rp

	do j=1,7
		
		if (j==7) then
			vort_flag = 1 !return vorticity and second invariant Q
		else
			vort_flag = 0 !do not return vorticity and Q
		endif
		
		trk=(tc-1.0_rp+crk(j))*dt
		call evaluate_rhs_incompressible(nx,ny,nz,nyp,nzp,kx,kyp,kzp,k_2p,icommy,icommz,npuy,npuz,mu,urk,vrk,wrk,imag,&
		rhsx,rhsy,rhsz,ns_x_1,ns_y_1,ns_z_1,vort_1,vort_2,vort_3,vort_mag,qcrit,vort_flag,tc)
					              
						
		big_u      = rhsx      
		big_v      = rhsy      
		big_w      = rhsz      

		w_u        = w_u*alphark(j)  + dt*big_u
		w_v        = w_v*alphark(j)  + dt*big_v
		w_w        = w_w*alphark(j)  + dt*big_w

		urk        = urk      + betark(j)*w_u
		vrk        = vrk      + betark(j)*w_v
		wrk        = wrk      + betark(j)*w_w

		!trk=(tc-1.0d0+crk(j))*dt
	enddo
	
	u   = urk
	v   = vrk
	w   = wrk
	

	!Evaluate pressure field
	do k=1,nzp
		do j=1,nyp
			do i=1,nx
				pressure(i,j,k)=imag*1.0_rp/k_2p(i,j,k)*( (ns_x_1(i,j,k))*&
				kx(i) + (ns_y_1(i,j,k))*kyp(j) + (ns_z_1(i,j,k))*kzp(k) ) 
			enddo
		enddo
	enddo
	
	
	!Controls for AB4 case (4 first solutions are calculated using RK46)
	if (flag_AB==1 .and. tc<=3) then
	
		call evaluate_rhs_incompressible(nx,ny,nz,nyp,nzp,kx,kyp,kzp,k_2p,icommy,icommz,npuy,npuz,mu,u,v,w,imag,&
		rhsx,rhsy,rhsz,ns_x_1,ns_y_1,ns_z_1,vort_1,vort_2,vort_3,vort_mag,qcrit,vort_flag,tc)
					              
		ab_u(:,:,:,tc)      = rhsx(:,:,:)
		ab_v(:,:,:,tc)      = rhsy(:,:,:)
		ab_w(:,:,:,tc)      = rhsz(:,:,:)
		
		if (tc==3) then
			temp_integration_set=1
		endif
		
		
	endif
	
	
	!===============================================================
	!evaluate fields in Fourier space (F) using Euler Scheme       =
	!===============================================================
	else if (temp_integration_set==3) then
	
		call evaluate_rhs_incompressible(nx,ny,nz,nyp,nzp,kx,kyp,kzp,k_2p,icommy,icommz,npuy,npuz,mu,u,v,w,imag,&
		rhsx,rhsy,rhsz,ns_x_1,ns_y_1,ns_z_1,vort_1,vort_2,vort_3,vort_mag,qcrit,1,tc)

		k1_u       = rhsx      
		k1_v       = rhsy      
		k1_w       = rhsz      		

		u       = u   + dt*k1_u 
		v       = v   + dt*k1_v
		w       = w   + dt*k1_w
		

	!Evaluate pressure field
	do k=1,nzp
		do j=1,nyp
			do i=1,nx
				pressure(i,j,k)=imag*1.0_rp/k_2p(i,j,k)*( (ns_x_1(i,j,k))*&
				kx(i) + (ns_y_1(i,j,k))*kyp(j) + (ns_z_1(i,j,k))*kzp(k) ) 
			enddo
		enddo
	enddo



	!===============================================================
	!evaluate fields in Fourier space (F) using RK2				   =
	!===============================================================
	else if (temp_integration_set==4) then
	
		!K1
		!call forcing terms
		call evaluate_rhs_incompressible(nx,ny,nz,nyp,nzp,kx,kyp,kzp,k_2p,icommy,icommz,npuy,npuz,mu,u,v,w,imag,&
		rhsx,rhsy,rhsz,ns_x_1,ns_y_1,ns_z_1,vort_1,vort_2,vort_3,vort_mag,qcrit,0,tc)

		k1_u      = rhsx
		k1_v      = rhsy
		k1_w      = rhsz
	
	
		urk       = u       + k1_u*dt
		vrk       = v       + k1_v*dt
		wrk       = w       + k1_w*dt
			
	
		!K2
		!call forcing terms
		call evaluate_rhs_incompressible(nx,ny,nz,nyp,nzp,kx,kyp,kzp,k_2p,icommy,icommz,npuy,npuz,mu,urk,vrk,wrk,imag,&
		rhsx,rhsy,rhsz,ns_x_1,ns_y_1,ns_z_1,vort_1,vort_2,vort_3,vort_mag,qcrit,1,tc)


		k2_u      = rhsx      
		k2_v      = rhsy      
		k2_w      = rhsz      
	
		!RK approx
		u      = u      + dt*( k1_u      + k2_u )/2.0_rp
		v      = v      + dt*( k1_v      + k2_v )/2.0_rp
		w      = w      + dt*( k1_w      + k2_w )/2.0_rp

	
	!Evaluate pressure field
	do k=1,nzp
		do j=1,nyp
			do i=1,nx
				pressure(i,j,k)=imag*1.0_rp/k_2p(i,j,k)*( (ns_x_1(i,j,k))*&
				kx(i) + (ns_y_1(i,j,k))*kyp(j) + (ns_z_1(i,j,k))*kzp(k) ) 
			enddo
		enddo
	enddo

	!===============================================================
	!evaluate fields in Fourier space (F) using RK4				   =
	!===============================================================
	else if (temp_integration_set==5) then

		!K1
		!call forcing terms
		call evaluate_rhs_incompressible(nx,ny,nz,nyp,nzp,kx,kyp,kzp,k_2p,icommy,icommz,npuy,npuz,mu,urk,vrk,wrk,imag,&
		rhsx,rhsy,rhsz,ns_x_1,ns_y_1,ns_z_1,vort_1,vort_2,vort_3,vort_mag,qcrit,0,tc)
		
		k1_u       = rhsx      
		k1_v       = rhsy      
		k1_w       = rhsz      
		
	
		urk      = u      + k1_u*dt/2.0_rp
		vrk      = v      + k1_v*dt/2.0_rp
		wrk      = w      + k1_w*dt/2.0_rp
				
	
		!K2
		!call forcing terms
		call evaluate_rhs_incompressible(nx,ny,nz,nyp,nzp,kx,kyp,kzp,k_2p,icommy,icommz,npuy,npuz,mu,urk,vrk,wrk,imag,&
		rhsx,rhsy,rhsz,ns_x_1,ns_y_1,ns_z_1,vort_1,vort_2,vort_3,vort_mag,qcrit,0,tc)

		k2_u      = rhsx      
		k2_v      = rhsy      
		k2_w      = rhsz      

		urk      = u + k2_u*dt/2.0_rp
		vrk      = v + k2_v*dt/2.0_rp
		wrk      = w + k2_w*dt/2.0_rp
			
	
		!K3
		!call forcing terms
		call evaluate_rhs_incompressible(nx,ny,nz,nyp,nzp,kx,kyp,kzp,k_2p,icommy,icommz,npuy,npuz,mu,urk,vrk,wrk,imag,&
		rhsx,rhsy,rhsz,ns_x_1,ns_y_1,ns_z_1,vort_1,vort_2,vort_3,vort_mag,qcrit,0,tc)
	
		k3_u      = rhsx      
		k3_v      = rhsy      
		k3_w      = rhsz      
	
		urk      = u + k3_u*dt
		vrk      = v + k3_v*dt
		wrk      = w + k3_w*dt		
		
		
		!K4
		!call forcing terms
		call evaluate_rhs_incompressible(nx,ny,nz,nyp,nzp,kx,kyp,kzp,k_2p,icommy,icommz,npuy,npuz,mu,urk,vrk,wrk,imag,&
		rhsx,rhsy,rhsz,ns_x_1,ns_y_1,ns_z_1,vort_1,vort_2,vort_3,vort_mag,qcrit,1,tc)
	
		k4_u      = rhsx      
		k4_v      = rhsy      
		k4_w      = rhsz      
	
	
		!RK approximation
		u      = u + dt*( k1_u + 2.0_rp*k2_u + 2.0_rp*k3_u + k4_u )/6.0_rp
		v      = v + dt*( k1_v + 2.0_rp*k2_v + 2.0_rp*k3_v + k4_v )/6.0_rp
		w      = w + dt*( k1_w + 2.0_rp*k2_w + 2.0_rp*k3_w + k4_w )/6.0_rp	
		
			
	!Evaluate pressure field
	do k=1,nzp
		do j=1,nyp
			do i=1,nx
				pressure(i,j,k)=imag*1.0_rp/k_2p(i,j,k)*( (ns_x_1(i,j,k))*&
				kx(i) + (ns_y_1(i,j,k))*kyp(j) + (ns_z_1(i,j,k))*kzp(k) ) 
			enddo
		enddo
	enddo

	
endif


end subroutine temporal_integration_incompressible
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 subroutine projection(nx,ny,nz,nyp,nzp,kx,kyp,kzp,k_2p,rhsx,rhsy,rhsz)

 implicit none
 integer,parameter:: rp=selected_real_kind(8,200),sp=selected_int_kind(r=8)

 !intent in
 integer(sp),intent(in)::nx,ny,nz,nyp,nzp	
 complex(rp),dimension(nx)::kx
 complex(rp),dimension(nyp)::kyp
 complex(rp),dimension(nzp)::kzp
 complex(rp),dimension(nx,nyp,nzp)::k_2p

 !intent in/out
 integer(sp)::i,j,k
 complex(rp),dimension(nx,nyp,nzp)::rhsx,rhsy,rhsz,proj_tensor
 
 
 
do k=1,nzp
	do j=1,nyp
		do i=1,nx
			proj_tensor(i,j,k)=(rhsx(i,j,k)*kx(i)+rhsy(i,j,k)*kyp(j)+rhsz(i,j,k)*kzp(k))/k_2p(i,j,k)
			rhsx(i,j,k)=rhsx(i,j,k)-proj_tensor(i,j,k)*kx(i)
			rhsy(i,j,k)=rhsy(i,j,k)-proj_tensor(i,j,k)*kyp(j)
			rhsz(i,j,k)=rhsz(i,j,k)-proj_tensor(i,j,k)*kzp(k)
		 end do
	 end do
 enddo
 
 end subroutine projection
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++




!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 subroutine filter_field(nx,ny,nz,nyp,nzp,filter,u)

 implicit none
 integer,parameter:: rp=selected_real_kind(8,200),sp=selected_int_kind(r=8)
	
 !intent in
 integer(sp),intent(in)::nx,ny,nz,nyp,nzp	
 complex(rp),dimension(nx,nyp,nzp)::filter
	
 !intent in/out
 complex(rp),dimension(nx,nyp,nzp)::u


 u       = filter*u

 end subroutine filter_field
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++




!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 subroutine screen_info(max_u,max_v,max_w,tc,t,dt,solver_set)
 
 implicit none
 integer,parameter:: rp=selected_real_kind(8,200),sp=selected_int_kind(r=8)

 !intent in
 integer(sp)::tc,solver_set,id
 real(rp)::max_u,max_v,max_w,t,dt
 
 if (id == 0) then
 write(*,*)'*********************************************'
 write(*,*)'MAX U = ',max_u
 write(*,*)'MAX V = ',max_v
 write(*,*)'MAX W = ',max_w
 write(*,*)'TIME STEP tc = ',tc
 write(*,*)'PHYSICAL TIME = ',t+dt
 write(*,*)'*********************************************'
 endif
 
 end subroutine screen_info
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++




!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   
!subroutine: random_generator
 subroutine random_generator(nx,ny,nz,nyp,nzp,randomsetu,randomsetv,randomsetw)

implicit none
integer,parameter:: rp=selected_real_kind(8,200),sp=selected_int_kind(r=8)

integer(sp) :: nx,ny,nz,nyp,nzp
real(rp),dimension(nx,nyp,nzp) :: randomsetu
real(rp),dimension(nx,nyp,nzp) :: randomsetv
real(rp),dimension(nx,nyp,nzp) :: randomsetw
real(rp),dimension(nx,nyp,nzp) :: randomsetrho

call random_seed()

call random_number(randomsetu)
randomsetu=(2.0_rp*randomsetu-1.0_rp)

call random_number(randomsetv)
randomsetv=(2.0_rp*randomsetv-1.0_rp)

call random_number(randomsetw)
randomsetw=(2.0_rp*randomsetw-1.0_rp)

call random_number(randomsetrho)
randomsetrho=(2.0_rp*randomsetrho-1.0_rp)


!write(*,*) (arr(n),n=1,10)
end subroutine random_generator 
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   





!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   
subroutine grid_to_print(field,prfield,nx,nyp,nzp,node)
!> Parallel output: generates a new scalar field for printing, based on the scalar field.
! "field" the original field (used by the solver)
! "prfield" the field to be generated for impression prfield is a copy of the original scalar field,
! with additional data to fill the gaps generated by the parallel output

	complex(rp), dimension ( nx, nyp, nzp ) :: field
	real(rp), dimension ( 0:nx, 0:nyp, 0:nzp ) :: prfield
	
	real(rp), dimension ( 0:nx ) :: x_exch
	real(rp), dimension ( (nx+1)*nzp ) :: y_exch
	real(rp), dimension ( (nx+1)*nyp ) :: z_exch
    type (mapnode) :: node
	integer:: status(MPI_STATUS_SIZE)

	integer(sp) :: i, j, k, gp_idx
	
	call MPI_BARRIER (MPI_COMM_WORLD, ierr)
	
	! field interpolation
	prfield(1:nx, 1:nyp, 1:nzp) = real ( field )
	
	! x-direction: periodicity
	prfield (0, 1:nyp, 1:nzp) = prfield (nx, 1:nyp, 1:nzp)
	
	! y-direction: takes data from process (idy, idz-1). Periodicity also taken into account
	if (npuy > 1) then
		if (idy == 0) then 
			
			gp_idx = 1
			do k = 1, nzp
				do i = 0, nx
					y_exch (gp_idx) = prfield (i, nyp, k)
					gp_idx = gp_idx + 1
				enddo
			enddo
			
			call MPI_SEND (y_exch, (nx+1)*nzp, MPI_DOUBLE_PRECISION, node%E, id, MPI_COMM_WORLD, ierr)
			call MPI_RECV (y_exch, (nx+1)*nzp, MPI_DOUBLE_PRECISION, node%W, node%W, MPI_COMM_WORLD, status, ierr)
			
			gp_idx = 1
			do k = 1, nzp
				do i = 0, nx
					prfield (i, 0, k) = y_exch (gp_idx)
					gp_idx = gp_idx+1
				enddo
			enddo
			
		else 
		
			call MPI_RECV (y_exch, (nx+1)*nzp, MPI_DOUBLE_PRECISION, node%W, node%W, MPI_COMM_WORLD, status, ierr)
			
			gp_idx = 1
			do k = 1, nzp
				do i = 0, nx
					prfield (i, 0, k) = y_exch (gp_idx)
					gp_idx = gp_idx + 1
				enddo
			enddo
			
			gp_idx = 1
			do k = 1, nzp
				do i = 0, nx
					y_exch (gp_idx) = prfield (i, nyp, k)
					gp_idx = gp_idx + 1
				enddo
			enddo
			
			call MPI_SEND (y_exch, (nx+1)*nzp, MPI_DOUBLE_PRECISION, node%E, id, MPI_COMM_WORLD, ierr)
		
		endif
	else
		
		prfield (0:nx, 0, 1:nzp) = prfield (0:nx, nyp, 1:nzp)

	endif
	
	
	! z-direction: takes data from process (idy, idz-1). Periodicity also taken into account
	if (npuz > 1) then
		if (idz == 0) then 
			
			gp_idx = 1
			do j = 1, nyp
				do i = 0, nx
					z_exch (gp_idx) = prfield (i, j, nzp)
					gp_idx = gp_idx + 1
				enddo
			enddo
			
			call MPI_SEND (z_exch, (nx+1)*nyp, MPI_DOUBLE_PRECISION, node%N, id, MPI_COMM_WORLD, ierr)
			call MPI_RECV (z_exch, (nx+1)*nyp, MPI_DOUBLE_PRECISION, node%S, node%S, MPI_COMM_WORLD, status, ierr)
			
			gp_idx = 1
			do j = 1, nyp
				do i = 0, nx
					prfield (i, j, 0) = z_exch (gp_idx)
					gp_idx = gp_idx + 1
				enddo
			enddo
			
		else
		
			call MPI_RECV (z_exch, (nx+1)*nyp, MPI_DOUBLE_PRECISION, node%S, node%S, MPI_COMM_WORLD, status, ierr)
			
			gp_idx = 1
			do j = 1, nyp
				do i = 0, nx
					prfield (i, j, 0) = z_exch (gp_idx)
					gp_idx = gp_idx + 1
				enddo
			enddo
			
			gp_idx = 1
			do j = 1, nyp
				do i = 0, nx
					z_exch (gp_idx) = prfield(i, j, nzp)
					gp_idx = gp_idx + 1
				enddo
			enddo
			
			call MPI_SEND (z_exch, (nx+1)*nyp, MPI_DOUBLE_PRECISION, node%N, id, MPI_COMM_WORLD, ierr)
		
		endif
	else
		prfield (0:nx, 1:nyp, 0) = prfield(0:nx, 1:nyp, nzp)
	endif
	
	! corners: takes data from process (idy-1, idz-1). Periodicity also taken into account
	if (npu > 1) then
	
		if (idz == 0) then
		
			x_exch(0:nx) = prfield (0:nx, nyp, nzp)
				
			call MPI_SEND (x_exch, nx+1, MPI_DOUBLE_PRECISION, node%NE, id, MPI_COMM_WORLD, ierr)
			call MPI_RECV (x_exch, nx+1, MPI_DOUBLE_PRECISION, node%SW, node%SW, MPI_COMM_WORLD, status, ierr)
			
			prfield (0:nx, 0, 0) = x_exch (0:nx)
			
		else
			call MPI_RECV (x_exch, nx+1, MPI_DOUBLE_PRECISION, node%SW, node%SW, MPI_COMM_WORLD, status, ierr)
			
			prfield (0:nx, 0, 0) = x_exch (0:nx)
			
			
			x_exch (0:nx) = prfield (0:nx, nyp, nzp)
			
			call MPI_SEND (x_exch, nx+1, MPI_DOUBLE_PRECISION, node%NE, id, MPI_COMM_WORLD, ierr)
		endif
	else

		prfield(0:nx, 0, 0) = prfield(0:nx, nyp, nzp)
		
	endif
	
	return
	
end subroutine grid_to_print
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   





!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   
subroutine evaluate_spectrum(nx,nyp,nzp,n3,k0,dk,kmax,k_2p,freq_spectrum,tc,ax1,ay1,az1)
	
	implicit none
	integer,parameter:: rp=selected_real_kind(8,200),sp=selected_int_kind(r=8)
	integer(sp)::ik,i3,i4,k_size,i,j,k,freq_spectrum,tc,ierr,nx,nyp,nzp,n3
	integer::status(MPI_STATUS_SIZE)
	real(rp)::sum_k,sum2_k,mean_k,kik,kmax,k0,dk
	real(rp),allocatable,dimension(:)::ek, kt
	complex(rp),intent(in),dimension(nx,nyp,nzp)::ax1,ay1,az1,k_2p
	complex(rp),allocatable,dimension(:,:,:)::axn,ayn,azn
	312	format(1(es12.5),',',1(1x,es17.10))
	200	format(a10,1x,a40)
	character(80)::str1,str2,str3,str4,filename,filename2

	str1 = '../RESULTS/STATISTICS/'

	allocate(axn(nx,nyp,nzp),ayn(nx,nyp,nzp),azn(nx,nyp,nzp))

	axn=ax1/n3
	ayn=ay1/n3
	azn=az1/n3
	
	k_size = int((kmax-k0)/dk)
	allocate(ek(k_size), kt(k_size))
	
	ek = 0.0_rp

	do ik = 1,k_size
		i3=0
		sum_k  = 0.0_rp
		mean_k = 0.0_rp
		kt(ik)=ik*dk
		do k=1,nzp
			do j=1,nyp
				do i=1,nx
					kik=sqrt(k_2p(i,j,k))
					if ((kik>=(kt(ik)-dk/2.0_rp)).and.(kik<=(kt(ik)+dk/2.0_rp))) then
						i3=i3+1
						sum_k=sum_k+(axn(i,j,k)*conjg(axn(i,j,k))+ayn(i,j,k)*conjg(ayn(i,j,k))+azn(i,j,k)*conjg(azn(i,j,k)))
					endif
				enddo
			enddo
		enddo
		call MPI_BARRIER (MPI_COMM_WORLD, ierr)
		call MPI_ALLREDUCE(sum_k,sum2_k,1, MPI_DOUBLE_PRECISION,MPI_SUM, MPI_COMM_WORLD, ierr)
		call MPI_ALLREDUCE(i3,i4,1, MPI_INTEGER,MPI_SUM, MPI_COMM_WORLD, ierr)
		if (i4 .ne. 0) then
			mean_k=sum2_k/i4
		endif
		ek(ik)=2.0_rp*pi*(kt(ik)**2)*mean_k
	enddo
	deallocate(axn,ayn,azn)
	
	if (id==0) then
		str2 = 'TKE_spectrum'
		str4 = 'Kolmogorov_slope'
		write(str3,*) tc
		str3=adjustl(str3)
		filename = trim(str1)//trim(str2)//trim(str3)//'.csv'
		open(30, file = trim(filename),status='unknown')
		write(30,*)'variables= "k", "ek",\'
		if (tc==freq_spectrum) then
			do ik=1,k_size
				filename2 = trim(str1)//trim(str4)//trim(str3)//'.csv'
				open(40, file = trim(filename2),status='unknown')
				write(40,312) kt(ik),(kt(ik))**(-5.0d0/3.0d0) !Kolmogorov slope -5/3
			enddo
			close(40)
		endif
		do ik=1,k_size
			write(30,312) kt(ik),ek(ik)
		enddo
		close(30)
	endif
	
	deallocate(ek,kt)

end subroutine evaluate_spectrum
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   






!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   
!SUBROUTINE spatial_average

!	!REAL(rp), DIMENSION(nx, nyp, nzp) :: axn, axf
!  	integer(sp) :: j
! 	real(rp), allocatable, dimension (:) :: us, vs, ws, usp, vsp, wsp, ustot, vstot, wstot
 	
!	allocate ( usp (nyp), vsp(nyp), wsp(nyp) )
	
!	if (idz == 0) then
!		allocate ( us (nyp), vs(nyp), ws(nyp) )
!	endif

 	
! 	if (id == 0) then
!		allocate ( ustot (ny), vstot (ny), wstot (ny) )
!	endif
 	
! 	usp = 0.0_rp
! 	vsp = 0.0_rp
! 	wsp = 0.0_rp
 	
! 	do j = 1, nyp
!		usp(j) = sum( real(axx(1:nx,j,1:nzp)) )
!		vsp(j) = sum( real(axy(1:nx,j,1:nzp)) )
!		wsp(j) = sum( real(axz(1:nx,j,1:nzp)) )
! 	enddo
 	
! 	!call MPI_BARRIER
! 	call MPI_REDUCE ( usp , us, nyp, MPI_DOUBLE_PRECISION, MPI_SUM, 0, icommz, ierr )
!	call MPI_REDUCE ( vsp , vs, nyp, MPI_DOUBLE_PRECISION, MPI_SUM, 0, icommz, ierr )
!	call MPI_REDUCE ( wsp , ws, nyp, MPI_DOUBLE_PRECISION, MPI_SUM, 0, icommz, ierr )
	
	
!	deallocate ( usp, vsp, wsp )
	
!	if (idz == 0) then
!		us = us/(nx*nz)
!		vs = vs/(nx*nz)
!		ws = ws/(nx*nz)
!		call MPI_GATHER ( us, nyp, MPI_DOUBLE_PRECISION, ustot, nyp, MPI_DOUBLE_PRECISION, 0, icommy, ierr ) 
!		call MPI_GATHER ( vs, nyp, MPI_DOUBLE_PRECISION, vstot, nyp, MPI_DOUBLE_PRECISION, 0, icommy, ierr ) 
!		call MPI_GATHER ( ws, nyp, MPI_DOUBLE_PRECISION, wstot, nyp, MPI_DOUBLE_PRECISION, 0, icommy, ierr ) 
!		deallocate ( us, vs, ws)
!	endif
	
!	call MPI_BARRIER ( MPI_COMM_WORLD, ierr)
 	
 	
 	
! 	! Printing
! 	if (id == 0) then
!		str2 = 'mean_profiles'
!		write(str3,*) nstep2
!		str3=ADJUSTL(str3)
!		filename = TRIM(str1)//TRIM(str2)//TRIM(str3)//'.csv'
!		open(30, file = TRIM(filename),status='unknown')
!		!WRITE(30,*)'variables= "y", "ek",\'
!		write(30,*)'y,us,vs,ws'
!		!WRITE(30,*)'zone i=',INT((kmax-k0)/dk)
		
!		do j = 1, ny
!			write(30,*) y(j), ustot(j), vstot(j), wstot(j) 
!		enddo
!		close(30)
!		deallocate ( ustot, vstot, wstot )
! 	endif
 	
	
!end subroutine spatial_average
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   


end module incompressible_solver

