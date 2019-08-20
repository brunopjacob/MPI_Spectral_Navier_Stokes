module declare_variables

use mpi

implicit none

integer,parameter:: rp=selected_real_kind(8,200),sp=selected_int_kind(r=8)  !Precision settings

integer(sp)::nx,ny,nz,nt,freqsave,freqsave_screen,freq_restart,freq_spectrum,vtk_counter,max_it_mdf

integer(sp)::i,j,k,l,tc,tc_restart,temp_integration_set,l2_case,ic_set,vort_flag,domain_set,temp_input_set,filter_set,solver_set,&
buffer_set,forcing_set,weight_function_set,n3,n3p,turbulence_modeling_set

!MPI-related variables
integer(sp)::npu,id,ierr,nyp,nzp,npuy,npuz,idy,idz,namelength,icommy,icommz,lnpu(3)
character (MPI_MAX_PROCESSOR_NAME)::procname
integer:: status(MPI_STATUS_SIZE)
real(rp), allocatable, dimension(:)::yp,zp,xsave,ypsave,zpsave
real(rp), allocatable, dimension(:,:,:)::usave,vsave,wsave,psave,&
vort_1_save,vort_2_save,vort_3_save,qcrit_save
complex(rp), allocatable, dimension(:)::kyp,kzp
complex(rp), allocatable, dimension(:,:,:)::k2_p
type :: mapnode
	integer (sp) :: n, ne, e, se, s, sw, w, nw
	real (rp) :: xmin, xmax, ymin, ymax, zmin, zmax
end type
type (mapnode) :: node



real(rp)::var_read,pi,Lx,Ly,Lz,dx,dy,dz,ds,L2,dt,dt1,dt2,t,mu,trk,t0,small,Re,&
max_u,max_v,max_w,max_u0,max_v0,max_w0,max_rho0,CFL,trestart,&
t_elapsed_start,t_elapsed_end,amplitude_perturb,k0,kmax,dk

real(rp),dimension(7):: alphark,betark,crk !Runge-Kutta HALE46 components

complex(rp), allocatable, dimension(:)::kx,ky,kz  !Wavenumbers

real(rp), allocatable, dimension(:)::x,y,z     !Mesh components

real(rp),allocatable,dimension(:,:,:)::fc !Filter for Smagorinsky

complex(rp), allocatable, dimension(:,:,:)::k_2p,u,ua,v,va,w,wa,pressure,unew,vnew,wnew,vort_1,vort_2,vort_3,vort_mag,qcrit,filter,aux,nueff

complex(rp),allocatable,dimension(:,:,:,:)::ab_u,ab_v,ab_w !Adams-Bashforth variables

complex(rp)::imag !Imaginary unit

character(25)::str1,str2
character(80)::line,line2,line3,line4,line5
character(80)::str3
character(10)::rank

logical:: restart_flag !If true, restart is activated


end module declare_variables
