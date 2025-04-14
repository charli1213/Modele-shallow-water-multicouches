      module data_initial
      integer exec_loc
      character(88) fftw_loc
      integer nx, ny, nz, nnx, nny
      integer nx_cou, ny_cou, nnx_cou, nny_cou
      integer i_diags
      double precision pi, twopi, Lx, Ly, dx, dy, H1, H2, H3, H4, H5, H6, Htot
      real f0, beta, r_drag, Ah2, Ah4, r_invLap, rf, g, alpha, fraction
      real tau0, tau1, wind_t0, variance, rho_atm
      real fileperday, daysperrestart
      integer nsteps,start_spec, cut_days
      real ndays,totaltime,dt
      real restart_from
      integer subsmprto,itape,ispechst,iout,itlocal,itsrow,ntsrow,nspecfile,endx,endy
      integer szsubx,szsuby
      integer ftsubsmprto,forcingtype, iou_method
      logical restart, use_ramp, ifsteady, gaussian_bump_eta
      logical calc1Dspec,save_movie,save2dfft
      real c_theta, c_mu, c_sigma, c_tauvar

      !
      ! I/O instruction for diognostics, to override, change parameters.f90
      logical  IO_field,   IO_forcing,  IO_QGAG
      logical  IO_psivort, IO_coupling, IO_RHS_uv
      logical  IO_psimodes
      logical  thickness_correction
      
      !
      ! >>> Defining WAVEWATCH III coupling variables >>>
      LOGICAL cou, ustar, waves, stokes, HS_fixed
      REAL step
      INTEGER :: ierror, numprocs, procid, err_class, err_len, iproc, numprocs_sec, procid_sec
      CHARACTER(80) :: err_str
      INTEGER :: MPI_SECOND
      INTEGER :: nxcou, nycou, mpiratio
      ! <<< Defining WAVEWATCH III coupling variables (End) <<<
      !
      
      include 'parameters.f90'
      parameter( ntsrow=itape/ispechst  )! how many lines for a time series table (e.g. spectrum)   
      !random number
      integer iseed,values(8)
      !forcing subroutine
      real delta_omega,f_thrhld
      integer n_omega,Fws
      real,allocatable :: omega(:),A_n(:),phi(:)
      !for complex number i
      double complex eye


      end module data_initial

      program main
      use data_initial
      use fishpack
      !!! MPI-coupling
      USE MPI
      !!! MPI-coupling
      implicit none
      !


      
    !!! ---------- Physical quantities definition ---------- !!!
      ! Mailles/Sides :
      REAL :: u(0:nnx,0:ny,nz,3), v(0:nx,0:nny,nz,3)
      REAL :: u_ag(0:nnx,0:ny,nz), v_ag(0:nx,0:nny,nz)
      REAL :: u_qg(0:nnx,0:ny,nz), v_qg(0:nx,0:nny,nz)
      REAL :: u_ag_p(0:nnx,0:ny,2), v_ag_p(0:nx,0:nny,2)
      REAL :: uBT(0:nnx,0:ny), vBT(0:nx,0:nny) 
      REAL :: uBC(0:nnx,0:ny,nz), vBC(0:nx,0:nny,nz)
      REAL :: grad2u(0:nnx,0:ny), grad4u(0:nnx,0:ny)
      REAL :: grad2v(0:nx,0:nny), grad4v(0:nx,0:nny)
      REAL :: taux(0:nnx,0:ny), tauy(0:nx,0:nny)
      REAL :: taux_steady(0:nx,0:nny), taux_var(0:nx,0:nny)
      REAL :: rhs_u(0:nnx,0:ny,nz), rhs_v(0:nx,0:nny,nz)
      REAL :: uu(0:nnx,0:ny), vv(0:nx,0:nny)
      REAL :: uu1(0:nnx,0:ny), vv1(0:nx,0:nny)
      REAL :: uu_old(0:nnx,0:ny), vv_old(0:nx,0:nny)
      REAL :: uh(0:nnx,0:ny), vh(0:nx,0:nny)
      REAL :: dissi_u(0:nnx,0:ny), dissi_v(0:nx,0:nny)
      REAL :: invLap_u(0:nnx,0:ny), invLap_v(0:nx,0:nny)
      REAL :: array_x(0:nnx,0:ny), array_y(0:nx,0:nny) ! dummy
      
      ! Centres/Centers :
      REAL :: eta(0:nx,0:ny,nz,3)
      REAL :: rhs_eta(0:nx,0:ny,nz)
      REAL :: div(0:nx,0:ny)
      REAL :: div1(0:nx,0:ny),div2(0:nx,0:ny)
      REAL :: divBT(0:nx,0:ny)
      REAL :: Psurf(0:nx,0:ny), rhs_Psurf(0:nx,0:ny)
      REAL :: B(0:nx,0:ny), B_nl(0:nx,0:ny), BS(0:nx,0:ny)
      REAL :: pressure(0:nx,0:ny), thickness(0:nx,0:ny)
      REAL :: eta_ag(0:nx,0:ny), eta_qg(0:nx,0:ny)
      REAL :: eta_ag_p(0:nx,0:ny,2)
      REAL :: p_out(0:nx,0:ny)
      REAL :: faces_array(0:nx,0:ny) ! dummy
      
      ! Noeuds/Nodes :
      REAL :: zeta(1:nx,1:ny)
      REAL :: zetaBT(1:nx,1:ny) ! fishpack
      REAL :: psiBT(1:nx,1:ny)  ! fishpack
      REAL :: zeta_G(0:nnx,0:nny,nz), zeta_AG(0:nnx,0:nny,nz) !AGdecomp.f90
      REAL :: q(1:nx,1:ny,nz), psi(1:nx,1:ny,nz) ! BCdecomp.f90
      REAL :: qmode(1:nx,1:ny,nz), psimode(1:nx,1:ny,nz) ! BCdecomp.f90
      REAL :: zeta_k(1:nx,1:ny,nz), zetamode(1:nx,1:ny,nz) ! BCdecomp.f90
      
      REAL :: array(0:nnx,0:nny) ! dummy
      

      
      !!! ----------- WW3 Coupling qties ----------- !!!
      ! Send/receive quantities
      REAL :: large_cur2WW3   (0:nxm1, 0:nym1, 2) ! Large current sent to WW3 (before bining)
      REAL :: large_WW3Ustokes(0:nxm1, 0:nym1, 2) ! Large current sent to WW3 (before bining)
      REAL :: large_WW3tauUst (0:nxm1, 0:nym1, 2) ! Large current sent to WW3 (before bining)
      REAL :: large_WW3tauDS  (0:nxm1, 0:nym1, 2) ! Large current sent to WW3 (before bining)
      REAL :: large_WW3tauIN  (0:nxm1, 0:nym1, 2) ! Large current sent to WW3 (before bining)
      REAL :: cur2WW3   (1:nxcou,1:nycou,2) ! Current sent to WW3.
      REAL :: WW3Ustokes(1:nxcou,1:nycou,2) ! Stokes' transport received.
      REAL :: WW3tauUst (1:nxcou,1:nycou,2) ! Directly received from WW3.
      REAL :: WW3tauDS  (1:nxcou,1:nycou,2) ! Directly received from WW3.
      REAL :: WW3tauIN  (1:nxcou,1:nycou,2) ! Directly received from WW3.

      ! Interpolated quantities from A-grid to C-grid.
      REAL :: taux_ust(0:nnx,0:ny),   tauy_ust(0:nx,0:nny)
      REAL :: taux_IN(0:nnx,0:ny),    tauy_IN(0:nx,0:nny)
      REAL :: taux_DS(0:nnx,0:ny),    tauy_DS(0:nx,0:nny)
      REAL :: taux_oc(0:nnx,0:ny,2),  tauy_oc(0:nx,0:nny,2) ! WW3 Coupling
      REAL :: UStokes(0:nnx,0:ny,2),  VStokes(0:nx,0:nny,2) ! WW3 Coupling
      REAL :: RHSu_SC(0:nnx,0:ny), RHSv_SC(0:nx,0:nny) ! WW3 Coupling (See rhs.f90 and dump_bin.f90)
      REAL :: RHSu_CL(0:nnx,0:ny), RHSv_CL(0:nx,0:nny) ! WW3 Coupling (See rhs.f90 and dump_bin.f90)
      REAL :: RHSu_BS(0:nnx,0:ny), RHSv_BS(0:nx,0:nny) ! WW3 Coupling (See rhs.f90 and dump_bin.f90)
      REAL :: rhsu_SW(0:nnx,0:ny,nz), rhsv_SW(0:nx,0:nny,nz) ! WW3 Coupling (See rhs.f90 and dump_bin.f90)
      !
      INTEGER :: mpi_grid_size
      real    :: HS(0:nx,0:ny) ! Stokes' drift thickness layer.
      
      
     !!! ---------- I/O qties definition ---------- !!!
      ! *** same as the old FFTW model 
      ! Sides :
      REAL :: u_out(1:szsubx,1:szsuby,nz)
      REAL :: uBT_out(1:szsubx,1:szsuby)
      REAL :: UStokes_out(1:szsubx,1:szsuby)
      REAL :: taux_ust_out(1:szsubx,1:szsuby)
      REAL :: taux_IN_out(1:szsubx,1:szsuby)
      REAL :: taux_DS_out(1:szsubx,1:szsuby)
      ! >
      REAL :: v_out(1:szsubx,1:szsuby,nz)
      REAL :: vBT_out(1:szsubx,1:szsuby)
      REAL :: VStokes_out(1:szsubx,1:szsuby)
      REAL :: tauy_ust_out(1:szsubx,1:szsuby)
      REAL :: tauy_IN_out(1:szsubx,1:szsuby)
      REAL :: tauy_DS_out(1:szsubx,1:szsuby)

      ! Centers :
      REAL :: eta_out(1:szsubx,1:szsuby,nz)
      REAL :: div_out(1:szsubx,1:szsuby)
      REAL :: divBT_out(1:szsubx,1:szsuby)
      REAL :: thickness_out(1:szsubx,1:szsuby)
      
      ! Nodes : 
      REAL :: zeta_out(1:szsubx,1:szsuby)
      REAL :: psiBT_out(1:szsubx,1:szsuby)
      REAL :: zetaBT_out(1:szsubx,1:szsuby)
      REAL :: psimode_out(1:szsubx,1:szsuby,nz)
      REAL :: zetamode_out(1:szsubx,1:szsuby,nz)

      
      !!! ---------- Other quantities ---------- 
      REAL :: sl, ed
      REAL :: ran2
      REAL :: mean_rhsuBT, mean_rhsvBT
      DOUBLE PRECISION :: ff(nx, ny)     ! Fishpack (Input field)
      DOUBLE PRECISION :: xa, xb, yc, yd ! Fishpack (Domain)
      DOUBLE PRECISION :: elmbda, pertrb ! Fishpack (Parameters)
      DOUBLE PRECISION :: bda(1), bdb(1), bdc(1), bdd(1) ! Fishpack (Newmann bndy)
      INTEGER          :: mbdcnd, nbdcnd ! Fishpack (boundary type)

      
      
      ! Baroclinic/Barotropic modes/solutions with LAPACK (see initialise.f90) : 
      REAL :: F_layer(1:nz,1:nz), A(1:nz,1:nz), A2(1:nz,1:nz), Fmodes(nz)
      REAL :: WI(1:nz), VL(1:nz,1:nz)
      REAL :: L2M(1:nz,1:nz),M2L(1:nz,1:nz)
      INTEGER WORKL(1:4*nz), INFO
      CHARACTER(8) :: ministr
      
      REAL :: forcing_qg(0:nnx,0:nny), forcing_ag(0:nnx,0:nny)
      REAL :: forcing_total(0:nnx,0:nny)
      REAL :: tmp(0:10)
      
      ! Initialize qties (see initialise.f90)
      real f(0:nny)
      real gprime(nz), Hsum, Hmin, H(nz), rho(nz+1) ! +1 because see initialise.f90
      real exp_coef, exp_amp, ratio
      real top(nz), bot(nz)
      real pdf(-100:100)
      real x, y, z, ramp, ramp0, today, cut_its, r, time

      ! real amp_matrix(864000) !3000 days
      real amp_forcing, amp, rms_amp, ampfactor, Omgrange !initialize_forcing
      real amp_matrix(nsteps+1),time_tmp(nsteps+1),amp_save,amp_load !nsteps+1 days
      real,allocatable:: amp_matrix_rand(:)
      real ke1, ke2, ke1_qg, ke2_qg, pe, pe_qg, etot, etot_qg
      real*4 tmp_out(10)

      ! Thickness correction (thickness_correction.f90)
      real hgap
      integer ijposition(2), ipos, jpos
      
      

      ! *** On bloque tout ce qui est en lien avec fftw parce que ça marchera pas, anyway.
      !2-D FFT spectra reduced
      !real*4 ke1_spec(0:nx), ke2_spec(0:nx)
      !real*4 ke1_qg_spec(0:nx), ke2_qg_spec(0:nx)
      !real*4 ke1_ag_spec(0:nx), ke2_ag_spec(0:nx)
      !real*4 for_to_spec(0:nx), for_ag_spec(0:nx)
      !1-D spectra
      !real*4 kex1_spec(0:nx/2), kex2_spec(0:nx/2)
      !real*4 key1_spec(0:ny/2), key2_spec(0:ny/2)
      !real*4 kex1_spec_tb(0:nx/2,1:ntsrow),key1_spec_tb(0:ny/2,1:ntsrow)
      !real*4 kex2_spec_tb(0:nx/2,1:ntsrow),key2_spec_tb(0:ny/2,1:ntsrow)
      !
      !double complex,dimension(nx/2+1,ny,nz) :: ufft,vfft,etafft
      !double complex,dimension(nx/2+1,ny) :: ftotalfft,fagfft,div_fft
      !
      !double complex,dimension(nx/2+1,ny,nz) :: u_agfft,v_agfft,u_qgfft,v_qgfft
      !double complex,dimension(nx/2+1,ny) :: u_agfft_bc,v_agfft_bc,u_gfft_bc,v_gfft_bc
      !double complex,dimension(nx/2+1,ny) :: eta_agfft,eta_qgfft

      ! Poincare modes
      !double complex,dimension(nx/2+1,ny,2) :: u_agfft_p,v_agfft_p,eta_agfft_p,div_fft_p 
      !real,dimension(nx/2+1,ny,2) :: omega_p ! omega field for poincaire plus and minus
      !real sgn1,tmp2,tmp3
      !double complex,dimension(nx/2+1,ny) :: kappa_ijsq,M !kappa**2 at (i,j) 
      integer i, j, k, ii, jj, kk, im, jp, jm, kp, km
      integer ip1, im1, jp1, jm1 
      integer ilevel, itt,  it, its, imode, ntimes, inkrow      
      real*4 tstime(1:ntsrow)
      
      !subsampling arrays (I/O)
      integer,allocatable:: isubx(:),isuby(:)!,iftsubkl(:,:)
      !integer rdsubk,rdsubl !temporary variables for reading (k,l) pair 
      
      !subsampling size (I/O)
      integer szftrdrow,szftrdcol
      
      !I/O info
      integer icount ,iftcount, count_specs_1, count_specs_2
      integer icount_srt,iftcount_srt
      integer count_specs_to, count_specs_ag

      character(2)  k_str
      character(88) scrap, which,which2,which3, spectbszx,spectbszy,specform, pathspects
      character(88) string1, string2, string3, string4, string5
      character(88) string1i, string2i, string3i, string4i, string5i
      character(88) string6, string7, string8, string9, string10
      character(88) string6i, string7i, string8i, string9i, string10i
      character(88) string11, string12, string13, string14, string15
      character(88) string11i, string12i, string13i, string14i, string15i
      character(88) string16, string17, string18, string19, string20
      character(88) string21, string22, string23, string24, string25
      character(88) string26, string27, string28, string29, string30
      character(88) string31, string32, string33, string34, string35
      character(88) string36, string37, string38, string39, string40
      character(88) string41, string42, string43, string44, string45
      character(88) string46, string47, string48, string49, string50
      character(88) string51, string52, string53, string54, string55
      character(88) string56, string57, string58, string59, string60
      character(88) string99,string98,fmtstr,fmtstr1

      INTEGER :: dummy_int
      REAL ::  dummy !Usep this for anything

     !!! Lowpass variables (see subs/Lowpass_filter.f90)
      REAL :: tcenter, tstart, tstop, param
      !! Filtered :   
      ! Curl
      REAL :: curlRHS_filtered(1:nx,1:ny)
      REAL :: curlRHS_BS_filtered(1:nx,1:ny)
      REAL :: curlRHS_CL_filtered(1:nx,1:ny)
      REAL :: curlRHS_SC_filtered(1:nx,1:ny)
      REAL :: zeta1_filtered(1:nx,1:ny)
      REAL :: curlTauUST_filtered(1:nx,1:ny)
      REAL :: curlTauIN_filtered(1:nx,1:ny)
      REAL :: curlTauDS_filtered(1:nx,1:ny)
      REAL :: curlUStokes_filtered(1:nx,1:ny)

      ! Div
      REAL :: divRHS_filtered(0:nx,0:ny)
      REAL :: divRHS_BS_filtered(0:nx,0:ny)
      REAL :: divRHS_CL_filtered(0:nx,0:ny)
      REAL :: divRHS_SC_filtered(0:nx,0:ny)
      REAL :: div1_filtered(0:nx,0:ny)
      REAL :: divTauIN_filtered(0:nx,0:ny)
      REAL :: divTauUST_filtered(0:nx,0:ny)
      REAL :: divTauDS_filtered(0:nx,0:ny)
      REAL :: divUStokes_filtered(0:nx,0:ny)

      !! Snap
      ! Curl
      REAL :: curlRHS_snap(1:nx,1:ny)
      REAL :: curlRHS_BS_snap(1:nx,1:ny)
      REAL :: curlRHS_CL_snap(1:nx,1:ny)
      REAL :: curlRHS_SC_snap(1:nx,1:ny)
      REAL :: zeta1_snap(1:nx,1:ny)
      REAL :: curlTauUST_snap(1:nx,1:ny)
      REAL :: curlTauIN_snap(1:nx,1:ny)
      REAL :: curlTauDS_snap(1:nx,1:ny)
      REAL :: curlUStokes_snap(1:nx,1:ny)

      ! Div
      REAL :: divRHS_snap(0:nx,0:ny)
      REAL :: divRHS_BS_snap(0:nx,0:ny)
      REAL :: divRHS_CL_snap(0:nx,0:ny)
      REAL :: divRHS_SC_snap(0:nx,0:ny)
      REAL :: div1_snap(0:nx,0:ny)
      REAL :: divTauIN_snap(0:nx,0:ny)
      REAL :: divTauUST_snap(0:nx,0:ny)
      REAL :: divTauDS_snap(0:nx,0:ny)
      REAL :: divUStokes_snap(0:nx,0:ny)

     ! ==================================== ! 

      
      
      !include 'fftw_stuff/fft_params.f90'
      !include 'fftw_stuff/fft_init.f90'
      
      ! >>> Modification CEL >>>
      IF (cou) THEN
         mpi_grid_size = nxcou*nycou
         
         !!!  --- Starting MPI --- !!!
         CALL MPI_INIT(ierror)
         CALL MPI_COMM_RANK(MPI_COMM_WORLD, procid, ierror)
         CALL MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierror)
         PRINT *, "SW  (COMM_WORLD) : Je suis le proc :", procid, "sur", numprocs
         
         !!! --- Splitting MPI because WW3 is dumb
         CALL MPI_COMM_SPLIT( MPI_COMM_WORLD, 2, 0, MPI_SECOND, ierror)
         CALL MPI_COMM_RANK(MPI_SECOND, procid_sec, ierror)
         CALL MPI_COMM_SIZE(MPI_SECOND, numprocs_sec, ierror)
         PRINT *, "SW (COMM_SECOND) : Je suis le proc :", procid_sec, "sur", numprocs_sec
      END IF
      ! <<< Modification CEL (END) <<<

      
      ! === Complex number i
       eye=(0.0,1.0)

      ! === Set random number seeds
      call date_and_time(VALUES=values)
      iseed = -(values(8)+values(7)+values(6))

      ! === Allocate variables
      write(*,*) 'Physical field subsmpling res.',szsubx,'by',szsuby
      ! szftrdrow=ceiling((nx/2+1)/(ftsubsmprto+1e-15))
      ! szftrdcol=ceiling(ny/(ftsubsmprto+1e-15))
      allocate(isubx(szsubx),isuby(szsuby))
      ! === Define subsampling range in spatial space
      isubx=(/(i, i=1,nx, subsmprto)/)
      isuby=(/(i, i=1,ny, subsmprto)/)
      !!!do i=1,size(isubx)
      !!!write(*,*) 'isubx',i,'sub-x indices',isubx(i)
      !!!end do
      !!!do j=1,size(isuby)
      !!!write(*,*) 'isuby',j,'sub-y indices',isuby(j)
      !!!end do
      write(*,*) 'Subsampled file size: 4X',size(isubx),'X',size(isubx),'X',nz


      ! === FFT subsmpling, fed with (k,l) pair
      ! include 'subs/read_kxky_subsmp.f90'

      ! === 
      !set icount and time based on if restart is true or false !moved to initialize
      !  icount = 0 !for  output file index
      !  time = 0.  !in second
      write(*,*) 'total step, nstep= ', nsteps
      write(*,*) 'write output data for every',iout, 'steps, for ndays=', ndays
      write(*,*) 'one day = ', 86400/dt, 'time steps'
      write(*,*) 'Spectrum time series has', ntsrow, 'lines'
      write(*,*) 'dx = ', dx
      write(*,*) 'CFL = ', dx/dt
      
      ! --- Initializing each fields and rossby radii
      include 'subs/initialize.f90'

      its=1
      itlocal=1
      ! (***) Faut le mettre si on veut un forçage transient.
      !include 'subs/initialize_forcing.f90'


      if(restart .eqv. .false.)   nspecfile=0
      write(*,*) 'iout,ispechst,ntsrow',iout,ispechst,ntsrow
      
      !call get_taux(taux_steady,amp_matrix(its),taux)

      ! Lowpass filter (see subs/Lowpass_filter.f90)
      if (cou) then
         ! 100 days after coupling
         tstart  = 100.*86400.
      else
         ! Two years after begining : 
         tstart  = 2.*365.*86400.
      endif
      tstop   = tstart + 8.*86400.
      tcenter = (tstart + tstop)/2
      ! Gaussian parameter for each timestep.
      param = (time-tcenter)/1.1/86400.
      param = param**2
      !normalization 2304 for dt = 300, window = 8 days width = 1.1


      
      !==============================================================
      !
      !      1st time step
      !
      !==============================================================
      write(*,*) 'First time step'
      
      !ramp = 1.
      !if ( use_ramp .eqv. .true. ) then
      !   tmp(1) = 10*twopi/f0/dt
      !   ramp0  =  1./tmp(1)
      !   ramp = ramp0
      !endif


      !
      ! >>> Modification CEL >>>
      ! --- FIRST MPI CALL HERE. 
      ! Receiving first  Wavewatch atmospheric stresses here.
      IF (cou) THEN
      ! --- Coupling
         !its = 1
         uu(:,:) = u(:,:,1,1) 
         vv(:,:) = v(:,:,1,1)
         include 'subs/coupling_ww3.f90'
      
      ! ... then forgetting them to keep our restart files. 
         UStokes(:,:,2) = UStokes(:,:,1) 
         VStokes(:,:,2) = VStokes(:,:,1) 
         taux_oc(:,:,2) = taux_oc(:,:,1)
         tauy_oc(:,:,2) = tauy_oc(:,:,1)
      END IF
      ! <<< Modification CEL (END) <<<.
      !

      pressure(:,:) =  0.
      do k = 1,nz
         uu(:,:) = u(:,:,k,1) 
         vv(:,:) = v(:,:,k,1)
         uu_old(:,:) = u(:,:,k,1) 
         vv_old(:,:) = v(:,:,k,1)

         ! Adding random noise from a psiBT
         if (restart .eqv. .false.) then
            psiBT(:,:) = 0.
            !CALL RANDOM_NUMBER(psiBT(2:nx-1,2:ny-1))
            do j = 1,ny
               jm=j-1
            do i = 1,nx
               im=i-1
               psiBT(i,j) = SIN(2*twopi*im/(nx-1))*SIN(2*twopi*jm/(ny-1))/1e10
            enddo
            enddo

            psiBT(1,:)  = 0.
            psiBT(nx,:) = 0.
            psiBT(:,1)  = 0.
            psiBT(:,ny) = 0.
            
            do j = 1,ny-1
            do i = 1,nx-1
               uu(i,j) =  - (psiBT(i,j+1) - psiBT(i,j))/dy  ! barotropic part-x
               vv(i,j) =    (psiBT(i+1,j) - psiBT(i,j))/dx  ! barotropic part-y
            enddo
            enddo
         endif
         !

         ! Ramp :
         ! ramp0 is the slope of the ramp. Such that ramp = ramp0*its
         ramp = 1
         cut_its = cut_days*86400/dt
         if (use_ramp) then
            ramp0 = real(dt/(31*86400)) ! Slope per iteration (extended on 31 days).
            ramp = max(0.,ramp0*(float(its)-cut_its))
         endif

         ! Pas besoin d'update les contours de zeta ou psi.
         do j = 2,ny-1
            jp = j+1
            jm = j-1
         do i = 2,nx-1
            ip1 = i+1
            im = i-1
            
            zetaBT(i,j) = (psiBT(ip1,j)+psiBT(im,j)-2.*psiBT(i,j))/dx/dx   &
            &             + (psiBT(i,jp)+psiBT(i,jm)-2.*psiBT(i,j))/dy/dy
            
         enddo
         enddo

         
         ! Boundaries on the random noise
         array_x = uu
         array_y = vv
         include 'subs/no_normal_flow.f90'
         include 'subs/free_or_partial_slip.f90'
         uu = array_x
         vv = array_y

         



         ! Finding thickness locally
         if (k.eq.1) then
            thickness(:,:) = H(k) - eta(:,:,k+1,1) 
         else if (k.eq.nz) then
            pressure(:,:) = pressure(:,:) + rho(1)*gprime(k)*eta(:,:,k,1) 
            thickness(:,:) = H(k) + eta(:,:,k,1)   
         else
            pressure(:,:) =  pressure(:,:) + rho(1)*gprime(k)*eta(:,:,k,1) 
            thickness(:,:) =  H(k) + eta(:,:,k,1) - eta(:,:,k+1,1)
         endif
         include 'subs/rhs.f90'
      enddo  ! end of the do k = 1,nz loop


      ! ================================================== !
      !                  RHS (begining)                    !
      ! ================================================== !

      ! eta-loop : Starts from the bottom, because RHS eta_k = RHS h_k + RHS eta_k-1
      rhs_eta(:,:,1) = 0. ! First layer is a rigid lid.
      faces_array(:,:) = rhs_eta(:,:,nz) ! bottom layer. faces_array = \Delta \eta_{k+1}
      eta(:,:,nz,2) = eta(:,:,nz,1) + dt*faces_array(:,:)
      do k = nz-1, 2, -1
         faces_array(:,:)   = rhs_eta(:,:,k) + faces_array(:,:)
         eta(:,:,k,2) = eta(:,:,k,1) + dt*faces_array(:,:)
      end do ! end k-loop
      
      u(:,:,:,2) = u(:,:,:,1) + dt*rhs_u(:,:,:)
      v(:,:,:,2) = v(:,:,:,1) + dt*rhs_v(:,:,:)

      
      ! ================================================== !
      !                      RHS (END)                     !
      ! ================================================== !

      
      ! -------------------------------------------------- !
      !        ->>>- Barotropic psi-correction ->>>-       !
      ! -------------------------------------------------- !
        !
        !     need to correct RHS_u,v for surface pressure 
        !
      write(*,*) 'First barotropic psi-correction'
      ilevel = 2
      p_out(:,:) = 0.
      ! include 'subs/psiBT_correction.f90'
      include 'subs/init_fishpack.f90'
      include 'subs/psiBT_correction_fishpack.f90'

      ! -------------------------------------------------- !
      !    -<<<- Barotropic psi-correction (End) -<<<-     !
      ! -------------------------------------------------- !

      ! --- updating time ---
      time = dt
      its = its + 1
      !call get_taux(taux_steady,amp_matrix(its),taux)

      ! --- Updating fields indicators ---
      u(:,:,:,3) = u(:,:,:,2)
      v(:,:,:,3) = v(:,:,:,2)
      eta(:,:,2:nz,3) = eta(:,:,2:nz,2)
      ! eta1 is barotropic streamfunction (psiBT).
      eta(1:nx,1:ny,1,3) = PsiBT(1:nx,1:ny)

      ! --- Print output --- 
      include 'subs/dump_bin.f90'
      !


      
      !FFTW   !include 'subs/IOheader.f90' 
      !==============================================================
      !
      !      subsequent time steps
      !
      !==============================================================

      do its = 2, nsteps
         ! determine local time step !tstime
         itlocal=mod(its,itape)       ! for 1-D time series, relative its for itape
         itsrow=int(itlocal/ispechst) ! xxth row in a 1-D time series file
         ! forced to set 0-mod to the last row
         if (itsrow==0.and.its.gt.ispechst) itsrow=ntsrow

         if(itlocal==0) then
            itlocal = itape       ! for 0-mod
            nspecfile=nspecfile+1 ! increase the file index by 1
         end if

         if(mod(itlocal,ispechst)==0) then
            tstime(itsrow)=time/86400. !record each time when print out 
            ! write(*,*) 'its,itsrow,time[day]',its,itsrow,tstime(itsrow)
         end if

         !        =================
         if ( mod(its,iout*10).eq.0 ) then
         write(*,*) 'current time step, iout', its, iout
         !shuffle forcing terms
         ! call shuffle_phi
         end if
         !        =================
         ramp = 1
         if (use_ramp) then
            if (its.le.(cut_its)) then 
               ramp = max(0.,ramp0*(float(its) - cut_its))
            else
               ramp = min(1.,ramp0*(float(its) - cut_its))
            endif
         endif
         !        =================



         !
         ! --- SUBSEQUENT MPI CALL HERE ---
         ! Fetching coupled fields
         IF (cou) THEN
            uu(:,:) = u(:,:,1,2) 
            vv(:,:) = v(:,:,1,2)
            include 'subs/coupling_ww3.f90'
         END IF
         !
         

         ! --- Thickness correction (begin) ---
         !  n.b. this mecanism balance isopicnals of the model before it explodes
         !       or before layers get null thicknesses.
         !
         if (.true.) then
            do k = 1,nz-1
               ! 1. Finding thickness of current layer. 
               if (k.eq.1) then
                  thickness(:,:) = H(k) - eta(:,:,k+1,2) 
               else if (k.eq.nz) then
                  thickness(:,:) = H(k) + eta(:,:,k,2)
               else
                  thickness(:,:) = H(k) + eta(:,:,k,2) - eta(:,:,k+1,2)
               end if
            
               ! 2. Finding thickness ratio 
               tmp(1) = minval(thickness)/H(k)

               ! 3. Threshold with while loop
               do while (tmp(1).lt.0.02) ! 2%
                  !
                  print *, " ------------------------------------------- "
                  print *, " > Thickness before correction :: its=", its
                  print *, "  -> minval(h) = ", minval(thickness)
                  print *, "  -> Ratio = minval(h)/H(k) = ", tmp(1)
                  print *, "  -> H(k) =", H(k)
                  !
                  include 'subs/thickness_correction.f90'
                  tmp(1) = minval(thickness)/H(k)
                  !
                  print *, " > Thickness after correction :"
                  print *, "  -> minval(h) = ", minval(thickness)
                  print *, "  -> Ratio = minval(h)/H(k) = ", tmp(1)
                  ! 
               enddo
            enddo  ! k
         endif
         ! Thickness correction (End)
         !
         
         pressure(:,:) =  0.
         do k = 1,nz
            uu(:,:) = u(:,:,k,2)
            vv(:,:) = v(:,:,k,2)
            uu_old(:,:) = u(:,:,k,1)
            vv_old(:,:) = v(:,:,k,1)

            if (k.eq.1) then
               ! p = 0, because eta(:,:,1,:) = 0.
               thickness(:,:) = H(k) - eta(:,:,k+1,2) 
            else if (k.eq.nz) then
               pressure(:,:)  = pressure(:,:) + rho(1)*gprime(k)*eta(:,:,k,2) 
               thickness(:,:) = H(k) + eta(:,:,k,2)
            else
               pressure(:,:)  = pressure(:,:) + rho(1)*gprime(k)*eta(:,:,k,2) 
               thickness(:,:) = H(k) + eta(:,:,k,2) - eta(:,:,k+1,2)
            end if
            include 'subs/rhs.f90'
         enddo  ! k
         
         
      ! ================================================== !
      !                  RHS (begining)                    !
      ! ================================================== !
         !
         ! -> We apply Leapfrog timestep qty(3) = qty(1) + RHS_qty(2)*dt
         !
         ! eta-loop : Starts from the bottom, because RHS eta_k = RHS h_k + RHS eta_k-1
         rhs_eta(:,:,1) = 0. ! Rigid lid
         faces_array(:,:) = rhs_eta(:,:,nz)
         eta(:,:,nz,3) = eta(:,:,nz,1) + 2.*dt*faces_array(:,:)
         do k = nz-1, 2, -1
            faces_array(:,:) = rhs_eta(:,:,k) + faces_array(:,:)
            eta(:,:,k,3) = eta(:,:,k,1) + 2.*dt*faces_array(:,:)
         end do ! end k-loop
         
         u(:,:,:,3) = u(:,:,:,1) + 2.*dt*rhs_u(:,:,:)
         v(:,:,:,3) = v(:,:,:,1) + 2.*dt*rhs_v(:,:,:)
         
      ! ================================================== !
      !                      RHS (End)                     !
      ! ================================================== !


         
      ! -------------------------------------------------- !
      !        ->>>- Barotropic psi-correction ->>>-       !
      ! -------------------------------------------------- !
         !
         !     stuff for barotropic psi-correction
         !
         ilevel = 3
         p_out(:,:) = 0.
         include 'subs/psiBT_correction_fishpack.f90'
         
      ! -------------------------------------------------- !
      !    -<<<- Barotropic psi-correction (End) -<<<-     !
      ! -------------------------------------------------- !

         
         ! ->- Lowpass filter output ->-
         include "subs/Lowpass_filter.f90"
         
         ! --- Updating time ---
         time = time + dt
         today = time/86400.
         !call get_taux(taux_steady,amp_matrix(its),taux)

         ! --- robert filter ---
         u(:,:,:,2) = u(:,:,:,2) + rf*(u(:,:,:,1)+u(:,:,:,3)-2*u(:,:,:,2))
         v(:,:,:,2) = v(:,:,:,2) + rf*(v(:,:,:,1)+v(:,:,:,3)-2*v(:,:,:,2))
         eta(:,:,:,2) = eta(:,:,:,2)   &
            &        + rf*(eta(:,:,:,1)+eta(:,:,:,3)-2*eta(:,:,:,2))
         u(:,:,:,1) = u(:,:,:,2)
         v(:,:,:,1) = v(:,:,:,2)
         eta(:,:,:,1) = eta(:,:,:,2)
         u(:,:,:,2) = u(:,:,:,3)
         v(:,:,:,2) = v(:,:,:,3)
         eta(:,:,:,2) = eta(:,:,:,3)

         ! --- mean kinetic energy diagnostics --- 
         ke1 = 0.
         ke2 = 0.
         do j = 1, ny-1
         do i = 1, nx-1
            ke1 = ke1+(u(i,j,1,3)**2+v(i,j,1,3)**2)/2.
            ke2 = ke2+(u(i,j,2,3)**2+v(i,j,2,3)**2)/2.
         enddo
         enddo      

         write(300,'(i8,1f12.4,3e12.4)') its, time/86400.,taux(nx/2,ny/2), ke1/nx/ny, ke2/nx/ny
         call flush(300)

         !if(nsteps.lt.1.and.save_movie) then
         !   if ( mod(its,iout).eq.0 ) then  ! output 
         !      include 'subs/div_vort.f90'
         !      include 'subs/dump_gnu1a.f90'
         !   end if
         !end if


         !start of diognostics en savefiles.
         
         ! Calculating diognostic only when outputting physical or Fourier fields
         if(mod(its,ispechst).eq.0.or.mod(its,iout).eq.0) then 
            !  include 'subs/div_vort.f90' 
            !  include 'subs/tmp_complex.f90'
            !  include 'subs/calc_q.f90'
            include 'subs/diags.f90'
         end if

         ! Printing/saving fields in /data/. 
         if (save_movie.and. mod(its,iout).eq.0 ) then  ! output 
            icount = icount + 1

            ! *** Eta(z=1) is the barotropic streamfunction (since eta(0) = 0 because of rigid lid)
            eta(1:nx,1:ny,1,3) = PsiBT(1:nx,1:ny)

            include 'subs/dump_bin.f90'

            if(mod(its,iout).eq.0)write(*,*) 'current its',its
            print*, 'writing data No.', icount
         end if

         !!! FFTW stuff that CÉ removed.
         !if ( mod(its,ispechst).eq.0 ) then
         !   count_specs_1 = count_specs_1 + 1
         !   count_specs_2 = count_specs_2 + 1 
         !   count_specs_to = count_specs_to + 1
         !   count_specs_ag = count_specs_ag + 1
         !
         !   include 'subs/calc_1Dspec.f90'
         !   include 'subs/calc_2Dspec.f90'              
         !   
         !   
         !   ! AG+, AG- decomposition
         !   ! 1. convert eta_A and div to FFT space
         !   datr(:,:) = eta_ag(1:nx,1:ny)
         !   include 'fftw_stuff/spec1.f90'
         !   eta_agfft(:,:)=datc ! BC mode
         !
         !   k=1
         !   include 'subs/div_vort.f90'
         !   div1 = div
         !
         !   k=2
         !   include 'subs/div_vort.f90'
         !   div2 = div
         !
         !   datr(:,:) = div2(1:nx,1:ny)-div1(1:nx,1:ny)
         !   include 'fftw_stuff/spec1.f90'
         !   div_fft(:,:)=datc ! BC mode
         !
         !   ! for each k,l pair, calculate M, then eta_+ and eta_-
         !   ! In polarization relations, (2*pi) in FT is vanished
         !   kappa_ijsq=nkx**2+nky**2 !nkx,nky defined in fft_params.f90
         !   M=eye*sqrt(c_bc**2*kappa_ijsq+f0**2)*gprime(2)/c_bc**2 !only bc mode
         !   ! two AG frequencies omega_+ and omega_- (BC only for rigid lid)
         !   do imode = 1,2
         !      sgn1=-(-1.)**imode; ! sgn1=1 when imode ==1, sgn1=-1 when imode==2
         !      omega_p(:,:,imode)=sgn1*sqrt(c_bc**2*kappa_ijsq+f0**2)
         !      eta_agfft_p(:,:,imode)=1/(2*M)*(M*eta_agfft+sgn1*div_fft)
         !      u_agfft_p(:,:,imode)=gprime(2)*eta_agfft(:,:)*(omega_p(:,:,imode)*nkx+eye*f0*nky)&
         !           /(c_bc**2*kappa_ijsq)
         !      v_agfft_p(:,:,imode)=gprime(2)*eta_agfft(:,:)*(omega_p(:,:,imode)*nky-eye*f0*nkx)&
         !           /(c_bc**2*kappa_ijsq)
         !   end do ! imode
         !   ! AG+-mode decomp finished
         !
         !   ! Write 2-D FFT fields
         !   if(save2dfft) then
         !      iftcount=iftcount+1
         !      include 'subs/dump_bin_spec2d.f90'
         !   endif !itsrow==ntsrow
         !endif


        ! writing back-up restarting files for
        if(mod(its,int(daysperrestart*86400/dt))==0) then
                include 'subs/writing_restart_files.f90'
        end if
        if(its==nsteps) then
         include 'subs/writing_restart_files.f90'
      end if
      
     enddo ! its
     !===== time loop ends here

     ! MPI-COUPLING
     if (cou) then 
        CALL MPI_FINALIZE(ierror)
     endif
     ! MPI-COUPLING
     !include 'fftw_stuff/fft_destroy.f90'
    end program main




    
function ran2(idum)
   use data_initial
   integer :: idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
   real :: ran2,AM,EPS,RNMX
   PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,&
         IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,&
         IR2=3791,NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
   INTEGER :: idum2,j1,k1,iv(NTAB),iy
   SAVE iv,iy,idum2
   DATA idum2/123456789/,iv/NTAB*0/,iy/0/

   if (idum .le. 0) then
         idum=max(-idum,1)
         idum2 = idum
         do j1 = NTAB+8,1,-1
            k1=idum/IQ1
            idum=IA1*(idum-k1*IQ1)-k1*IR1
            if (idum .lt. 0) idum=idum+IM1
            if (j1 .le. NTAB) iv(j1) = idum
         end do
         iy=iv(1)
   end if
   k1=idum/IQ1
   idum=IA1*(idum-k1*IQ1)-k1*IR1
   if (idum .lt. 0) idum=idum+IM1
   k1=idum2/IQ2
   idum2=IA2*(idum2-k1*IQ2)-k1*IR2
   if (idum2 .lt. 0) idum2=idum2+IM2
   j1 = 1+iy/NDIV
   iy = iv(j1) - idum2
   iv(j1) = idum
   if (iy .lt. 1) iy = iy+IMM1
   ran2=min(AM*iy,RNMX)
   return
end function ran2

function ran1(idum)
   !
   ! ----------- stolen from numerical recipes,p. 271
   !
   integer idum, ia, im, iq, ir, ntab, ndiv
   real ran1, am, eps, rnmx
   parameter (ia=16807,im=2147483647,am=1.0/im,iq=127773,ir=2836.0,&
               ntab=32,ndiv=1+(im-1)/ntab,eps=1.2e-07,rnmx=1.0-eps)
   integer j1, k1, iv(ntab), iy
   save iv, iy
   data iv /ntab*0/, iy /0/
   if(idum .le. 0 .or. iy .eq. 0) then
      idum = max(-idum,1)
      do j1=ntab+8,1,-1
         k1 = idum/iq
         idum = ia*(idum - k1*iq) - ir*k1
         if(idum .lt. 0) idum = idum + im
         if(j1 .le. ntab) iv(j1) = idum
      enddo
      iy = iv(1)
   endif
   k1     = idum/iq
   idum  = ia*(idum - k1*iq) - ir*k1
   if(idum .lt. 0) idum = idum + im
   j1     = 1 + iy/ndiv
   iy    = iv(j1)
   iv(j1) = idum
   ran1  = min(am*iy, rnmx)
   !
   return
end function ran1


subroutine get_tau_amp_AR(time1,amp1)
   use data_initial
   real,intent(in) :: time1
   real,intent(out) :: amp1
   real amp_forcing,ran2
   allocate(phi(n_omega))
   amp1=0.0
   amp1=sum(A_n*sin(omega*time1 + phi))
end subroutine

subroutine shuffle_phi
   use data_initial
   real ran2
   integer iw
   do iw = 1, n_omega !1e*5
      phi(iw)= ran2(iseed)*2.0*pi
   end do
   write(*,*) 'phi for stochastic forcing is shuffled from 1 to',n_omega
end subroutine 

subroutine Euler_Maruyama(amp1,amp2)
   use data_initial
   real ran2,volsqrdt
   integer i,j
   ! set dt to be unity - could add as input...
   ! set in data_initial's parameters.f90
   
   ! pre set vol*sqrt(dt) vol to reduce run time
   volsqrdt = vol*sqrt(dt);
   amp2 = amp1 + theta*(mu - s0)*dt + ran2(iseed)*volsqrdt;
end subroutine

FUNCTION gasdev(idum) 
   INTEGER idum
   REAL gasdev
   !Returns a normally distributed deviate with zero mean and unit variance, using ran1(idum)
   !as the source of uniform deviates.
   INTEGER iset
   REAL fac,gset,rsq,v1,v2,ran1 
   SAVE iset,gset
   DATA iset/0/
   if (idum.lt.0) iset=0
   if (iset.eq.0) then
1     v1=2.*ran1(idum)-1. 
      v2=2.*ran1(idum)-1.
      rsq=v1**2+v2**2 
      if(rsq.ge.1..or.rsq.eq.0.) goto 1 
      fac=sqrt(-2.*log(rsq)/rsq) 
      gset=v1*fac
      gasdev=v2*fac
      iset=1
   else 
      gasdev=gset
      iset=0 
   endif
   return 
 END FUNCTION gasdev

 SUBROUTINE get_taux(taux_steady_in,amp_in,taux_out)
   use data_initial
   implicit none
   real,intent(in):: taux_steady_in(0:nnx,0:nny),amp_in
   real,intent(out)::taux_out(0:nnx,0:nny)
   if (forcingtype==0) then
      taux_out(:,:) = taux_steady_in(:,:) +  c_tauvar*tau0*amp_in
   else if(forcingtype==1) then
      taux_out(:,:) = taux_steady_in(:,:)*(1.+amp_in)
   endif
 END SUBROUTINE get_taux
    
 !SUBROUTINE interp_matrix(matin,matout,dirx,diry,irrt)
 !  use data_initial
 !  real,dimension(0:nnx,0:nny),intent(in) :: matin
 !  real,dimension(0:nnx,0:nny),intent(out) :: matout
 !  integer,intent(in) :: irrt,dirx,diry
 !  real array(0:nnx,0:nny),array1(0:nnx,0:nny),array2(0:nnx,0:nny)
 !  integer ii,i,j,ip1,jp
 !  if (abs(dirx)<=1.and.abs(diry)<=1) then
 !     array1=matin
 !     array2=0.
 !     do ii=1,irrt
 !        if(dirx*diry.eq.0) then ! 1-D interpolate
 !           array(1:nx,1:ny)= 0.5*(array1(1:nx,1:ny)+array1(1+dirx:nx+dirx,1+diry:ny+diry))
 !           include 'subs/bndy.f90'
 !           array2=array2+array
 !           ! inferred v at eta-grid from v-grid (going upward)
 !           array(1:nx,1:ny)= 0.5*(array1(1:nx,1:ny)+array1(1-dirx:nx-dirx,1-diry:ny-diry))
 !           include 'subs/bndy.f90'
 !        else if (dirx*diry.ne.0) then
 !           array(1:nx,1:ny)= 0.25*(array1(1:nx,1:ny)+array1(1+dirx:nx+dirx,1:ny)+ &
 !                &              array1(1:nx,1+diry:ny+diry)+array1(1+dirx:nx+dirx,1+diry:ny+diry))
 !           include 'subs/bndy.f90'
 !           array2=array2+array
 !           array(1:nx,1:ny)= 0.25*(array1(1:nx,1:ny)+array1(1-dirx:nx-dirx,1:ny)+ &
 !                &              array1(1:nx,1-diry:ny-diry)+array1(1-dirx:nx-dirx,1-diry:ny-diry))
 !           include 'subs/bndy.f90'
 !        end if
 !        ! get the residual
 !        array1 = matin - array
 !     end do
 !  else
 !     write(*,*) 'Interpolation is set wrong'
 !  end if
 !  matout=array2
 !
 !END SUBROUTINE interp_matrix


 
 SUBROUTINE geometric_interpolation(mat_in, mat_out, n_in, n_out)
   ! This subroutine takes a small matrix (mat_in) and geometricly interpolate an
   ! a big matrix out of it. Interpolation is done with stencil ratio over the
   ! small grid. 
   integer, intent(in) :: n_in, n_out
   real, dimension(n_in,n_in),   intent(in)  :: mat_in
   real, dimension(n_out,n_out), intent(out) :: mat_out
   real, dimension(n_out,n_out) :: mat_mid
   real, allocatable :: weight(:,:), stencil_mat(:,:)
   integer :: iin, jin, iout, jout, i, j
   integer :: ips, ims, jps, jms
   integer :: ratio, stencil_size, sarea ! Ratio is also the stencil size!
   integer :: i0,j0
   real    :: delta, square_sum
   
   ! 1. Finding stencil size. 
   delta = (real(n_out)/n_in)
   if (mod(delta,1.).lt.0.001) then
      ratio = nint(delta)
      sarea = ratio**2 ! Stencil area [squares]
   else
      print *,"geometric interpolation :: Size problem. n_out must be a multiple of n_in"
      stop
   end if

   
   ! 2. i0 and j0 are stencil radii in terme of indices.
   if (mod(ratio,2) .eq. 0) then
      ! If ratio is even :
      stencil_size = ratio+1      
      i0 = ratio/2
      j0 = i0
      allocate(stencil_mat(stencil_size,stencil_size))
      allocate(weight(stencil_size,stencil_size))
      weight(:,:) = 1.
      weight(1,:)            = 0.5*weight(1,:)
      weight(stencil_size,:) = 0.5*weight(stencil_size,:)
      weight(:,1)            = 0.5*weight(:,1)
      weight(:,stencil_size) = 0.5*weight(:,stencil_size)
   else
      ! If ratio is odd (way easier) :
      stencil_size = ratio
      i0 = (ratio+1)/2
      j0 = i0
      allocate(stencil_mat(stencil_size,stencil_size))
      allocate(weight(stencil_size,stencil_size))
      weight(:,:) = 1.
   endif

   
   ! 3. Filling intermediate matrix : mat_mid
   ! N.B. now general for odd and even ratios.
   do jin = 1,n_in
      jout = 1 + (jin-1)*ratio   
   do iin = 1,n_in
      iout = 1 + (iin-1)*ratio
      ! Creating square bins on mat_mid.
      mat_mid(iout:ratio*iin,jout:ratio*jin) = mat_in(iin,jin)
   enddo
   enddo

   
   ! 4. Applying interpolation stencil.
   mat_out = mat_mid
   do jout = 1+j0,n_out-j0
      jms = jout - j0
      jps = jout + j0
   do iout = 1+i0,n_out-i0
      ims = iout - i0
      ips = iout + i0
      
      ! Summing values on the stencil.
      square_sum = 0.
      stencil_mat = mat_mid(ims:ips, jms:jps)
      do j = 1,stencil_size
      do i = 1,stencil_size
         square_sum = square_sum + weight(i,j)*stencil_mat(i,j)
      enddo
      enddo
      mat_out(iout,jout) = square_sum/sarea
      
   enddo
   enddo

   ! 5. Boundaries
   mat_out(1,:) = mat_out(2,:)/2
   mat_out(:,1) = mat_out(:,2)/2
   mat_out(:,n_out) = mat_out(:,n_out-1)/2
   mat_out(n_out,:) = mat_out(n_out-1,:)/2
   

 END SUBROUTINE geometric_interpolation
