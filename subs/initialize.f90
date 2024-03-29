

     ! Coupling quantities (WAVEWATCH III) >>>
       taux_oc(:,:,:) = 0.
       tauy_oc(:,:,:) = 0.
       taux_ust(:,:) = 0.
       tauy_ust(:,:) = 0.
       taux_IN(:,:) = 0.
       tauy_IN(:,:) = 0.
       taux_DS(:,:) = 0.
       tauy_DS(:,:) = 0.
       UStokes(:,:,:) = 0.
       VStokes(:,:,:) = 0.

       
       ! Initialising Coupled RHS as nul.
       ! N.B. Will stay unmodified is cou = .false.
       RHSu_SC(:,:) = 0.
       RHSv_SC(:,:) = 0.
       RHSu_CL(:,:) = 0.
       RHSv_CL(:,:) = 0.
       RHSu_BS(:,:) = 0.
       RHSv_BS(:,:) = 0.
       ! This on will always be modified (see rhs.f90)
       rhsu_SW(:,:,:) = 0.
       rhsv_SW(:,:,:) = 0.

       
     ! Lowpass_filter (see subs/Lowpass_filter.f90)
      ! Snap
      ! curl
       curlRHS_snap(:,:)     = 0.
       curlRHS_BS_snap(:,:)  = 0.
       curlRHS_CL_snap(:,:)  = 0.
       curlRHS_SC_snap(:,:)  = 0.
       zeta1_snap(:,:)       = 0.
       curlTauUST_snap(:,:)  = 0.
       curlTauIN_snap(:,:)   = 0.
       curlTauDS_snap(:,:)   = 0.
       curlUStokes_snap(:,:) = 0.   
      ! Div
       divRHS_snap(:,:)     = 0.
       divRHS_BS_snap(:,:)  = 0.
       divRHS_CL_snap(:,:)  = 0.
       divRHS_SC_snap(:,:)  = 0.
       div1_snap(:,:)       = 0.
       divTauIN_snap(:,:)   = 0.
       divTauUST_snap(:,:)  = 0.
       divTauDS_snap(:,:)   = 0.
       divUStokes_snap(:,:) = 0.
       
     ! Filtered
      ! curl
       curlRHS_filtered(:,:)     = 0.
       curlRHS_BS_filtered(:,:)  = 0.
       curlRHS_CL_filtered(:,:)  = 0.
       curlRHS_SC_filtered(:,:)  = 0.
       zeta1_filtered(:,:)       = 0.
       curlTauUST_filtered(:,:)  = 0.
       curlTauIN_filtered(:,:)   = 0.
       curlTauDS_filtered(:,:)   = 0.
       curlUStokes_filtered(:,:) = 0.   
      ! Div
       divRHS_filtered(:,:)     = 0.
       divRHS_BS_filtered(:,:)  = 0.
       divRHS_CL_filtered(:,:)  = 0.
       divRHS_SC_filtered(:,:)  = 0.
       div1_filtered(:,:)       = 0.
       divTauIN_filtered(:,:)   = 0.
       divTauUST_filtered(:,:)  = 0.
       divTauDS_filtered(:,:)   = 0.
       divUStokes_filtered(:,:) = 0. 

       
     ! Model arrays >>>
       u(:,:,:,:) = 0.
       v(:,:,:,:) = 0.
       eta(:,:,:,:) = 0.
       u_ag(:,:,:) = 0.
       v_ag(:,:,:) = 0.
       eta_ag(:,:) = 0.
       u_qg(:,:,:) = 0.
       v_qg(:,:,:) = 0.
       eta_qg(:,:) = 0.
       taux(:,:) = 0.
       tauy(:,:) = 0.
       zeta(:,:) = 0.
       div(:,:) = 0.
       B(:,:) = 0.
       BS(:,:) = 0.
       B_nl(:,:) = 0.
       grad2u(:,:) = 0.
       grad2v(:,:) = 0.
       grad4u(:,:) = 0.
       grad4v(:,:) = 0.
       thickness(:,:) = H1
       rhs_u(:,:,:) = 0.
       rhs_v(:,:,:) = 0.
       rhs_eta(:,:,:) = 0.
     ! Dummy arrays
       array(:,:) = 0.
       array_x(:,:) = 0.
       array_y(:,:) = 0.
       faces_array(:,:) = 0.
  
     ! Kronecker deltas
       top(:) = 0.
       bot(:) = 0.
       top(1) = 1.
       bot(nz) = 1.

       !count_specs_1 = 0
       !count_specs_2 = 0
       !count_specs_to = 0
       !count_specs_AG = 0
       !ke1_spec(:) = 0.
       !ke2_spec(:) = 0.
       !for_to_spec(:) = 0.
       !for_ag_spec(:) = 0.

       ! Baroclinic variables :
       Fmodes(:)    = 0.
       F_layer(:,:) = 0.
       A(:,:)       = 0. 


       ! Thicknesses parameters :

       WRITE (*,*) ''
       print *, '| --------> Diagno. from initialize.f90 <-------- |'

       print *, "Setting layer's thicknesses :"
       Hmin = 200
       exp_coef = 3*log(1.*7/2)
       exp_amp  = 1600/7
       Hsum = 0.
       H(:) = 0.
       
       do k = 1,nz 
          ratio = real(k)/(nz+1)
          print *, 'ratio',ratio
          print *, 'exp_coef', exp_coef
          H(k) = exp_amp*exp(exp_coef*ratio) + 200
          print *, exp_amp*exp(exp_coef*ratio) + 200
          Hsum = Hsum + H(k)

       end do
       do k = 1,nz
          H(k) = nint(H(k)*Htot/Hsum)
          WRITE (k_str,'(I0)') k
          print *, ' H(',TRIM(k_str),')   =  ',H(k)
       end do

       ! Stokes' thickness :
       !HS = Htot
       HS(:,:) = H(1)
       
       ! >>> Setting densities :
       print *, 'Setting densities :'
       ! Logarithmic stratification : 
       z = 0.
       do k = 1,nz
          z = z - H(k)/2
          rho(k) = 1028 - 2*exp(z/1000)
          WRITE (k_str,'(I0)') k
          print *, ' rho(',TRIM(k_str),')   =  ', rho(k), ' at z=', z
          z = z - H(k)/2
       end do
       
       ! >>> Printing Diagnostics :

       ! gprime
       g = 9.81
       write (*,*) 'Reduced gravities'
       gprime(1) = g
       print *, ' gprime(1) = ', gprime(1)
       if (nz.gt.1) then
       do k = 2,nz
          WRITE (k_str,'(I0)') k
          gprime(k) = g*(rho(k) - rho(k-1))/rho(1)
          print *, ' gprime(',TRIM(k_str),') = ', gprime(k)
       enddo
       end if
       
       ! --- Finding baroclinic streamfunc :
       ! > "Stolen" from Louis-Philippe "Naydo" (See David for reference to the name).
       ! > ***N.B. Here g(k) is defined for ceiling of layer k (NOT THE BOTTOM of layer k)
       ! > (Which is why it's also different than LP's code)
       if (nz.gt.1) then
       do k=1,nz-1
          F_layer(k,k)   = f0**2/( H(k)*gprime(k)  )
          F_layer(k,k+1) = f0**2/( H(k)*gprime(k+1))
       end do
       end if
       F_layer(nz,nz) = f0**2/( H(nz)*gprime(nz)  ) 

       ! > Creating matrix of streamfunc linear operator A (From left to right):
       ! > First column
       A(1,1) = (F_layer(1,2) ) !+ F_layer(1,1)) Rigid lid means we remove that
       A(2,1) = -  F_layer(2,2)
       do k=2,nz-1,1
          A(k-1,k) = -  F_layer(k-1,k)
          A(k,k)   =   (F_layer(k,k+1) + F_layer(k,k))
          A(k+1,k) = -  F_layer(k+1,k+1)
       end do
       ! > Last column
       A(nz-1,nz) = - F_layer(nz-1,nz) 
       A(nz,nz)   =   F_layer(nz,nz)

       ! Printing diagnostic :
       print *, 'A matrix (buyancy linear operator into matrix form) :'
       do k = 1,nz
          print *, ' [', A(k,:),']'
       end do

       ! > Solving Eigenvalues problem : A.psi = lambda.psi
       ! > N.B. Fmodes = eigenvalues of "mode" system of equation :
       !        q_mode = laplacian(psi_mode)+F*psi_mode
       CALL SGEEV('N','V', nz,A,nz, Fmodes,WI, VL,nz, L2M,nz, WORKL, size(workl,1), INFO )

       
       ! LAPACK eigenvalues
       PRINT *, 'LAPACK Eigenvalues of A matrix : '
       do k=1,nz
          WRITE (k_str,'(I0)') k
          IF (Fmodes(k) .lt. 1e-15) THEN
             Fmodes(k) = 0.
          ENDIF
          print *, ' lambda_',TRIM(k_str),' =', Fmodes(k)
       end do


       
       ! Printing eigenvectors (baroclinic modes)
       print *, "Eigenvectors (barotropic and baroclinic modes) :"
       do k=1,nz
          write (k_str,'(I0)') k
          write (*,*) " v(",trim(k_str),")  =  [",L2M(:,k),"]"
       end do

       ! Analytic equal 3-layers eigenvalues (test)
       ! print *, ' > Analytic Eigenvalues (3 equal layer, linear density): '
       ! print *, '   lambda_1 =', 1*f0**2/gprime(2)/H(1)
       ! print *, '   lambda_2 =', 3*f0**2/gprime(2)/H(1)
       ! print *, '   lambda_3 =', 0.                    
       
       ! Deformation radii
       print *, 'Deformation Radii :'
       do k=1,nz
          WRITE (k_str,  '(I0)'  ) k
          WRITE (ministr,'(F8.3)') 1/SQRT(Fmodes(k))/1000
          PRINT *, ' Ld(', TRIM(k_str), ') = 1/SQRT[lambda(',trim(k_str),')] = ', ministr, ' km'
       end do
       print *, 'Deformation Radii over dx:'
       do k=1,nz
          WRITE (k_str,  '(I0)'  ) k
          WRITE (ministr,'(F8.3)') 1/SQRT(Fmodes(k))/dx
          PRINT *, ' Ld(', TRIM(k_str), ') over dx = ', ministr
       end do
       print *, 'Baroclinic wave speed'
       do k=1,nz
          WRITE (k_str,  '(I0)'  ) k
          WRITE (ministr,'(F8.1)') f0/SQRT(Fmodes(k))
          PRINT *, ' c_bc(',TRIM(k_str),') = f0/Ld(', TRIM(k_str), ') = ', ministr
       end do
       
       print *, '| -------------------------------------------------- |'
       WRITE (*,*) ''

       do j = 0,nny
          f(j) = f0 + beta*(j-1)*dy
       end do
       
       ! Stohastic wind
       !do j = 1,ny
       !   y = -Ly/2. + (j-1)*dy
       !   taux_steady(:,j) = tau0*cos(twopi*y/Ly)
       !enddo
       !array = taux_steady
       !include 'subs/bndy.f90'
       !taux_steady = array
       !
       !taux_var(:,:) = tau1
       !array = taux_var
       !include 'subs/bndy.f90'
       !taux_var = array


!   --- Restart
       icount = 0 !for  output file index
       iftcount =0

       if ( restart .eqv. .false. ) then
          time = 0.  !in second
          restart_from=time
          print*,'Restart from',restart_from, 'day','icount,iftcount',icount,iftcount

       else !if restart
          open(0,file='restart')
          ! Two layers 
          !#!if (nz.eq.2) then
             do j = 1,ny
             do i = 1,nx             
                read(0,*)   u(i,j,:,1), &
                     &      v(i,j,:,1), &
                     &      eta(i,j,:,1),                         &
                     &      UStokes(i,j,1),VStokes(i,j,1),        &
                     &      taux_oc(i,j,1),tauy_oc(i,j,1)
             enddo
             enddo
          ! 5 layers
          !#!else if (nz.eq.5) then
          !#!   do j = 1,ny
          !#!   do i = 1,nx             
          !#!      read(0,*)   u(i,j,1,1),u(i,j,2,1),u(i,j,3,1),u(i,j,4,1),u(i,j,5,1),    &
          !#!           &      v(i,j,1,1),v(i,j,2,1),v(i,j,3,1),v(i,j,4,1),v(i,j,5,1),    &
          !#!           &      eta(i,j,2,1),eta(i,j,3,1),eta(i,j,4,1),eta(i,j,5,1),       &
          !#!           &      UStokes(i,j,1),VStokes(i,j,1),                             &
          !#!           &      taux_oc(i,j,1),tauy_oc(i,j,1)
          !#!   enddo
          !#!   enddo
          !#!else ! Number of layers doesn't fit
          !#!   print*, 'Cannot read restart file : need to implement for ',k,'layers.'
          !#!   stop
          !#!end if
          
          read(0,*) icount_srt,time,nspecfile,iftcount_srt
          close(0)
          restart_from=time/86400
          print*, 'Restart from', restart_from, 'day'
 
!         if (restart_from == 999 ) then
!            time = 0
!            icount = 0
!         else 
!            time = (restart_from-100)*iout* dt
!            icount = restart_from - 100
!         end if

!        print*, 'time =', time/86400. , 'days'

!         WRITE(which,'(I3)') restart_from

!         string1 = 'data/u'  // '_' // which(1:3)
!         string2 = 'data/v'  // '_' // which(1:3)
!         string3 = 'data/eta'  // '_' // which(1:3)
! 
!         open(unit=14,file=string1,access='DIRECT',&
!              & form='UNFORMATTED',status='UNKNOWN',RECL=4*(nx+2)*(ny+2)*2)
!         read(14,REC=1) u(:,:,:,1)
!         close(14)
!         
!         open(unit=15,file=string2,access='DIRECT',&
!              & form='UNFORMATTED',status='UNKNOWN',RECL=4*(nx+2)*(ny+2)*2)
!         read(15,REC=1) v(:,:,:,1)
!         close(15)
! 
!         open(unit=16,file=string3,access='DIRECT',&
!              & form='UNFORMATTED',status='UNKNOWN',RECL=4*(nx+2)*(ny+2)*2)
!         read(16,REC=1) eta(:,:,:,1)
!         close(16)

!         eta(:,:,1,:) = 0.

      endif  ! --- restart


!
!    Set bndy conditions
!
       do k = 1,nz

          array_x(:,:) = u(:,:,k,1)
          array_y(:,:) = v(:,:,k,1)
          include 'subs/no_normal_flow.f90'
          include 'subs/free_or_partial_slip.f90'
          u(:,:,k,1) = array_x(:,:)
          v(:,:,k,1) = array_y(:,:)

       enddo ! k = 1,nz
