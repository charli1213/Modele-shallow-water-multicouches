!
!     need to correct u,v with barotropic streamfunction found with MUDPACK. 
!
  
       ! Initialising qties.
       rhs_u_BT(:,:) = 0.
       rhs_v_BT(:,:) = 0.
       rhs_u_BC(:,:,:) = 0.
       rhs_v_BC(:,:,:) = 0.
       
       ! delta_psiBT(:,:)   = 0.
       ! Finding dh/dt foreach layer. First layer must fit with dh1 + dh2 + dh3 ... = 0, cause dH=0
       ! rhs_eta is the derivative of thickness and we find the one for h1, here.
       rhs_eta(:,:,1) = 0.
       do k = 2, nz
          rhs_eta(:,:,1) =  rhs_eta(:,:,1) - rhs_eta(:,:,k)
       enddo

       ! Barotropic loop : 
       do k = 1, nz
          uu(:,:) = u(:,:,k,ilevel)
          vv(:,:) = v(:,:,k,ilevel)
          
          if (k.eq.1) then
             thickness(:,:) =  H(k) - eta(:,:,k+1,ilevel) 
          elseif(k.eq.nz) then
             thickness(:,:) =  H(k) + eta(:,:,k,ilevel)
          else
             thickness(:,:) =  H(k) + eta(:,:,k,ilevel)  &
              &             -  eta(:,:,k+1,ilevel)
          endif

          ! Barotropic RHS of u and v.
          do j=1,ny-1
          do i=1,nx-1

          rhs_u_BT(i,j) = rhs_u_BT(i,j)                                             & 
          &             + rhs_u(i,j,k)*(thickness(i,j) + thickness(i-1,j))/Htot/2   &
          &             + uu(i,j)*(rhs_eta(i,j,k) + rhs_eta(i-1,j,k))/Htot/2
          rhs_v_BT(i,j) = rhs_v_BT(i,j)                                             &
          &             + rhs_v(i,j,k)*(thickness(i,j) + thickness(i,j-1))/Htot/2   &
          &             + vv(i,j)*(rhs_eta(i,j,k) + rhs_eta(i,j-1,k))/Htot/2

          enddo
          enddo

          ! Bndy
          array_x = rhs_u_BT
          array_y = rhs_v_BT
          include 'subs/no_normal_flow.f90'
          include 'subs/free_slip.f90'
          rhs_u_BT = array_x
          rhs_v_BT = array_y

          
       enddo !end k-loop


       
       ! baroclinic RHS_u,v
       ! note : no need for bndy conditions.
       do k = 1, nz
          rhs_u_BC(:,:,k) = rhs_u(:,:,k) - rhs_u_BT(:,:)
          rhs_v_BC(:,:,k) = rhs_v(:,:,k) - rhs_v_BT(:,:)
       enddo

       
       ! finding curl of rhs_u_BT OR RHS_zetaBT (same thing)
       ! note : no need for bndy conditions here : Boundaries are set to 0.
       do j = 2,ny-1
          jm = j-1
       do i = 2,nx-1
          im = i-1
          
          RHS_zetaBT(i,j,2) =  (rhs_v_BT(i,j)-rhs_v_BT(im,j))/dx    &
          &                 -  (rhs_u_BT(i,j)-rhs_u_BT(i,jm))/dy           
       enddo
       enddo

       ! See MUDPACK documentation for this. It makes solutions
       ! more trustworthy. 
       correction_RHSzetaBT(:,:) = RHS_zetaBT(:,:,2) - RHS_zetaBT(:,:,1)

       
    ! ######################################################## !
    !                                                          !
    !        RHS_zetaBT is RHS of the Poisson equation :       !
    !                                                          !
    !             nabla^2(d_psi_BT) = d_zeta_BT                !
    !                                                          !
    !      we solve for d_psi_BT instead of pressure gradient  !
    !                                                          !
    ! ######################################################## !

       call mud2(iparm,fparm,workm,coef,bndyc,correction_RHSzetaBT, & 
            &    correction_deltaPsiBT,mgopt,ierror)
       
    ! ######################################################## !
    !                                                          !
    !                -- DELTA_PSI_BT SOLVED --                 !
    !                                                          !   
    ! ######################################################## !

       delta_PsiBT(:,:,2) = correction_deltaPsiBT(:,:) + delta_PsiBT(:,:,1)
       
       
       ! Finding updated velocities (With TRUE barotropic RHS now)
       !
       !     rhs_u = rhs_u_BT + rhs_u_BC
       !
       ! Note : u = - \curl(\psi \kvec) = k \times \gradient(\psi)
       do j = 1,ny-1
       do i = 1,nx-1
          rhs_u_BT(i,j) =  - (delta_psiBT(i,j+1,2) - delta_psiBT(i,j,2))/dy  ! barotropic part-x
          rhs_v_BT(i,j) =    (delta_psiBT(i+1,j,2) - delta_psiBT(i,j,2))/dx  ! barotropic part-y
       enddo
       enddo

       ! Bndy
       array_x = rhs_u_BT
       array_y = rhs_v_BT
       include 'subs/no_normal_flow.f90'
       include 'subs/free_slip.f90'
       rhs_u_BT = array_x
       rhs_v_BT = array_y

       
       
       do k = 1,nz
          rhs_u(:,:,k) = rhs_u_BC(:,:,k) + rhs_u_BT(:,:)
          rhs_v(:,:,k) = rhs_v_BC(:,:,k) + rhs_v_BT(:,:)
       enddo
       ! Note : no need for boundary correction here.
       ! RHS_u,v are now updated! (Cheers!)

       
       ! Updating quantities for next timestep
       RHS_zetaBT(:,:,1)  = RHS_zetaBT(:,:,2)
       delta_PsiBT(:,:,1) = delta_PsiBT(:,:,2)