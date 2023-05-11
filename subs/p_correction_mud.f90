!
!     need to correct u,v for surface pressure 
!
       ! Restarting qties.
       zeta_BT(:,:) = 0.
       u_BT(:,:)    = 0.
       v_BT(:,:)    = 0.
       
       ! Calculating thickness
       do k = 1, nz
          if (k.eq.1) then
             thickness(:,:) =  H(k) - eta(:,:,k+1,ilevel) 
          elseif(k.eq.nz) then
             thickness(:,:) =  H(k) + eta(:,:,k,ilevel)
          else
             thickness(:,:) =  H(k) + eta(:,:,k,ilevel)  &
              &             -  eta(:,:,k+1,ilevel)
          endif

          ! Finding barotropic zeta (zeta_BT)
          ! N.B. boundaries are included in zetaBT.f90
          uu(:,:) = u(:,:,k,ilevel)
          vv(:,:) = v(:,:,k,ilevel)
          include 'subs/zetaBT.f90' ! array = curl(u*h) or barotropic vorticity.
          zeta_BT(:,:) = zeta_BT(:,:) + array(:,:)
          u_BT(:,:)    = u_BT(:,:)    + uh(:,:)
          v_BT(:,:)    = v_BT(:,:)    + vh(:,:)
          ! (***) Don't we also put Stokes transport (Ust) here too? 
          
       enddo !end k-loop
       
       ! barotropic qty : 
       zeta_BT = zeta_BT/Htot
       u_BT    = u_BT/Htot
       v_BT    = v_BT/Htot
       
!
!      zeta_BT is RHS of Poisson equation 
!
!          nabla^2 psi_BT = zeta_BT
!
!      we solve for psi_BT
!
       include 'subs/mudpack_solver.f90'
       array = psi_BT
       include 'subs/bndy.f90'
       psi_BT = array
 
       ! Finding updated velocities
       !
       !     u = u_BT + u_BC = psi_y + (\tilde{u} - \tilde{u}_BT)
       !
       do j = 1,ny
       do i = 1,nx
          u(i,j,:,ilevel) = u(i,j,:,ilevel)                    &     ! \tilde{u(k)}
               &          + (psi_BT(i,j+1) - psi_BT(i,j))/dy   &     ! mudpack psi_y
               &          - u_BT(i,j)                                ! \tilde{u}_BT
          v(i,j,:,ilevel) = v(i,j,:,ilevel)                    &     ! \tilde{v(k)}
               &          - (psi_BT(i+1,j) - psi_BT(i,j))/dx   &     ! mudpack -psi_x
               &          - v_BT(i,j)                                ! \tilde{v}_BT
       enddo
       enddo

       do k = 1,nz
       array(:,:) = u(:,:,k,ilevel)
       include 'subs/bndy.f90'
       u(:,:,k,ilevel) = array(:,:)
       array(:,:) = v(:,:,k,ilevel)
       include 'subs/bndy.f90'
       v(:,:,k,ilevel) = array(:,:)
       enddo ! end k-loop

       ! Velocities are now updated! (Cheers)