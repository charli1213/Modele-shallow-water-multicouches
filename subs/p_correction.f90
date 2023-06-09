!
!     need to correct u,v for surface pressure 
!
       rhs_Psurf(:,:)= 0.
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
          include 'subs/divBT.f90' ! array = div(u*h) Barotropic Divergence
          rhs_Psurf(:,:) = rhs_Psurf(:,:) + array(:,:)
       enddo !end k-loop


       ! N.B. UStokes\VStokes = 0 (see initialize.f90) if coupling = .false.
       do j = 1, ny
       do i = 1, nx
          rhs_Psurf(i,j) = rhs_Psurf(i,j)                        &
               &         + (UStokes(i+1,j,2)-UStokes(i,j,2))/dx  &
               &         + (VStokes(i,j+1,2)-VStokes(i,j,2))/dy
          ! UEkman was here
       enddo
       enddo
       
       rhs_Psurf(:,:) = rhs_Psurf(:,:)/Htot !/dt
!
!      rhs_Psurf is vert average of d/dt (div u) 
!      -grad^2 p = -rhs_Psurf so that the corrected flow in non-divergent
!      N.B. the multiplication and division by dt is omitted
!      
!      solve nabla^2 P_surf = rhs_Psurf
!
       datr(:,:) = rhs_Psurf(1:nx,1:ny)
       include 'fftw_stuff/invLaplacian.f90'
       Psurf(1:nx,1:ny) = datr(:,:)
       array = Psurf
       include 'subs/bndy.f90'
       Psurf = array

       do j = 1,ny
       do i = 1,nx
          u(i,j,:,ilevel) = u(i,j,:,ilevel)  &
          &               - (Psurf(i,j)-Psurf(i-1,j))/dx ! *dt 
          v(i,j,:,ilevel) = v(i,j,:,ilevel)  &
          &               - (Psurf(i,j)-Psurf(i,j-1))/dy ! *dt
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

       ! We re-apply the p_correction,
       ! because we can have as much p-correction as we want (we just took 2 because it's enough)
       p_out(:,:) = p_out(:,:) + Psurf(:,:)
