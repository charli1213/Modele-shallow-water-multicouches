!
!      calc pv
!      called after time steps switched, so 2 and 3 are both at "n+1"
!      time level
!
     do k = 1,nz

        uu(:,:) = u(:,:,k,3)  
        vv(:,:) = v(:,:,k,3)

        !!! Here thickness is more of a 'relative' thickness.
        if (k.eq.1) then
           thickness(:,:) =  -eta(:,:,k+1,3)/H(1)
        elseif (k.eq.nz) then
           thickness(:,:) =  eta(:,:,k,3)/H(k)
        else
           thickness(:,:) =  (eta(:,:,k,3)-eta(:,:,k+1,3))/H(k)
        endif
        !!! À modifier (end)

        array = uu
        include 'subs/bndy.f90'
        uu = array

        array = vv
        include 'subs/bndy.f90'
        vv = array

        array = thickness
        include 'subs/bndy.f90'
        thickness = array

        do j = 1,ny
        do i = 1,nx
        !!! N.B. zeta is always calculated "on the spot" 
           zeta(i,j) =  (vv(i,j)-vv(i-1,j))/dx    &
                &    -  (uu(i,j)-uu(i,j-1))/dy 
        enddo
        enddo

        array = zeta
        include 'subs/bndy.f90'
        zeta = array

        ! q : Quasi-Geostrophic Potential vorticity (where beta = 0)
        ! This is the definition, works for any number of layers
        do j = 1, ny
        do i = 1, nx
           q(i,j,k) = zeta(i,j) -0.25*f(j)*(thickness(i,j)     +    &
                &                           thickness(i-1,j)   +    &
                &                           thickness(i,j-1)   +    &
                &                           thickness(i-1,j-1) )
        enddo
        enddo

        array(:,:) = q(:,:,k)
        include 'subs/bndy.f90'
        q(:,:,k) = array(:,:)
        
     enddo ! k loop (end)
     

     !  invert to get psimode     
     !do k=1,nz
     do k=1,2
        datr(:,:) = qmode(1:nx,1:ny,k)
        include 'fftw_stuff/invertQ.f90'
        psimode(1:nx,1:ny,k) = datr(:,:)

        array(:,:) = psimode(:,:,k)
        include 'subs/bndy.f90'
        psimode(:,:,k) = array(:,:)
     enddo

     psi(:,:,1) = psimode(:,:,1) - H(2)*psimode(:,:,2)/Htot
     psi(:,:,2) = psimode(:,:,1) + H(1)*psimode(:,:,2)/Htot

     array(:,:) = psi(:,:,1)
     include 'subs/bndy.f90'
     psi(:,:,1) = array(:,:)
     array(:,:) = psi(:,:,2)
     include 'subs/bndy.f90'
     psi(:,:,2) = array(:,:)

     !!! Begining of k-loop
     do k = 1,nz

        do j = 1,ny
        do i = 1,nx
           v_qg(i,j,k) =  (psi(i+1,j,1)-psi(i,j,1))/dx
           u_qg(i,j,k) = -(psi(i,j+1,1)-psi(i,j,1))/dy
        enddo
        enddo
        
        
        array(:,:) = v_qg(:,:,k)
        include 'subs/bndy.f90'
        v_qg(:,:,k) = array(:,:)

        array(:,:) = u_qg(:,:,k)
        include 'subs/bndy.f90'
        u_qg(:,:,k) = array(:,:)

        array(:,:) = v(:,:,k,2)
        include 'subs/bndy.f90'
        v(:,:,k,2) = array(:,:)

        array(:,:) = u(:,:,k,2)
        include 'subs/bndy.f90'
        u(:,:,k,2) = array(:,:)

        v_ag(:,:,k) = v(:,:,k,2)-v_qg(:,:,k)
        u_ag(:,:,k) = u(:,:,k,2)-u_qg(:,:,k)   

        
        array(:,:) = v_ag(:,:,k)
        include 'subs/bndy.f90'
        v_ag(:,:,k) = array(:,:)

        array(:,:) = u_ag(:,:,k)
        include 'subs/bndy.f90'
        u_ag(:,:,k) = array(:,:)
        end do ! end of k-loop.

        ! Finding eta-qg
        do j = 1,ny
        do i = 1,nx
           ! gprime(nk)
           ! psi(nx,ny,nz)
           ! psi
           eta_qg(i,j) = (f0/gprime(2))*0.25*(                        &
                &                 psi(i,j,2)-psi(i,j,1)               &
                &               + psi(i+1,j,2)-psi(i+1,j,1)           &
                &               + psi(i,j+1,2)-psi(i,j+1,1)           &
                &               + psi(i+1,j+1,2)-psi(i+1,j+1,1))
        enddo
        enddo
        array(:,:) = eta_qg(:,:)
        include 'subs/bndy.f90'
        eta_qg(:,:) = array(:,:)

        array(:,:) = eta(:,:,2,3)
        include 'subs/bndy.f90'
        eta(:,:,2,3) = array(:,:) 

        eta_ag(:,:) = eta(:,:,2,3)-eta_qg(:,:)

        array(:,:) = eta_ag(:,:)
        include 'subs/bndy.f90'
        eta_ag(:,:) = array(:,:)

        !do k=1,nz
        do k=1,2
           uu=u_qg(:,:,k)
           vv=v_qg(:,:,k)
           ! zeta_G
           do j = 1,ny
           do i = 1,nx
              zeta_G(i,j,k) =  (vv(i,j)-vv(i-1,j))/dx    &
                   &         -  (uu(i,j)-uu(i,j-1))/dy 
           enddo ! i-loop
           enddo ! j-koop
           array = zeta_G(:,:,k)
           include 'subs/bndy.f90'
           zeta_G(:,:,k)= array
        enddo ! k-loop
        
        !do k=1,nz
        do k=1,2
           uu=u_ag(:,:,k)
           vv=v_ag(:,:,k)
           ! zeta_AG
           do j = 1,ny
           do i = 1,nx
              zeta_AG(i,j,k) =  (vv(i,j)-vv(i-1,j))/dx    &
                   &         -  (uu(i,j)-uu(i,j-1))/dy 
           enddo ! i-loop
           enddo ! j-koop
           array = zeta_AG(:,:,k)
           include 'subs/bndy.f90'
           zeta_AG(:,:,k)= array
        enddo ! k-loop
        ! [***] À modifier (end)


        
!     ke1 = 0.
!     ke2 = 0.
!     do j = 1, ny
!     do i = 1, nx
!        ke1 = ke1 + (u(i,j,1,2)**2 + v(i,j,1,2)**2)/2.
!        ke2 = ke2 + (u(i,j,2,2)**2 + v(i,j,2,2)**2)/2.
!     enddo
!     enddo

!     ke1_qg = 0.
!     ke2_qg = 0.
!     do j = 1, ny
!     do i = 1, nx
!        ke1_qg = ke1_qg + ((psi(i+1,j,1)-psi(i,j,1))/dx)**2   &
!        &      + ((psi(i,j+1,1)-psi(i,j,1))/dy)**2 
!        ke2_qg = ke2_qg + ((psi(i+1,j,2)-psi(i,j,2))/dx)**2   &
!        &      + ((psi(i,j+1,2)-psi(i,j,2))/dy)**2 
!     enddo
!     enddo
!     ke1_qg = ke1_qg/2.
!     ke2_qg = ke2_qg/2.

!     pe_qg = 0.
!     pe = 0.
!     do j = 1,ny
!     do i = 1,nx
!        pe_qg = pe_qg + (f0/c_bc)**2 * (psi(i,j,2)-psi(i,j,1))**2 
!        pe = pe + gprime(2)*eta(i,j,2,2)**2
!     enddo
!     enddo
!     pe_qg = pe_qg/2.
!     pe = pe/2.
!     pe = pe*Htot/H(1)/H(2)
!     etot = (H(1)*ke1 + H(2)*ke2)/Htot + pe
!     etot_qg = (H(1)*ke1_qg + H(1)*ke2_qg)/Htot + pe_qg

!      write(300,*), time/86400., ke1/nx/ny, ke2/nx/ny, pe/nx/ny   &
!          &                   , etot/nx/ny
!      write(301,*), time/86400., ke1_qg/nx/ny, ke2_qg/nx/ny   &
!          &                   , pe_qg/nx/ny, etot_qg/nx/ny 
!      write(302,*), time/86400., eta(1,ny/2,2,1), eta(1,ny/4,2,1)   &
!          &                   , eta(1,1,2,1)
!      call flush(300)
!      call flush(301)
!      call flush(302)
!
!
