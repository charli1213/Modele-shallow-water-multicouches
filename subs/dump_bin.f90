
  eta(:,:,1,3) = Psurf(:,:)
  
  WRITE(which,'(I6)') 100000 + icount
     
! Note indices for (u,v,eta ...) starting with 0, useful part is 1:256
  !  real u_out(0:nx/subsmprto+1,0:ny/subsmprto+1,nz), v_out(0:nx/subsmprto+1,0:ny/subsmprto+1,nz)
  !  real div_ek_out(0:nx/subsmprto+1,0:ny/subsmprto+1), eta_out(0:nx/subsmprto+1,0:ny/subsmprto+1,nz)
  if (IO_field) then
    ! U and V field
    string1 = 'data/u_o'  // '_' // trim(which)
    string2 = 'data/v_o'  // '_' // trim(which)
    
    do k = 1,nz
      u_out(:,:,k)   = u(isubx,isuby,k,3)
      v_out(:,:,k)   = v(isubx,isuby,k,3)
      eta_out(:,:,k) = eta(isubx,isuby,k,3)
    enddo

    zeta1_out(:,:) = zeta1(isubx,isuby)
    div1_out(:,:)  = div1(isubx,isuby)
    
    open(unit=101,file=string1,access='DIRECT',&
          & form='BINARY',status='UNKNOWN',RECL=4*(size(isubx)*size(isubx)*nz))
    write(101,REC=1) (((u_out(i,j,k),i=1,nx/subsmprto),j=1,ny/subsmprto),k=1,nz)
    close(101)
    
    open(unit=102,file=string2,access='DIRECT',&
          & form='BINARY',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)*nz))
    write(102,REC=1) (((v_out(i,j,k),i=1,nx/subsmprto),j=1,ny/subsmprto),k=1,nz)
    close(102)

    ! ETA field
    !string3 = 'data/eta_o'  // '_' // trim(which)
    !open(unit=103,file=string3,access='DIRECT',&
    !& form='BINARY',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)*nz))
    !write(103,REC=1) (((eta_out(i,j,k),i=1,nx/subsmprto),j=1,ny/subsmprto),k=1,nz)
    !close(103)

    ! DIV_EK (W_EK)
    !string4 = 'data/div_o'  // '_' // trim(which)
    !div_ek_out=div_ek(isubx,isuby)
    !open(unit=104,file=string4,access='DIRECT',&
    !& form='BINARY',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)))
    !write(104,REC=1) ((div_ek_out(i,j),i=1,nx/subsmprto),j=1,ny/subsmprto)
    !close(104)


     ! Divergence AND curl of first layer [1/s]
     string39 = 'data/div1'  // '_' // trim(which)
     open(unit=139,file=string39,access='DIRECT',&
          & form='BINARY',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)))
     write(139,REC=1) ((div1_out(i,j),i=1,nx/subsmprto),j=1,ny/subsmprto)
     close(139)

     string40 = 'data/zeta1'  // '_' // trim(which)
     open(unit=140,file=string40,access='DIRECT',&
          & form='BINARY',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)))
     write(140,REC=1) ((zeta1_out(i,j),i=1,nx/subsmprto),j=1,ny/subsmprto)
     close(140)

     ! Divergence AND curl of second layer [1/s]
     !string41= 'data/div2'  // '_' // trim(which)
     !open(unit=141,file=string41,access='DIRECT',&
     !     & form='BINARY',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)))
     !write(141,REC=1) ((div2(i,j),i=1,nx/subsmprto),j=1,ny/subsmprto)
     !close(141)

     !string42= 'data/zeta2'  // '_' // trim(which)
     !open(unit=142,file=string42,access='DIRECT',&
     !     & form='BINARY',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)))
     !write(142,REC=1) ((zeta2(i,j),i=1,nx/subsmprto),j=1,ny/subsmprto)
     !close(142)

    
    
  end if !IO_field
  
  if (IO_coupling) then

     Ust_out(:,:) = Ust(isubx,isuby,2)
     Vst_out(:,:) = Vst(isubx,isuby,2)
     taux_ocean_out(:,:) = taux_ocean(isubx,isuby,2)
     tauy_ocean_out(:,:) = tauy_ocean(isubx,isuby,2)
     
     ! U Stokes
     string25 = 'data/Ust'  // '_' // trim(which)
     open(unit=125,file=string25,access='DIRECT',&
          & form='BINARY',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)))
     write(125,REC=1) ((Ust_out(i,j),i=1,nx/subsmprto),j=1,ny/subsmprto)
     close(125)

     string26 = 'data/Vst'  // '_' // trim(which)
     open(unit=126,file=string26,access='DIRECT',&
          & form='BINARY',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)))
     write(126,REC=1) ((Vst_out(i,j),i=1,nx/subsmprto),j=1,ny/subsmprto)
     close(126)

     !string27 = 'data/taux_eff'  // '_' // trim(which)
     !open(unit=127,file=string27,access='DIRECT',&
     !     & form='BINARY',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)))
     !write(127,REC=1) ((taux_eff(i,j),i=1,nx/subsmprto),j=1,ny/subsmprto)
     !close(127)

     !string28 = 'data/tauy_eff'  // '_' // trim(which)
     !open(unit=128,file=string28,access='DIRECT',&
     !     & form='BINARY',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)))
     !write(128,REC=1) ((tauy_eff(i,j),i=1,nx/subsmprto),j=1,ny/subsmprto)
     !close(128)

     string29 = 'data/taux_ocean'  // '_' // trim(which)
     open(unit=129,file=string29,access='DIRECT',&
          & form='BINARY',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)))
     write(129,REC=1) ((taux_ocean_out(i,j),i=1,nx/subsmprto),j=1,ny/subsmprto)
     close(129)
     
     string30 = 'data/tauy_ocean'  // '_' // trim(which)
     open(unit=130,file=string30,access='DIRECT',&
          & form='BINARY',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)))
     write(130,REC=1) ((tauy_ocean_out(i,j),i=1,nx/subsmprto),j=1,ny/subsmprto)
     close(130)


     
  endif ! IO_coupling
  
  if (IO_ek) then ! IO_ek
     div_ek_out(:,:)  = div_ek(isubx,isuby)
     zeta_ek_out(:,:) = zeta_ek(isubx,isuby)
     Uek_out(:,:)     = Uek(isubx,isuby,2)
     Vek_out(:,:)     = Vek(isubx,isuby,2)
     
     ! Divergence AND curl of Ekman layer (of Ekman transport) [m/s]
     string37 = 'data/div_ek'  // '_' // trim(which)
     open(unit=137,file=string37,access='DIRECT',&
          & form='BINARY',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)))
     write(137,REC=1) ((div_ek_out(i,j),i=1,nx/subsmprto),j=1,ny/subsmprto)
     close(137)

     string38 = 'data/zeta_ek'  // '_' // trim(which)
     open(unit=138,file=string38,access='DIRECT',&
          & form='BINARY',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)))
     write(138,REC=1) ((zeta_ek_out(i,j),i=1,nx/subsmprto),j=1,ny/subsmprto)
     close(138)

     
     ! Transport
     string5 = 'data/Uek'  // '_' // trim(which)
     open(unit=105,file=string5,access='DIRECT',&
          & form='BINARY',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)))
     write(105,REC=1) ((Uek_out(i,j),i=1,nx/subsmprto),j=1,ny/subsmprto)
     close(105)
     
     string6 = 'data/Vek'  // '_' // trim(which)
     open(unit=106,file=string6,access='DIRECT',&
          & form='BINARY',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)))
     write(106,REC=1) ((Vek_out(i,j),i=1,nx/subsmprto),j=1,ny/subsmprto)
     close(106)

     
  end if ! IO_ek


  if (IO_forcing) then
    ! Forcing-AG
      string7 = 'data/forci_ag'  // '_' // trim(which)
      open(unit=107,file=string7,access='DIRECT',&
      & form='BINARY',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)))
      write(107,REC=1) ((forcing_ag(i,j),i=1,nx/subsmprto),j=1,ny/subsmprto)
      close(107)

    ! Forcing
      string8 = 'data/forci_to'  // '_' // trim(which)
      open(unit=108,file=string8,access='DIRECT',&
      & form='BINARY',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)))
      write(108,REC=1) ((forcing_total(i,j),i=1,nx/subsmprto),j=1,ny/subsmprto)
      close(108)

  !  string9 = 'data/dissi_u'  // '_' // trim(which)
  !  string10 = 'data/dissi_v'  // '_' // trim(which)
  !  string11 = 'data/q'  // '_' // trim(which)

      string12 = 'data/taux'  // '_' // trim(which)
      open(unit=112,file=string12,access='DIRECT',&
      & form='BINARY',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)*nz))
      write(112,REC=1) ((taux(i,j),i=1,nx/subsmprto),j=1,ny/subsmprto)
      close(112)
  end if !IO_forcing
  
  if (IO_QGAG) then
      string13 = 'data/u_qg'  // '_' // trim(which)
  !  string14 = 'data/u_ag'  // '_' // trim(which)
      string15 = 'data/v_qg'  // '_' // trim(which)
  !  string16 = 'data/v_ag'  // '_' // trim(which)
  end if !IO_QGAG

  if(IO_psivort) then
    ! ETA-G
      string17 = 'data/eta_qg'  // '_' // trim(which)
      open(unit=117,file=string17,access='DIRECT',&
      & form='BINARY',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)*nz))
      write(117,REC=1) ((eta_qg(i,j),i=1,nx/subsmprto),j=1,ny/subsmprto)
      close(117)

  !  string18 = 'data/eta_ag'  // '_' // trim(which)

    ! ZETA-G
    string20 = 'data/zeta_G' // '_' // trim(which)
    open(unit=120,file=string20,access='DIRECT',&
    & form='BINARY',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)*nz))
    write(120,REC=1) (((zeta_G(i,j,k),i=1,nx/subsmprto),j=1,ny/subsmprto),k=1,nz)
    close(120)

    ! ZETA-AG
    string21 = 'data/zeta_AG' // '_' // trim(which)
    open(unit=121,file=string21,access='DIRECT',&
    & form='BINARY',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)*nz))
    write(121,REC=1) (((zeta_AG(i,j,k),i=1,nx/subsmprto),j=1,ny/subsmprto),k=1,nz)
    close(121)

    
    ! PSI
    string22 = 'data/PSImode' // '_' // trim(which)
    open(unit=122,file=string22,access='DIRECT',&
    & form='BINARY',status='UNKNOWN',RECL=4*(size(isubx)*size(isuby)*nz))
    write(122,REC=1) (((psimode(i,j,k),i=1,nx/subsmprto),j=1,ny/subsmprto),k=1,nz)
    close(122)
  end if !IO_psivort
