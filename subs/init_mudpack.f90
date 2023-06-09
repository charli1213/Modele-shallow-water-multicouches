   
  curl_of_RHS_u_BT(:,:) = 0.
  delta_psi_BT(:,:)  = 0.

  
  ! ************************************************** !
  !                                                    !
  !                    Init mudpack                    !
  !                                                    !
  ! ************************************************** !
  PRINT *, ">>> Initialising MUDPACK parameters/inputs"
  
  ! iparm : integer vector of length 17
  iparm(1)  = 0    ! intl : Initializing {0,1}={Yes,No}.
  iparm(2)  = 0    ! nxa  : Flag for boundary conditions.
  iparm(3)  = 0    ! nxb  ! {0}=periodic boundaries
  iparm(4)  = 0    ! nyc  ! {1}=Dirichlet boundary
  iparm(5)  = 0    ! nyd  ! {2}=mixed boundary condition (Neumann)

  iparm(6)  = ixp ! ixp
  iparm(7)  = jyq ! jyq Plus grand commun diviseur de nx et ny (512 ou 256)
  iparm(8)  = iex
  iparm(9)  = jey ! 2^[8] = nx = ny  >>>>  iparm(8 ou 9) = log2(nx ou ny)
  iparm(10) = nnx          ! nx : # de points dans l'intervalle [xa,xb] (incluant la frontière)
                          ! nx = ixp*(2**(iex-1)) + 1
  iparm(11) = nny          ! ny
  iparm(12) = 0           ! iguess : Prendre un guess (0=False)
  
  !     DOCUMENTATION :
  !
  !     Solution p(t) satisfies
  !
  !          l(p(t)) = r(t)
  !
  !     if the differential operator "l" has time dependence (either thru
  !     the coefficients in the pde or the coefficients in the derivative
  !     boundary conditions) then use p(t) as an initial guess to p(t+dt)
  !     when solving
  !
  !          l(t+dt)(p(t+dt)) = r(t+dt)
  !
  !     with intl = 0 for all time steps (the discretization must be repeated
  !     for each new "t"). either iguess = 0 (p(t) will then be an initial
  !     guess at the coarsest grid level where cycles will commence) or
  !     iguess = 1 (p(t) will then be an initial guess at the finest grid
  !     level where cycles will remain fixed) can be tried.
  !     > On va essayer les deux.
  
  iparm(13) = 4  ! maxcy  : the exact number of cycles executed between the finest and the coarsest
  iparm(14) = 0  ! method : Méthode conseillée si Cxx = Cyy partout. (Si ça chie, prendre 3)
  length = int(4*(nnx*nny*(10+0+0)+8*(nnx+nny+2))/3)
  iparm(15) = length ! Conseillé.

  
  write (*,100) (iparm(i),i=1,15)
  100 format(' > Integer input arguments ',/'      intl = ',I2,       &
           /'      nxa = ',I2,' nxb = ', I2,' nyc = ',I2, ' nyd = ',I2,  &
           /'      ixp = ',I2,' jyq = ',I2,' iex = ',I2,' jey = ',I2,    &
           /'      nnx = ',I3,' nny = ',I3,' iguess = ',I2,' maxcy = ',I2,  &
           /'      method = ',I2, ' work space estimate = ',I7)
  PRINT *, "     Calculated nx = ", iparm(6)*(2**(iparm(8)-1)) + 1  

  
  ! fparm : float point vector of length 6
  fparm(1) = 0.
  fparm(2) = Lx
  fparm(3) = 0.
  fparm(4) = Ly
  fparm(5) = 0.0 ! tolmax : Tolérance maximum.
                 ! Ils conseillent de mettre 0.0 et de passer à travers tous les cycles (maxcy)
  write(*,103) (fparm(i), i=1,5)
  103 format(/' > Floating point input parameters ',              &
           /'      xa = ',f9.1,' xb = ',f9.1,                     &
           /'      yc = ',f9.1,' yd = ',f9.1,                     &
           /'      tolerance (error control) =   ',e10.3)

  
  ! workm : one dimensionnal real save work space.
  ALLOCATE(workm(length))
  workm(:) = 0.0

  ! bndyc : Boundary conditions (Voir plus bas) :
  !
  !          a subroutine with  arguments (kbdy,xory,alfa,gbdy) which
  !          are used to input mixed boundary conditions to mud2. bndyc
  !          must be declared "external" in the program calling mud2.
  !          the boundaries are numbered one thru four and the mixed
  !          derivative boundary conditions are described below (see the
  !          sample driver code "tmud2.f" for an example of how bndyc is
  !          can beset up).
  !                                                                               
  !          * * * * * * * * * * * *  y=yd
  !          *       kbdy=4        *
  !          *                     *
  !          *                     *
  !          *                     *
  !          * kbdy=1       kbdy=2 *
  !          *                     *
  !          *                     *
  !          *                     *
  !          *       kbdy=3        *
  !          * * * * * * * * * * * *  y=yc
  !          x=xa                  x=xb

  ! coef : sous-routine des coefficient (Voir plus bas) :
  !
  !          cxx(x,y)*pxx + cyy(x,y)*pyy + cx(x,y)*px + cy(x,y)*py +
  ! 
  !          ce(x,y)*p(x,y) = r(x,y).
  
  ! RHS_mud : Vecteur nx par ny qui représente le RHS. 
  ! (Calculé plus haut)
  !


  ! phi : Vecteur nx par ny qui représente notre solution.
  !
  ! Si iparm(5) = nyc = 1, alors les conditions frontières doivent être incluses
  ! dans le vecteru de solution phi.
  ! Mais nous avons posé une solution Neumann qui s'applique sur la dérivée, donc ça
  ! peut être n'importe quoi.
  !
  ! Si iguess=0, alors phi doit quand même être initialisé à tous les points de grille.
  ! Ces valeurs vont être utilisées comme guess initial. Mettre tous à zéro si une
  ! solution approximative n'est pas illustrée.
  CALL RANDOM_NUMBER(delta_psi_BT)
  delta_psi_BT = delta_psi_BT/10
  ! Conditions Dirichlet
  !delta_psi_BT(1,1:nx) = phi(1,1:nx)
  !delta_psi_BT(nx,1:nx) = phi(nx,1:nx)
  !delta_psi_BT(1:nx,1) = phi(1:nx,1)
  !delta_psi_BT(1:nx,ny) = phi(1:nx,ny)


  ! mgopt
  !           an integer vector of length 4 which allows the user to select
  !           among various multigrid options. 
  mgopt(1) = 0 ! kcycle (Default)
  mgopt(2) = 2 ! iprer (Default)
  mgopt(3) = 1 ! ipost (Default)
  mgopt(4) = 3 ! intpol (Default)

  write (*,102) (mgopt(i),i=1,4)
  102 format(/' > Multigrid option arguments ', &
           /'      kcycle = ',i2, &
           /'      iprer = ',i2, &
           /'      ipost = ',i2, &
           /'      intpol = ',i2,/)

  
  
  ierror = 0 ! No error = 0 
  WRITE (*,*) " > Checking shapes :"
  WRITE (*,*) "     Shape iparm    =" ,SHAPE(iparm)
  WRITE (*,*) "     Shape fmarp    =" ,SHAPE(fparm)
  WRITE (*,*) "     Shape work     =" ,SHAPE(workm)
  WRITE (*,*) "     Shape rhs      =" ,SHAPE(curl_of_RHS_u_BT(1:nnx,1:nny))
  WRITE (*,*) "     Shape solution =" ,SHAPE(delta_psi_BT(1:nnx,1:nny))
  WRITE (*,*) "     Shape mgopt    =" ,SHAPE(mgopt)
  WRITE (*,*) " "

  ! initialising MUD2 function
  PRINT *, " > Initialising MUDPACK (iparm(1)=0)"
  call mud2(iparm,fparm,workm,coef,bndyc,curl_of_RHS_u_BT(1:nnx,1:nny),delta_psi_BT(1:nnx,1:nny),mgopt,ierror)
  PRINT *, "     ERROR =",ierror
  PRINT *, " "
  IF (ierror .gt. 0) THEN
     print *, "     Severe MUDPACK error, check documentation"
     stop
  ENDIF
  iparm(1) = 1 ! for subsequent calls
  
  ! ************************************************** !
  !                                                    !
  !                (END) Init MUDPACK                  !
  !                                                    !
  ! ************************************************** !
