 SUBROUTINE sc_solsh( sec )
 !            --------------------------------
 !            input material properties for multi-layered lamina
 !            Total Lagrangean Solid-shell elements (SOLSH)
 !  iprop: 1 :nlayr
 !         2 :nvars:  number of internal variables for the whole laminate (NOT YET)
 !         3 :0  elastic lamina,  1:inelastic   (0 ELASTIC ONLY)
 !         4 :
 !    5:4+nl :associated material label to each layer
 !
 !  rprop  1 : thick¿?
 !         2 : mean young modulus
 !         3 : mean density
 !         4 : shear factor (5/6 is the standard value)
 !         5 : surface density
 !      6-44 : (39) Standard integrated elasticity matrix
 !     45-86 : (42) in-plane integrated matrices associated with zig-zag function
 !     87-93 : ( 7) transverse integrated matrices associated with zig-zag function
 ! shifted 5 it is
 !             6  7  8  9 29 30 31 31        55 56 57 58
 !               10 11 12 33 34 35 36        59 60 61 62
 !                  13 14 37 38 39 40        63 64 65 66
 !                     15 41 42 43 44        67 68 69 70
 !                        16 17 18 19        71 72 73 74
 !                           20 21 22        75 76 77 78
 !                              23 24        79 80 81 82
 !                                 25        83 84 85 86
 !                                    26 27              90 91
 !                                       28              92 93
 !                                           45 46 47 48
 !                                              49 50 51
 !                                                 52 53
 !                                                    54
 !                                                       87 88
 !                                                          89
 !      for assumed shear strain
 !            94 95 100 101 107 108
 !               96 102 103 109 110
 !                   97  98 111 112
 !                       99 113 114
 !                          104 105
 !                              106
 !
 !
 !
 !    2 nl   : phi (zig-zag) functions
 !    2 nl   : Beta (zig-zag derivative) values
 !      For each layer (npt values)
 !          1: layer thickness
 !          2: layer angle
 !          3: layer z coordinate
 !        4-8: rotation coefficients
 !          9: shear factor

 IMPLICIT NONE
 TYPE (section), POINTER :: sec

 !local variables
 INTEGER (kind=4), PARAMETER :: npl = 9     !number of real values per layer
 INTEGER (kind=4), PARAMETER :: matsa(2) = (/1,5/)  ! material types allowed for this section
 REAL (kind=8), PARAMETER :: factor = 1.7453292519943295769236907684886d-2   !pi/180

 INTEGER (kind=4) matno,nlayr,ith,il,sizer,sizei,ii,imat,k
 LOGICAL :: found,read_integrated_matrix,fastr
 REAL(kind=8) :: thick,youn,dens,thickl,de(109),angle,zl,ce(13),sg,cg,celoc(9),gi(2),g1, &
                 f11,f22,f12,f1,f2,sf,im(0:2),thicks,ds(10)
 REAL(kind=8), ALLOCATABLE :: rprop(:,:),phi(:,:),beta(:,:)
 TYPE (mater), POINTER :: mat
 TYPE (postv), POINTER :: var1,var2,var3,var4,var5 !pointers to postprocess variables

 sec%secty =  5         ! Laminated solid-shell
 sec%mabas = -1         ! flag to indicate that no unique associated material exist

 CALL listen('SOLSH')  !read a line

 sizer = 117            ! minimum size of real array
 thick = getrea('THICK ',0.D0,' Thickness ........................')
 IF( thick > 0d0 ) THEN
   ith=1 ! thickl will be input as a percentage of thick
   thicks = 0d0 !initializes thickness to check
 ELSE
   ith=0 ! thickl will be input as a Real dimension
   thick=0d0   !Initializes thickness
 END IF
 nlayr  = getint('NLAYR ',   0,'!NUMBER OF LAYERS .................')
 IF ( nlayr <= 0) CALL runend('SOLSH: VERIFY NUMBER OF LAYERS ')
 sizei =  4+nlayr*2                !size of integer array
 ALLOCATE (sec%iprop(sizei), rprop(npl,nlayr) )
 sec%iprop = 0  !initializes
 rprop = 0d0    !initializes
 sf = getrea('SHEARF',0.83333D0,' Shear Factor .....................') !for FSDT
 read_integrated_matrix = exists('READMA')

 youn = 0d0            !Initializes mean young modulus
 dens = 0d0            !Initializes mean density
 ii  = 5               !pointer to integer array
 gi = 0d0              !average values of trasnverse modulus
 im = 0d0              !rot inertia

 DO il=1,nlayr         !for each layer
   ! read layer information
   CALL listen('SOLSH')  !read a line
   matno    = getint('NUMAT ', 1 ,'!Associated Material ..............')
   CALL mat_search(matno,found,mat,imat)        !search base material
   IF(.NOT.found) CALL runend ('SOLSH: NOT FOUND ASSOCIATED MATERIAL ')
   IF( .NOT.ANY( mat%mtype == matsa) ) &  !check if base material is correct
            CALL runend('SOLSH: Invalid material association')
   sec%iprop(ii)       = imat     !store associated material index (in PMATS)
   IF(ith == 0) THEN   !thicknes must be computed
     thickl = getrea('THICKL',1.D0,'!LAYER THICKNESS ..................')
     thick  = thick + thickl         ! sum layer thickness
   ELSE                !thickness is known
     thickl = getrea('THICKL',1.D+2/nlayr,' LAYER THICKNESS ( % ).............')
     thickl = thick * thickl *1.D-2   ! compute layer thickness
     thicks = thicks + thickl         !
   END IF
   angle = getrea('ANGLE ',0.D0, ' ANGLE FROM COMP COORDS TO LA MAT 1')
   youn = youn + thickl*mat%prope(1)  !mean young modulus
   dens = dens + thickl*mat%prope(5)  !surface density
   rprop(1,il) = thickl               !store layer thickness (real dimension)
   rprop(2,il) = angle                !store layer angle     (degrees)
   angle = angle*factor       !angle between orthotropic axes and global axes in rad
   sg = SIN(angle)            !sine
   cg = COS(angle)            !cosine
   ! matrix to transform global strains into orthotropic strains
   rprop(4,il) =  cg**2                ! [ c^2     s^2     c s       ]
   rprop(5,il) =  sg**2                ! [ s^2     c^2    -c s       ]
   rprop(6,il) =  cg*sg                ! [-2 c s   2 c s   c^2 - s^2 ]
   rprop(7,il) =  2d0*cg*sg
   rprop(8,il) =  cg**2 - sg**2
   !compute average values
   SELECT CASE (mat%mtype)
   CASE(1)   !isotropic
     gi(:) = gi(:)+thickl/mat%prope(3)    ! Int 1/G
   CASE(5)   !orthotropic
     ! use angle to compute local cartesyan values
     g1 = mat%prope(14) * rprop(4,il) + mat%prope(15) * rprop(5,il)  !G13 c^2 + G23 s^2
     gi(1) = gi(1)+thickl/g1
     g1 = mat%prope(14) * rprop(5,il) + mat%prope(15) * rprop(4,il)  !G13 s^2 + G23 c^2
     gi(2) = gi(2)+thickl/g1
   END SELECT
   ii = ii+1                           !updates pointer to integer array
 END DO
 gi(1) = thick/gi(1)
 gi(2) = thick/gi(2)

 IF( ith == 1 .AND. ABS(thicks/thick - 1d0) > 0.001 )THEN
   WRITE(lures,"(' layers thickness sums',f8.4,'%')")thicks/thick*1d2
   CALL runend('SOLSH: Sum of layers /= 100%')
 END IF
 ALLOCATE (phi(2,0:nlayr),beta(2,nlayr))
 phi(:,0) = 0d0
 !      For each layer (npt values)
 !          1: layer thickness
 !          2: layer angle
 !          3: layer z coordinate
 !        4-8: rotation coefficients
 !          9: shear factor             !what for?
 sizer = sizer+(nlayr+1)*2              ! include zig-zag functions
 sizer = sizer+nlayr*2                  ! include zig-zag functions derivatives
 sec%iprop(1) = nlayr                   !keep number of layers
 sec%iprop(2) = 0                       !set total number of internal variables
 sec%iprop(3) = 0                       !initializes to elastic lamina

 sizer = sizer+nlayr*npl                !update size of real array including RPROP
 sec%secdef = (/ sizei,sizer, 0,0,0 /)  !size of arrays
 ALLOCATE (sec%rprop(sizer) )           !get memory
 sec%rprop = 0d0                        !initializes
 sec%rprop(1) = thick                   !thickness
 sec%rprop(2) = youn/thick              !mean young modulus
 sec%rprop(3) = dens/thick              !mean density
 sec%rprop(4) = sf                      !keep shear factor
 sec%rprop(5) = dens                    !surface density

 IF(.NOT.read_integrated_matrix)THEN
   ! compute elastic integrated matrix
   de = 0d0                             !initializes integrated thickness
   ds = 0d0                             !initializes auxiliar values for transverse shear
   zl = -thick/2d0      !initializes Z position at bottom surface
   DO il=1,nlayr                           !for each layer
     imat= sec%iprop(4+il)                   !material index
     mat => pmats(imat)%p                    !point to material
     ! celoc(1:6) in plane values  celoc(7:9) shear values
     CALL elastiff(celoc,mat,5)              !get elastic stiffnes matrix 5=3D Solid Principal directions
     thickl = rprop(1,il)                        !layer thickness
     g1 = thickl*thickl/12d0                     !auxiliar h^2/12
     zl = zl+thickl/2d0                          !update Z position to layer mid-point z_k-bar
     rprop(3,il) = zl                            !keep
     CALL rotor3d(celoc,rprop(4:,il),ce,.TRUE.)  !rotate elastic in-plane matrix  11 22 12 33 13 23
     beta(1,il) = gi(1)/ce(11) - 1d0
     beta(2,il) = gi(2)/ce(13) - 1d0
     phi(1,il) = phi(1,il-1) + beta(1,il)*thickl
     phi(2,il) = phi(2,il-1) + beta(2,il)*thickl
     f1 = (phi(1,il) + phi(1,il-1))/2d0          !phi_k-bar 1
     f2 = (phi(2,il) + phi(2,il-1))/2d0          !phi_k-bar 2
     ! ¿¿¿
     im(0) = im(0) + thickl*mat%prope(5)*(zl**2+thickl**2/12d0)      !rotatory inertia
     im(1) = im(1) + thickl*mat%prope(5)*(f1**2+(beta(1,il)*thickl)**2/12d0)   !additional rotatory inertia
     im(2) = im(2) + thickl*mat%prope(5)*(f2**2+(beta(2,il)*thickl)**2/12d0)
     ! ???
     ce = ce*thickl                              !scale layer stiffness matrix
!
!         | 1  2  3  4       |
!         | 2  5  6  7       |
!      Ce=| 3  6  8  9       |
!         | 4  7  9 10       |
!         |            11 12 |
!         |               13 |
!
 !  components of the integrated stiffness matrix
 !             1  2  3  4 24 25 26 27        50 51 52 53            Exx
 !             2  5  6  7 28 29 30 31        54 55 56 57            Eyy
 !             3  6  8  9 32 33 34 35        58 59 60 61            Exy
 !             4  7  9 10 36 37 38 39        62 63 64 65            Ezz
 !            24 28 32 36 11 12 13 14        66 67 68 69            Xxx
 !            25 29 33 37 12 15 16 17        70 71 72 73            Xyy
 !            26 30 34 38 13 16 18 19        74 75 76 77            Xxy
 !            27 31 35 39 14 17 19 20        78 79 80 81            Alpha
 !                                    21 22              85 86      Gx
 !                                    22 23              87 88      Gy
 !            50 54 58 62 66 70 74 78        40 41 42 43            Psi_x,x
 !            51 55 59 63 67 71 75 79        41 44 45 46            Psi_y,y
 !            52 56 60 64 68 72 76 80        42 45 47 48            Psi_x,y
 !            53 57 61 65 69 73 77 81        43 46 48 49            Psi_y,x
 !                                                       82 83      Psi_x
 !                                                       83 84      Psi_y
 !            89  90  95  96 102 103  C13(1)
 !            90  91  97  98 104 105  C23(1)
 !            95  97  92  93 106 107  C13(2)
 !            96  98  93  94 108 109  C23(2)
 !           102 104 106 108  99 100  Psi_x
 !           103 105 107 109 100 101  Psi_y

     !standard part
     de( 1:10) = de( 1:10) + ce(1:10)              !membrane part
     de(11:20) = de(11:20) + ce(1:10)*(zl*zl+g1)   !bending part
     de(21:23) = de(21:23) + ce(11:13)             !shear part
     !coupled part
     de(24:39) = de(24:39) + ce((/1:4,2,5:7,3,6,8,9,4,7,9,10/))*zl
     !de(24:27) = de(24:27) + ce(1:4)*zl              !E11 - Curv
     !de(28:31) = de(28:31) + ce((/2,5,6,7/))*zl      !E22 - Curv
     !de(32:35) = de(32:35) + ce((/3,6,8,9/))*zl      !E12 - Curv
     !de(36:39) = de(36:39) + ce((/4,7,9,10/))*zl     !E33 - Curv
     ! Dp(3,3)  phi^T Dk phi                !40 41 42 43
     f11 = f1*f1+beta(1,il)**2*g1           !41 44 45 46
     f12 = f1*f2+beta(1,il)*beta(2,il)*g1   !42 45 47 48
     f22 = f2*f2+beta(2,il)**2*g1           !43 46 48 49
     de(40) = de(40) + ce(1)*f11
     de(41) = de(41) + ce(2)*f12
     de(42) = de(42) + ce(3)*f11
     de(43) = de(43) + ce(3)*f12
     de(44) = de(44) + ce(5)*f22
     de(45) = de(45) + ce(6)*f12
     de(46) = de(46) + ce(6)*f22
     de(47) = de(47) + ce(8)*f11
     de(48) = de(48) + ce(8)*f12
     de(49) = de(49) + ce(8)*f22
     ! Dp(1,3)                                     !50 51 52 53
     de(50) = de(50) + ce(1)*f1                    !54 55 56 57
     de(51) = de(51) + ce(2)*f2                    !58 59 60 61
     de(52) = de(52) + ce(3)*f1                    !62 63 64 65
     de(53) = de(53) + ce(3)*f2
     de(54) = de(54) + ce(2)*f1
     de(55) = de(55) + ce(5)*f2
     de(56) = de(56) + ce(6)*f1
     de(57) = de(57) + ce(6)*f2
     de(58) = de(58) + ce(3)*f1
     de(59) = de(59) + ce(6)*f2
     de(60) = de(60) + ce(8)*f1
     de(61) = de(61) + ce(8)*f2
     de(62) = de(62) + ce(4)*f1
     de(63) = de(63) + ce(7)*f2
     de(64) = de(64) + ce(9)*f1
     de(65) = de(65) + ce(9)*f2
     !Dp(2,3)                                      !66 67 68 69
     f11 = f1*zl+beta(1,il)*g1                     !70 71 72 73
     f22 = f2*zl+beta(2,il)*g1                     !74 75 76 77
     de(66) = de(66) + ce(1)*f11                   !78 79 80 81
     de(67) = de(67) + ce(2)*f22
     de(68) = de(68) + ce(3)*f11
     de(69) = de(69) + ce(3)*f22
     de(70) = de(70) + ce(2)*f11
     de(71) = de(71) + ce(5)*f22
     de(72) = de(72) + ce(6)*f11
     de(73) = de(73) + ce(6)*f22
     de(74) = de(74) + ce(3)*f11
     de(75) = de(75) + ce(6)*f22
     de(76) = de(76) + ce(8)*f11
     de(77) = de(77) + ce(8)*f22
     de(78) = de(78) + ce(4)*f11
     de(79) = de(79) + ce(7)*f22
     de(80) = de(80) + ce(9)*f11
     de(81) = de(81) + ce(9)*f22
     ! Dt(2,2)                                     !82 83
     f1 = beta(1,il)                               !   84
     f2 = beta(2,il)
     f11 = f1*f1
     f12 = f1*f2
     f22 = f2*f2
     de(82) = de(82) + ce(11)*f11
     de(83) = de(83) + ce(12)*f12
     de(84) = de(84) + ce(13)*f22
     ! Dt(1,2)                                     !85 86
     de(85) = de(85) + ce(11)*f1                   !87 88
     de(86) = de(86) + ce(12)*f2
     de(87) = de(87) + ce(12)*f1
     de(88) = de(88) + ce(13)*f2
     ! AUXILIAR ARRAYS FOR ANS TRANSVERSE SHEAR
     ! for assumed strain transverse shear D_tm(1)
     f1 = zl*2d0/thick
     ds(1:3) = ds(1:3) + ce(11:13)*f1
     ! for assumed strain transverse shear D_tm(2)
     f1 = (thickl**2/3d0+4d0*zl**2)/thick**2
     ds(4:6) = ds(4:6) + ce(11:13)*f1
     ! for assumed strain transverse shear D_tm(3)
     f1 = zl*beta(1,il)*2d0/thickl
     f2 = zl*beta(2,il)*2d0/thickl
     ds( 7) = ds( 7) + ce(11)*f1
     ds( 8) = ds( 8) + ce(12)*f2
     ds( 9) = ds( 9) + ce(12)*f1
     ds(10) = ds(10) + ce(13)*f2
     !inertia values
     zl = zl+thickl/2d0                          !update Z position to layer interface
   END DO
   de(89:91) = (de(21:23) - 2d0*ds(1:3) + ds(4:6))/4d0          !Dt_11
   de(92:94) = (de(21:23) + 2d0*ds(1:3) + ds(4:6))/4d0          !Dt_22
   de(95:98) = (de((/21:22,22:23/)) - ds((/4:5,5:6/)))/4d0      !Dt_12
   de(99:101) = de(82:84)                                       !Dt_33
   de(102:105) = (de(85:88)-ds(7:10))/2d0                       !Dt_13
   de(106:109) = (de(85:88)+ds(7:10))/2d0                       !Dt_23
   WRITE(lures,"(' phi and Beta Values')")
   WRITE(lures,"(2e15.8)")phi(:,0)
   DO il=1,nlayr
     WRITE(lures,"(5x,2e15.8)")beta(:,il)
     WRITE(lures,"(2e15.8)")phi(:,il)
   END DO
   WRITE(lures,"(' Integrated D matrices')")
   WRITE(lures,"(12e15.7)")de((/  1, 2, 3, 4,24,25,26,27,       50,51,52,53/))
   WRITE(lures,"(12e15.7)")de((/  2, 5, 6, 7,28,29,30,31,       54,55,56,57/))
   WRITE(lures,"(12e15.7)")de((/  3, 6, 8, 9,32,33,34,35,       58,59,60,61/))
   WRITE(lures,"(12e15.7)")de((/  4, 7, 9,10,36,37,38,39,       62,63,64,65/))
   WRITE(lures,"(12e15.7)")de((/ 24,28,32,36,11,12,13,14,       66,67,68,69/))
   WRITE(lures,"(12e15.7)")de((/ 25,29,33,37,12,15,16,17,       70,71,72,73/))
   WRITE(lures,"(12e15.7)")de((/ 26,30,34,38,13,16,18,19,       74,75,76,77/))
   WRITE(lures,"(12e15.7)")de((/ 27,31,35,39,14,17,19,20,       78,79,80,81/))
   WRITE(lures,"(12e15.7)")de((/                         21,22,             85,86 /))
   WRITE(lures,"(12e15.7)")de((/                         22,23,             87,88 /))
   WRITE(lures,"(12e15.7)")de((/ 50,54,58,62,66,70,74,78,       40,41,42,43/))
   WRITE(lures,"(12e15.7)")de((/ 51,55,59,63,67,71,75,79,       41,44,45,46/))
   WRITE(lures,"(12e15.7)")de((/ 52,56,60,64,68,72,76,80,       42,45,47,48/))
   WRITE(lures,"(12e15.7)")de((/ 53,57,61,65,69,73,77,81,       43,46,48,49/))
   WRITE(lures,"(12e15.7)")de((/                                            82,83/))
   WRITE(lures,"(12e15.7)")de((/                                            83,84/))
   WRITE(lures,"(12e15.7)")de((/ 89, 90, 95, 96,102,103/))
   WRITE(lures,"(12e15.7)")de((/ 90, 91, 97, 98,104,105/))
   WRITE(lures,"(12e15.7)")de((/ 95, 97, 92, 93,106,107/))
   WRITE(lures,"(12e15.7)")de((/ 96, 98, 93, 94,108,109/))
   WRITE(lures,"(12e15.7)")de((/102,104,106,108, 99,100/))
   WRITE(lures,"(12e15.7)")de((/103,105,107,109,100,101/))
   WRITE(lures,"(' Integrated Inertias')")
   WRITE(lures,"(5e15.8)")im
 ELSE
   fastr = .FALSE.
   CALL rdfrre('SC_SOLSH',de,il,109,fastr)
   CALL rdfrre('SC_SOLSH',im,il,3,fastr)
 END IF
 sec%rprop(6:114)   = de                         !keep integrated matrix
 sec%rprop(115:117) = im                         !keep integrated inertias
 CALL vecasi (2*nlayr+2,phi,sec%rprop(118))      !keep phi
 k = 120+2*nlayr
 CALL vecasi (2*nlayr,beta,sec%rprop(k))         !keep beta
 k = k+2*nlayr
 CALL vecasi (2*nlayr,beta,sec%rprop(k))         !keep rprop                           sec%rprop(118:sizer) = rprop                         !keep integrated inertias
 !IF( sec%iprop(4) > 0 ) THEN  !If shear stress is considered
 !  !- Compute shear factors for each layer --------------------------------------
 !  DO il = 1, nlayr                         ! Shear stress through the section - XMG
 !    sec%rprop(k) = - ( 3.0d0 / (2.0d0 * sec%rprop(1)) )                & ! Computed as: 1/2I * (h^2/4 - z^2) * V
 !                   + ( 6.0d0 * (rprop(3,il)**2) / (sec%rprop(1)**3) )    ! Where I = b * h^3 / 2 (with b=1)
 !    k = k+1                                ! h = section thickness
 !  END DO                                   ! z = layer midpoint
 !END IF
 DEALLOCATE (rprop,phi,beta)                 !release local memory
 !sec%secdef(3) =  5        !GENERAL orthotropic with COUPLING


 ! THIS COULD BE USED TO CONSIDER SHELL TYPE ELEMENTS IN A SOLID ELEMENT
 ! MAXIMUM Additional postprocess variables --------------------------------------------
 sec%secdef(4) = 5  ! variables
 sec%secdef(5) =14  ! values

 ALLOCATE(var1,var2,var3,var4,var5)  ! new variables

 var1%type = 2
 var1%dim = 3
 var1%name(:) = (/ 'Forces','N_xx','N_yy','N_xy',' ',' ',' ' /)

 var2%type = 2
 var2%dim = 3
 var2%name(:) = (/ 'Moments','M_xx','M_yy','M_xy',' ',' ',' ' /)

 var3%type = 1
 var3%dim = 2
 var3%name(:) = (/ 'Shear','Shear_X','Shear_y',' ',' ',' ',' ' /)

 var4%type = 2
 var4%dim = 4
 var4%name(:) = (/ 'Phi_Moment','Phi_Mxx','Phi_Myy','Phi_Mxy','Phi_Myx',' ',' ' /)

 var5%type = 1
 var5%dim = 2
 var5%name(:) = (/ 'Phi_Shear','Phi_Shear_X','Phi_Shear_y',' ',' ',' ',' ' /)

 sec%postp => var1   !keep pointer to first variable
 var1%next => var2
 var2%next => var3
 var3%next => var4
 var4%next => var5
 NULLIFY (var5%next)

 !------------------------------------------------------------------------------
 RETURN
 END SUBROUTINE sc_solsh
!********************************************************************
 SUBROUTINE rotor3d (c,rr,ce,flag)
 IMPLICIT NONE
 !
 !     rotate the stress-strain matrix to global/reference coordinates
 !
 !     set the coordinate transformation for rotation of properties
 !     INPUT
 !        | c(1)  c(2) c(3)       | E1
 !   Cm = | c(2)  c(4) c(5)       | E2    Cs = | c(8)      ! G13
 !        | c(3)  c(5) c(6)       | E3         |      c(9) ! G23
 !        |                  c(7) |2E12
 !  r = { cos^2  sin^2  cos*sin  2cos*sin  con^2-sin^2 }
 !
 !     OUTPUT
 !        | c(1)  c(2) c(3) c(4)  | EX
 !   Ce = | c(2)  c(5) c(6) c(7)  | EY     Cs = | c(11) c(12) ! GX3
 !        | c(3)  c(6) c(8) c(9)  |2EXY         |       c(13) ! GY3
 !        | c(4)  c(7) c(9) c(10) | E3
 REAL (kind=8) :: c(:),rr(:),ce(:)
 LOGICAL, OPTIONAL :: flag
 !local
 REAL (kind=8) c2,c4,s2,s4,c4s4,sc
 !     rotate the stress-strain matrix
 !     t(transpose) * c(material)

 IF( rr(1) < 0.9999999999 )THEN

   c2 = rr(1)          !cos^2(t)
   c4 = c2*c2          !cos^4(t)
   s2 = rr(2)          !sin^2(t)
   s4 = s2*s2          !sin^4(t)
   c4s4 = c4+s4        !cos^4(t)+sin^4(t)
   sc = rr(3)          !cos(t)*sin(t)

   ce(1) = c(1)*c4 + (c(2) +2*c(7))*2*s2*c2 + c(4)*s4                 !c(1,1)
   ce(2) = (c(1)+c(4)-4*c(7))*s2*c2 + c(2)*(s4+c4)                    !c(1,2)
   ce(3) = (c(1)-c(2)-2*c(7))*sc*c2 + (c(2)-c(4)+2*c(7))*s2*sc        !c(1,4)
   ce(4) =  c(3)*c2 + c(5)*s2                                         !c(1,3)
   ce(5) = c(1)*s4 + (c(2) +2*c(7))*2*s2*c2 + c(4)*c4                 !c(2,2)
   ce(6) = (c(1)-c(2)-2*c(7))*sc*s2 + (c(2)-c(4)+2*c(7))*c2*sc        !c(2,4)
   ce(7) =  c(5)*c2 + c(3)*s2                                         !c(2,3)
   ce(8) = (c(1)+c(4)-2*c(2)-2*c(7))*s2*c2 + c(7)*(s4+c4)             !c(4,4)
   ce(9) = (c(5)-c(3))*sc                                             !c(3,4)
   ce(10)=  c(6)                                                      !c(3,3)
   IF( PRESENT(flag) )THEN
     ce(11) = c(8)*c2+c(9)*s2        !c(5,5)
     ce(12) = (c(9)-c(8))*sc         !c(5,6)
     ce(13) = c(8)*s2+c(9)*c2        !c(6,6)
   END IF
 ELSE
   ce(1) = c(1)                                                       !c(1,1)
   ce(2) = c(2)                                                       !c(1,2)
   ce(3) = 0d0                                                        !c(1,4)
   ce(4) = c(3)                                                       !c(1,3)
   ce(5) = c(4)                                                       !c(2,2)
   ce(6) = 0d0                                                        !c(2,4)
   ce(7) = c(5)                                                       !c(2,3)
   ce(8) = c(7)                                                       !c(4,4)
   ce(9) = 0d0                                                        !c(3,4)
   ce(10)= c(6)                                                       !c(3,3)
   IF( PRESENT(flag) )THEN
     ce(11:13) = (/ c(8), 0d0, c(9) /)
   END IF
 END IF
 RETURN
 END SUBROUTINE rotor3d
