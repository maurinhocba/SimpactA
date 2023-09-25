 SUBROUTINE sc_shel5( sec )
 !            --------------------------------
 !            input material properties for multi-layered lamina
 !            Total Lagrangean Shear deformable elements (SHELQ,SHELT)
 !  iprop: 1 :nlayr
 !         2 :nvars:  number of internal variables for the whole laminate (NOT YET)
 !         3 :0  elastic lamina,  1:inelastic   (0 ELASTIC ONLY)
 !         4 :
 !    5:4+nl :associated material label to each layer
 !
 !  rprop  1 : thick
 !         2 : mean young modulus
 !         3 : mean density
 !         4 : shear factor (5/6 is the standard value)
 !         5 : surface density
 !      6-29 : (24) Standard integrated elasticity matrix
 !     30-63 : (34) in-plane integrated matrices associated with zig-zag function
 !     64-70 : ( 7) transverse integrated matrices associated with zig-zag function
 !
 !             6  7  8 21 22 23        40 41 42 43
 !                9 10 24 25 26        44 45 46 47
 !                  11 27 28 29        48 49 50 51
 !                     12 13 14        52 53 54 55
 !                        15 16        56 57 58 59
 !                           17        60 61 62 63
 !                              18 19              67 68
 !                                 20              69 70
 !                                     30 31 32 33
 !                                        34 35 36
 !                                           37 38
 !                                              39
 !                                                 64 65
 !                                                    66
 !                                                       71 72
 !                                                          73
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

 INTEGER (kind=4) matno,nlayr,ith,il,sizer,sizei,ir,ii,j,k,nvars,kvar,imat
 LOGICAL :: found,read_integrated_matrix,fastr
 REAL(kind=8) :: thick,youn,dens,thickl,de(65),angle,zl,ce(9),sg,cg,celoc(6),gi(2),g1, &
                 f11,f22,f12,f1,f2,sf,im(0:2),thicks, &
                 dbi(3,3),db(3,3),h(2,2),a(3,3),d1(2,2),ds(2,2),cz(3,3),det,csi(2,2)
 REAL(kind=8), ALLOCATABLE :: rprop(:,:),phi(:,:),beta(:,:)
 TYPE (mater), POINTER :: mat
 TYPE (postv), POINTER :: var1,var2 !pointers to postprocess variables

 sec%secty = 15         !shell type 5
 sec%mabas = -1         ! flag to indicate that no unique associated material exist

 CALL listen('SHELL5')  !read a line

 sizer = 76             ! minimum size of real array
 thick = getrea('THICK ',0.D0,' Thickness ........................')
 IF( thick > 0d0 ) THEN
   ith=1 ! thickl will be input as a percentage of thick
   thicks = 0d0 !initializes thickness to check
 ELSE
   ith=0 ! thickl will be input as a Real dimension
   thick=0d0   !Initializes thickness
 END IF
 nlayr  = getint('NLAYR ',   0,'!NUMBER OF LAYERS .................')
 IF ( nlayr <= 0) CALL runend('SHELL5: VERIFY NUMBER OF LAYERS ')
 sizei =  4+nlayr*2                !size of integer array
 ALLOCATE (sec%iprop(sizei), rprop(npl,nlayr) )
 sec%iprop = 0   !initializes
 rprop = 0d0     !initializes
 sf = getrea('SHEARF',0.00000D0,' Shear Factor .....................') !for FSDT
 read_integrated_matrix = exists('READMA')

 IF( exists('TSHEAR'))sec%iprop(4) = 1

 youn = 0d0            !Initializes mean young modulus
 dens = 0d0            !Initializes mean density
 ii  = 5               !pointer to integer array
 kvar=1                !pointer to layer internal variables
 gi = 0d0              !average values of trasnverse modulus
 im = 0d0              !rot inertia

 DO il=1,nlayr         !for each layer
   ! read layer information
   CALL listen('SHELL5')  !read a line
   matno    = getint('NUMAT ', 1 ,'!Associated Material ..............')
   CALL mat_search(matno,found,mat,imat)        !search base material
   IF(.NOT.found) CALL runend ('SHELL5: NOT FOUND ASSOCIATED MATERIAL ')
   IF( .NOT.ANY( mat%mtype == matsa) ) &  !check if base material is correct
            CALL runend('SHELL5: Invalid material association')
   sec%iprop(ii)       = imat     !store associated material index (in PMATS)
   sec%iprop(ii+nlayr) = kvar     !pointer to layer varin (internal variables)
   IF(ith == 0) THEN   !thicknes must be computed
     thickl = getrea('THICKL',1.D0,'!LAYER THICKNESS ..................')
     thick  = thick + thickl         ! sum layer thickness
   ELSE                !thickness is known
     thickl = getrea('THICKL',1.D+2/nlayr,' LAYER THICKNESS ( % ).............')
     thickl = thick * thickl *1.D-2   ! compute layer thickness
     thicks = thicks + thickl
   END IF
   angle = getrea('ANGLE ',0.D0, ' ANGLE FROM COMP COORDS TO LA MAT 1')
   youn = youn + thickl*mat%prope(1)  !mean young modulus
   dens = dens + thickl*mat%prope(5)  !surface density
   rprop(1,il) = thickl               !store layer thickness (real dimension)
   rprop(2,il) = angle                !store layer angle     (degrees)
   ! rprop(3,il) = getrea('POSIT',0.0d0,' LAYER POSITION ..................') !no sense
   angle = angle*factor       !angle between orthotropic axes and global axes in rad
   sg = SIN(angle)            !sine
   cg = COS(angle)            !cosine
   ! matrix to transform global strains into orthotropic strains
   rprop(4,il) =  cg**2                ! [ c^2     s^2     c s       ]
   rprop(5,il) =  sg**2                ! [ s^2     c^2    -c s       ]
   rprop(6,il) =  cg*sg                ! [-2 c s   2 c s   c^2 - s^2 ]
   rprop(7,il) =  2d0*cg*sg
   rprop(8,il) =  cg**2 - sg**2
   !get required internal variables increment
   kvar = kvar + get_nvar(mat,1)  !1: plane stress
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
   CALL runend('SHELL5: Sum of layers /= 100%')
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
 nvars = kvar-1                         !number of internal variables
 sec%iprop(1) = nlayr                   !keep number of layers
 sec%iprop(2) = nvars                   !set total number of internal variables
 sec%iprop(3) = 0                       !initializes to elastic lamina

 !IF( nvars > 0 )THEN      ! if internal variables exist
 !  sec%iprop(3) = 1                     !set to plastic lamina
 !END IF
 sizer = sizer+nlayr*npl                !update size of real array
 sec%secdef = (/ sizei,sizer, 0,0,0 /)  !size of arrays
 ALLOCATE (sec%rprop(sizer) )           !get memory
 sec%rprop = 0d0
 sec%rprop(1) = thick                   !thickness
 sec%rprop(2) = youn/thick              !mean young modulus
 sec%rprop(3) = dens/thick              !mean density
 sec%rprop(4) = sf                      !keep shear factor
 sec%rprop(5) = dens                    !surface density

 IF(.NOT.read_integrated_matrix)THEN
   ! compute elastic integrated matrix
   de = 0d0                               !initializes integrated thickness
   zl = -thick/2d0      !initializes Z position at bottom surface
   DO il=1,nlayr                           !for each layer
     imat= sec%iprop(4+il)                   !material index
     mat => pmats(imat)%p                    !point to material
     ! celoc(1:4) in plane values  celoc(5:6) shear values
     CALL elastiff(celoc,mat,4)              !get elastic stiffnes matrix 4=plain stress + transverse shear
     thickl = rprop(1,il)                          !layer thickness
     g1 = thickl*thickl/12d0                       !auxiliar h^2/12
     zl = zl+thickl/2d0                          !update Z position to layer mid-point z_k-bar
     rprop(3,il) = zl                            !keep
     CALL rotortm(celoc,rprop(4:,il),ce,.TRUE.)  !rotate elastic in-plane matrix
     beta(1,il) = gi(1)/ce(7) - 1d0
     beta(2,il) = gi(2)/ce(9) - 1d0
     phi(1,il) = phi(1,il-1) + beta(1,il)*thickl
     phi(2,il) = phi(2,il-1) + beta(2,il)*thickl
     f1 = (phi(1,il) + phi(1,il-1))/2d0          !phi_k-bar 1
     f2 = (phi(2,il) + phi(2,il-1))/2d0          !phi_k-bar 2
     im(0) = im(0) + thickl*mat%prope(5)*(zl**2+thickl**2/12d0)      !rotatory inertia
     im(1) = im(1) + thickl*mat%prope(5)*(f1**2+(beta(1,il)*thickl)**2/12d0)   !additional rotatory inertia
     im(2) = im(2) + thickl*mat%prope(5)*(f2**2+(beta(2,il)*thickl)**2/12d0)
     ce = ce*thickl                              !scale layer stiffness matrix
     !  components of the integrated stiffness matrix
     !             1  2  3 16 17 18        35 36 37 38
     !             2  4  5 19 20 21        39 40 41 42                 | 1  2  3  |
     !             3  5  6 22 23 24        43 44 45 46            Ce = | 2  4  5  |
     !            16 19 22  7  8  9        47 48 49 50                 | 3  5  6  |
     !            17 20 23  8 10 11        51 52 53 54
     !            18 21 24  9 11 12        55 56 57 58
     !                              13 14              62 63
     !                              14 15              64 65
     !            35 39 43 47 51 55        25 26 27 28
     !            36 40 44 48 52 56        26 29 30 31
     !            37 41 45 49 53 57        27 30 32 33
     !            38 42 46 50 54 58        28 31 33 34
     !                              62 64              59 60
     !                              63 65              60 61
     !standard part
     de( 1: 6) = de( 1: 6) + ce(1:6)               !membrane part
     de( 7:12) = de( 7:12) + ce(1:6)*(zl*zl+g1)    !bending part
     de(13:15) = de(13:15) + ce(7:9)               !shear part
     !! Barbero (pag 64 ), not consistent with present approach, includes the 5/6 factor
     !de(13:15) = de(13:15) + ce(7:9)*(1.25d0 - 5d0/thick**2*(zl*zl+thickl**2/12))   !shear part
     de(16:18) = de(16:18) + ce(1:3)*zl            !coupled part
     de(19:21) = de(19:21) + (/ce(2),ce(4:5)/)*zl  !coupled part
     de(22:24) = de(22:24) + (/ce(3),ce(5:6)/)*zl  !coupled part
     ! Dp(3,3)  phi^T Dk phi                !25 26 27 28
     f11 = f1*f1+beta(1,il)**2*g1           !   29 30 31
     f12 = f1*f2+beta(1,il)*beta(2,il)*g1   !      32 33
     f22 = f2*f2+beta(2,il)**2*g1           !         34
     de(25) = de(25) + ce(1)*f11
     de(26) = de(26) + ce(2)*f12
     de(27) = de(27) + ce(3)*f11
     de(28) = de(28) + ce(3)*f12
     de(29) = de(29) + ce(4)*f22
     de(30) = de(30) + ce(5)*f12
     de(31) = de(31) + ce(5)*f22
     de(32) = de(32) + ce(6)*f11
     de(33) = de(33) + ce(6)*f12
     de(34) = de(34) + ce(6)*f22
     ! Dp(1,3)                                     !35 36 37 38
                                                   !39 40 41 42
                                                   !43 44 45 46
     de(35) = de(35) + ce(1)*f1
     de(36) = de(36) + ce(2)*f2
     de(37) = de(37) + ce(3)*f1
     de(38) = de(38) + ce(3)*f2
     de(39) = de(39) + ce(2)*f1
     de(40) = de(40) + ce(4)*f2
     de(41) = de(41) + ce(5)*f1
     de(42) = de(42) + ce(5)*f2
     de(43) = de(43) + ce(3)*f1
     de(44) = de(44) + ce(5)*f2
     de(45) = de(45) + ce(6)*f1
     de(46) = de(46) + ce(6)*f2
     !Dp(2,3)                                      !47 48 49 50
     f11 = f1*zl+beta(1,il)*g1                     !51 52 53 54
     f22 = f2*zl+beta(2,il)*g1                     !55 56 57 58
     de(47) = de(47) + ce(1)*f11
     de(48) = de(48) + ce(2)*f22
     de(49) = de(49) + ce(3)*f11
     de(50) = de(50) + ce(3)*f22
     de(51) = de(51) + ce(2)*f11
     de(52) = de(52) + ce(4)*f22
     de(53) = de(53) + ce(5)*f11
     de(54) = de(54) + ce(5)*f22
     de(55) = de(55) + ce(3)*f11
     de(56) = de(56) + ce(5)*f22
     de(57) = de(57) + ce(6)*f11
     de(58) = de(58) + ce(6)*f22
     ! Dt(2,2)                                     !59 60
     f1 = beta(1,il)                               !   61
     f2 = beta(2,il)
     f11 = f1*f1
     f12 = f1*f2
     f22 = f2*f2
     de(59) = de(59) + ce(7)*f11
     de(60) = de(60) + ce(8)*f12
     de(61) = de(61) + ce(9)*f22
     ! Dt(1,2)                                     !62 63
     de(62) = de(62) + ce(7)*f1                    !64 65
     de(63) = de(63) + ce(8)*f2
     de(64) = de(64) + ce(8)*f1
     de(65) = de(65) + ce(9)*f2
     !inertia values
     zl = zl+thickl/2d0                          !update Z position to layer interface
   END DO
   WRITE(58,"(' phi and Beta Values')")
   WRITE(58,"(2e15.8)")phi(:,0)
   DO il=1,nlayr
     WRITE(58,"(5x,2e15.8)")beta(:,il)
     WRITE(58,"(2e15.8)")phi(:,il)
   END DO
   WRITE(58,"(' Integrated D matrices')")
   WRITE(58,"(5e15.8)")de
   WRITE(58,"(' Integrated Inertias')")
   WRITE(58,"(5e15.8)")im
 ELSE
   fastr = .FALSE.
   CALL rdfrre('SC_SHELL5',de,il,65,fastr)
   CALL rdfrre('SC_SHELL5',im,il,3,fastr)
 END IF
 sec%rprop(6:70) = de                          !keep integrated matrix
 sec%rprop(74:76) = im                         !keep integrated inertias
 ! NO ESTOY INCLUYENDO LAS FUNCIONES PHI NI SUS DERIVADAS BETA
 IF( nvars>0 ) THEN      !for inelastic behavior
   ir = 27            !initializes pointer to layer thickness
   ii = ir + nlayr    !initializes pointer to angle
   j  = ii + nlayr    !initializes pointer to z
   k  = j  + nlayr    !initializes pointer to rotation matrices
   DO il=1,nlayr    !for each layer
     sec%rprop(ir) = rprop(1,il)               !layer thickness
     sec%rprop(ii) = rprop(2,il)               !layer angle
     sec%rprop( j) = rprop(3,il)               !layer position z
     CALL vecasi( 5,rprop(4,il),sec%rprop(k))  !rotation matrices
     ! update pointers
     ir = ir+ 1   !layer thickness
     ii = ii+ 1   !angle
     j  = j + 1   !z
     k  = k + 5   !rotation matrices
   END DO
 END IF
 !  compute Shear correction factors
 db(:,1) = sec%rprop(12:14) ; db(:,2) = sec%rprop((/13,15:16/)) ; db(:,3) = sec%rprop((/14,16:17/))
 CALL invmtx(db,dbi,det,3)
 h = 0d0
 ds = 0d0
 zl = -thick/2d0      !initializes Z position at bottom surface
 DO il=1,nlayr                           !for each layer
   imat= sec%iprop(4+il)                   !material index
   mat => pmats(imat)%p                    !point to material
   thickl = rprop(1,il)                          !layer thickness
   CALL elastiff(celoc,mat,4)              !get elastic stiffnes matrix 4=plain stress
   CALL rotortm(celoc,rprop(4:,il),ce,.TRUE.)     !rotate elastic in-plane matrix
   cz(:,1) = ce((/1:3/)) ; cz(:,2) = ce((/2,4:5/)) ; cz(:,3) = ce((/3,5:6/))
   a = MATMUL(cz,dbi)*thickl
   d1 = RESHAPE( (/ a(1,1)+a(3,3), a(3,1)+a(2,3), a(1,3)+a(3,2), a(2,2)+a(3,3) /),(/2,2/))
   det = ce(7)*ce(9) - ce(8)*ce(8)
   csi = RESHAPE( (/ ce(9), -ce(8), -ce(8), ce(7) /),(/2,2/)) /det
   DO k=1,5  !use five-point integration for each layer
     zl = zl+thickl/10.                         !update Z position to layer mid-point z_k-bar
     ds = ds-zl/2d0*d1/5d0                      !
     h  = h + MATMUL(MATMUL(ds,csi),ds)*thickl/5d0
     zl = zl+thickl/10.                         !update Z position to layer mid-point z_k-bar
   END DO
 END DO                                   ! z = layer midpoint
 det = h(1,1)*h(2,2) - h(1,2)*h(2,1)
 csi = RESHAPE( (/ h(2,2), -h(2,1), -h(1,2), h(1,1) /),(/2,2/)) /det
 ce(7:9) = sec%rprop(18:20)
 sec%rprop(71) = csi(1,1)/ce(7)
 sec%rprop(73) = csi(2,2)/ce(9)
 IF(ce(8) > 1d-16 )THEN
     sec%rprop(72) = csi(1,2)/ce(8)
 ELSE
     sec%rprop(72) = (sec%rprop(71)+sec%rprop(73))/2d0
 END IF
 WRITE(lures,"(' SHEAR CORRECTION FACTORS ',/,' Kx: ',f10.5,/,' Kxy:',f10.5,/,' Ky: ',f10.5)")sec%rprop(71:73)
 !  end of computation of shear correction factors
 DEALLOCATE (rprop,phi,beta)                 !release local memory
 sec%secdef(3) =  5        !GENERAL orthotropic with COUPLING

 ! Additional postprocess variables --------------------------------------------
 sec%secdef(4) = 2  !7 variables
 sec%secdef(5) = 6  !7 values

 ALLOCATE(var1,var2)  ! new variables

 var1%type = 2
 var1%dim = 2
 var1%name(:) = (/ 'Phi_Moment','Phi_Mxx','Phi_Myy','Phi_Mxy','Phi_Myx',' ',' ' /)

 var2%type = 1
 var2%dim = 2
 var2%name(:) = (/ 'Phi_Shear','Phi_Shear_X','Phi_Shear_y',' ',' ',' ',' ' /)

 sec%postp => var1   !keep pointer to first variable
 var1%next => var2
 NULLIFY (var2%next)

 !------------------------------------------------------------------------------
 RETURN
 END SUBROUTINE sc_shel5
!********************************************************************
