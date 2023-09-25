 SUBROUTINE slum13 (nelem,heade,ymass,coord,dtscal)

 ! compute scaled lumped mass vector for element BST-LST (TLF) shell element

 IMPLICIT NONE
 !  dummy arguments
 INTEGER (kind=4), INTENT (IN) :: nelem  !number of elements in the set
 REAL (kind=8), INTENT(IN) :: dtscal
 REAL (kind=8) :: ymass(:,:),      & !nodal mass in each direction
                  coord(:,:)       !nodal coordinates
 TYPE (ele13), POINTER  :: heade     !pointer to first element
 ! local variables
 INTEGER (kind=4) :: ielem,in,jn,kn,isec,osec,mtype,i
 INTEGER (kind=4) :: lnode,idofn,inode
 REAL (kind=8) :: diagm,voli
 REAL (kind=8) :: lmax(4),t1(3),t2(3),young,densi,fac, &
                  a1,a2,b2,poiss,c,r,tcrit,thick,      &
                  fthic,x(3,3),l(3),d11,d22,d12,fact  !
 TYPE (ele13), POINTER :: elem
 TYPE (section), POINTER :: sec
 LOGICAL :: isotro

 fact = 2d0*dtscal**2               !factor
 osec = -1                            !previous material number
 elem => heade                        !point to first element
 DO ielem=1,nelem                     !for each element in the set
   isec = elem%matno                    !associated section
   IF( isec /= osec )THEN               !present and previous material differs
     sec => psecs(isec)%p                 !point to section
     thick = sec%rprop(1)                 !original thickness
     IF( sec%secty == 12)THEN             !standard section
       mtype  = sec%mtbas%mtype             !associated material
       isotro = mtype == 1 .OR. mtype == 6  !isotropic material => true
       IF( isotro )THEN                     !if material is isotropic !
         young  = sec%mtbas%prope(1)          !young modulus
         poiss  = sec%mtbas%prope(2)          !poisson ratio
         fac    = young/(1d0-poiss**2)        !plane stress factor
       ELSE
         densi = 0.75d0*thick               !surface mass density
         d11 = sec%rprop(5)/densi           !modified plane stress matrix
         d22 = sec%rprop(7)/densi
         d12 = 2d0*(sec%rprop(6)+sec%rprop(8)) /densi
         fthic = thick*thick                !thickness squared
         fac = fthic**3                     !thickness^6
         lmax(4) = 1d0                      !initializes
       END IF
     ELSE      !composite section
       isotro = .FALSE.             !anisotropic
       densi = 0.75d0*thick         !surface mass density
       d11 = sec%rprop(6)/densi     !modified plane stress matrix
       d22 = sec%rprop(9)/densi
       d12 = 2d0*(sec%rprop(7)+sec%rprop(11)) /densi
       fthic = thick*thick                !thickness squared
       fac = fthic**3                     !thickness^6
       lmax(4) = 1d0                      !initializes
     END IF
     osec  = isec                     !update previous material
   END IF
   IF( isotro )THEN             !isotropic material
     in = elem%lnods(1)                  !first node of the element
     jn = elem%lnods(2)                  !second node of the element
     kn = elem%lnods(3)                  !third node of the element
     t1 = coord(1:3,jn) - coord(1:3,in)  !first side 1-2
     t2 = coord(1:3,kn) - coord(1:3,jn)  !second side 1-3
     a1 = SQRT(DOT_PRODUCT(t1,t1))       !first side length (base length)
     a2 = DOT_PRODUCT(t2,t1)/a1          !proyection of second side over first
     b2 = SQRT(DOT_PRODUCT(t2,t2)-a2**2) !triangle heigth
     ! compute membrane eigenvalues
     c = (b2**2+a2**2+a1*a2+a1**2)/(a1*b2)**2            !center
     r = SQRT(2*a1**3*a2-a1**2*b2**2+a1**4+b2**4        &!radius
                +a2**4+3*poiss**2*a1**2*b2**2           &
                +2*b2**2*a1*a2+3*a2**2*a1**2+2*a2**3*a1 &
                +2*a2**2*b2**2)/(a1*b2)**2
     lmax(1) = (c-r)*fac                           !first eigenvalue
     lmax(2) = (c+r)*fac                           !second eigenvalue
     lmax(3) = c*young/(1d0+poiss)                 !third eigenvalue
     lmax(4) = elem%lb*thick/MIN(a1,b2)            !thickness influence
   ELSE    !anisotropic material
     x(:,1:3) = coord(:,elem%lnods(1:3)) ! get coordinates of the main triangle
     ! large strain version (does not include thickness change)
     !g(:,1) = -MATMUL(x,elem%b)     ! x(1)
     !g(:,2) = +MATMUL(x,elem%a)     ! x(2)
     !CALL vecpro(g(1,1),g(1,2),g(1,3))
     !CALL vecpro(g(1,2),g(1,3),t1(1))
     !CALL vecpro(g(1,3),g(1,1),t2(1))
     !t1 = t1*elem%lb
     !t2 = t2*elem%lb
     ! small strain version
     t1 = -MATMUL(x,elem%b)     ! x(1)
     t2 = +MATMUL(x,elem%a)     ! x(2)
     !
     DO i=1,3
       l = x(:,hh(2,i)) - x(:,hh(3,i))
       a1 = DOT_PRODUCT(t1,l)  !component along local coord 1
       a2 = DOT_PRODUCT(t2,l)  !component along local coord 2
       b2 = (a1**2+a2**2)**3   !length^6
       c  = (d11*a1**4+d12*a1**2*a2**2+d22*a2**4)/b2
       IF( b2 < fac ) c = c*fthic/b2**(0.3333d0)
       lmax(i) = c
     END DO
   END IF
   tcrit = fact/MAXVAL(lmax(1:3))     !critical time increment
   IF( lmax(4) > 1d0 )tcrit = tcrit/lmax(4)**2
   densi = 1d0/tcrit
   voli = elem%area1*densi*thick         !element mass
   diagm = voli/3.d0                     !nodal mass
   !
   ! ---generates lumped mass matrix proportional to diagonal
   !
   DO inode=1,3                          !for each element node
     lnode=elem%lnods(inode)             !node number
     DO idofn=1,3                        !translational DOFs only
       ymass(idofn,lnode)=ymass(idofn,lnode)+diagm
     END DO
   END DO
   elem => elem%next                     !point to next element
 END DO

 RETURN
 END SUBROUTINE slum13
