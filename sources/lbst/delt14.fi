 SUBROUTINE delt14(nelem ,heade ,deltc , coord )

 !     Compute Critical time increment for a set of elements BST-CST (TLF)

 USE ctrl_db, ONLY : mscal
 IMPLICIT NONE
 INTEGER (kind=4), INTENT (IN) :: nelem  !number of elements in the set
 REAL (kind=8), INTENT(IN) :: coord(:,:) !present mesh coordinates
 REAL (kind=8), INTENT(IN OUT) :: deltc     !critical time increment
 TYPE (ele14), POINTER :: heade          !pointer to first element

 INTEGER :: ielem,in,jn,kn,isec,osec,mtype,i
 REAL (kind=8) :: lmax(4),t1(3),t2(3),young,densi,fac, &
                  eigen,a1,a2,b2,poiss,c,r,tcrit,thick,&
                  x(3,3),l(3),d11,d22,d12
 TYPE (ele14), POINTER :: elem
 TYPE (section), POINTER :: sec
 LOGICAL :: isotro

 osec = -1                             !previous material number
 lmax = 0d0                            !initializes longest side
 lmax(4) = 1d0
 elem => heade                         !point to first element
 DO ielem=1,nelem                      !for each element in the set
   isec = elem%matno                  !associated material
   IF( isec /= osec )THEN             !present and previous material diff
     sec => psecs(isec)%p
     thick = sec%rprop(1)
     IF( sec%secty == 12)THEN
       mtype  = sec%mtbas%mtype
       isotro = mtype == 1 .OR. mtype == 6
       IF( isotro )THEN
         young  = sec%mtbas%prope(1)      !young modulus
         poiss  = sec%mtbas%prope(2)      !poisson ratio
         densi  = sec%mtbas%prope(5)*mscal      !mass density
         fac   = young/densi/(1d0-poiss**2)!factor
       ELSE
         isotro = .FALSE.
         densi = 0.75d0*sec%mtbas%prope(5)*thick*mscal     !surface mass density
         d11 = sec%rprop(5)/densi
         d22 = sec%rprop(7)/densi
         d12 = 2d0*(sec%rprop(6)+sec%rprop(8)) /densi
         thick = thick*thick
         fac = thick**3                 !thickness^6
       END IF
     ELSE
       isotro = .FALSE.
       densi = 0.75d0*sec%rprop(5)*mscal     !surface mass density
       d11 = sec%rprop(6)/densi
       d22 = sec%rprop(9)/densi
       d12 = 2d0*(sec%rprop(7)+sec%rprop(11)) /densi
       thick = thick*thick
       fac = thick**3                 !thickness^6
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
     eigen=(c-r)*fac                           !first eigenvalue
     IF(eigen > lmax(1)) lmax(1) = eigen       !compare with previous max.
     eigen=(c+r)*fac                           !second eigenvalue
     IF(eigen > lmax(2)) lmax(2) = eigen       !compare with previous max.
     eigen=c*young/densi/(1d0+poiss)           !third eigenvalue
     IF(eigen > lmax(3)) lmax(3) = eigen       !compare with previous max.
     !  thickness influence
     eigen = elem%lb*thick/MIN(a1,b2)
     IF( eigen > lmax(4) )lmax(4) = eigen
   ELSE    !anisotropic material
     x(:,1:3) = coord(:,elem%lnods(1:3)) ! get coordinates of the main triangle
     ! large strain version (does not include thickness change)
     !g(:,1) = -MATMUL(x,elem%b(:,0))     ! x(1)
     !g(:,2) = +MATMUL(x,elem%a(:,0))     ! x(2)
     !CALL vecpro(g(1,1),g(1,2),g(1,3))
     !CALL vecpro(g(1,2),g(1,3),t1(1))
     !CALL vecpro(g(1,3),g(1,1),t2(1))
     !t1 = t1*elem%lb
     !t2 = t2*elem%lb
     ! small strain version
     t1 = -MATMUL(x,elem%b(:,0))     ! x(1)
     t2 = +MATMUL(x,elem%a(:,0))     ! x(2)
     !
     DO i=1,3
       l = x(:,hh(1,i)) - x(:,hh(2,i))
       a1 = DOT_PRODUCT(t1,l)  !component along local coord 1
       a2 = DOT_PRODUCT(t2,l)  !component along local coord 2
       b2 = (a1**2+a2**2)**3   !length^6
       c  = (d11*a1**4+d12*a1**2*a2**2+d22*a2**4)/b2
       IF( b2 < fac ) c = c*thick/b2**(0.3333d0)
       IF( c > lmax(i) )lmax(i) = c
     END DO
   END IF
   elem => elem%next                         !point to next element
 END DO

! WRITE(55,*)SQRT(2d0/lmax(1)),SQRT(2d0/lmax(2)),SQRT(2d0/lmax(3))

 tcrit = SQRT(2d0/MAXVAL(lmax(1:3)))     !critical time increment
 IF( lmax(4) > 1d0 )tcrit = tcrit/lmax(4)
 IF( tcrit < deltc ) deltc = tcrit

 RETURN
 END SUBROUTINE delt14
