 SUBROUTINE delt24(nelem ,heade ,deltc , coord )

 !     Compute Critical time increment for a set of elements BST-CST (TLF)
 !     includes Young modulus temperature dependance

 USE ctrl_db, ONLY : mscal
 IMPLICIT NONE
 INTEGER (kind=4), INTENT (IN) :: nelem  !number of elements in the set
 REAL (kind=8), INTENT(IN) :: coord(:,:) !present mesh coordinates
 REAL (kind=8), INTENT(IN OUT) :: deltc     !critical time increment
 TYPE (ele24), POINTER :: heade          !pointer to first element

 INTEGER :: ielem,in,jn,kn,isec,osec,secty,i,nyoung
 REAL (kind=8) :: lmax(4),t1(3),t2(3),young,densi,fac, &
                  eigen,a1,a2,b2,poiss,c,r,tcrit,aux,fac1,thick
 TYPE (ele24), POINTER :: elem
 TYPE (section), POINTER :: sec
 TYPE (curve), POINTER :: cyoung

 osec = -1                             !previous material number
 lmax = 0d0                            !initializes longest side
 elem => heade                         !point to first element
 DO ielem=1,nelem                      !for each element in the set
   isec = elem%matno                  !associated material
   IF( isec /= osec )THEN             !present and previous material diff
     sec => psecs(isec)%p
     secty  = sec%secty
     poiss  = sec%mtbas%prope(2)      !poisson ratio (may be temperature dependant)
     densi  = sec%mtbas%prope(5)*mscal      !mass density
     thick = sec%rprop(1)             !original thickness
     nyoung= INT(sec%mtbas%props(1))        !associated curve for Young modulus
     IF( nyoung > 0 )THEN
       fac1  = 1d0/densi/(1d0-poiss**2)  !factor
       CALL cur_point( nyoung,sec%mtbas%chead,cyoung) !point to the curve
       i = 1
     ELSE
       young  = sec%mtbas%prope(1)      !young modulus (may be temperature dependant)
       fac   = young/densi/(1d0-poiss**2)!factor
     END IF
     osec  = isec                     !update previous material
   END IF
   IF( nyoung > 0 )THEN
     aux = SUM(tempe(1,elem%lnods(1:3)))/3d0         !average temperature (center or bottom)
     young = inte_cr(cyoung%val,nyoung,aux,i)
     fac = young*fac1
   END IF
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
   elem => elem%next                         !point to next element
 END DO

! WRITE(55,*)SQRT(2d0/lmax(1)),SQRT(2d0/lmax(2)),SQRT(2d0/lmax(3))

 tcrit = SQRT(2d0/MAXVAL(lmax(1:3)))     !critical time increment
 IF( lmax(4) > 1d0 )tcrit = tcrit/lmax(4)
 IF( tcrit < deltc ) deltc = tcrit

 RETURN
 END SUBROUTINE delt24
