 SUBROUTINE delt25(nelem ,heade ,deltc , coord )

 !     Compute Critical time increment for a set of elements BSQ (TLF)

 USE ctrl_db, ONLY : mscal
 IMPLICIT NONE
 INTEGER (kind=4), INTENT (IN) :: nelem  !number of elements in the set
 REAL (kind=8), INTENT(IN) :: coord(:,:) !present mesh coordinates
 REAL (kind=8), INTENT(IN OUT) :: deltc     !critical time increment
 TYPE (ele25), POINTER :: heade          !pointer to first element

 INTEGER :: ielem,i,isec,osec,mtype
 REAL (kind=8) :: lmin,ts(3,4),young,densi,fac,     &
                  poiss,tcrit,thick,t1(3),t2(3),    &
                  a1,a2,b2,x(3,4),ls(4),d11,d22,d12

 TYPE (ele25), POINTER :: elem
 TYPE (section), POINTER :: sec
 LOGICAL :: isotro


 osec = -1                             !previous material number
 lmin = 0d0                            !initializes longest side
 elem => heade                         !point to first element
 tcrit = 1e20
 DO ielem=1,nelem                      !for each element in the set
   isec = elem%matno                   !associated material
   IF( isec /= osec )THEN              !present and previous material diff
     sec => psecs(isec)%p              !point to section
     thick = sec%rprop(1)**2           !thickness squared
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
         densi = 0.75d0*sec%mtbas%prope(5)*sec%rprop(1)*mscal     !surface mass density
         d11 = sec%rprop(5)/densi
         d22 = sec%rprop(7)/densi
         d12 = 2d0*(sec%rprop(6)+sec%rprop(8)) /densi
         fac = thick**3                 !thickness^6
       END IF
     ELSE
       isotro = .FALSE.
       densi = 0.75d0*sec%rprop(5)*mscal     !surface mass density
       d11 = sec%rprop(6)/densi
       d22 = sec%rprop(9)/densi
       d12 = 2d0*(sec%rprop(7)+sec%rprop(11)) /densi
       fac = thick**3                 !thickness^6
     END IF
     osec  = isec                     !update previous material
   END IF
   x = coord(:,elem%lnods(1:4))
   DO i=1,4
     ts(:,i) = x(:,nextn(i)) - x(:,i)        !side
     ls(i) = DOT_PRODUCT(ts(:,i),ts(:,i))    !side length (squared)
   END DO
   IF( isotro )THEN             !isotropic material
     lmin = MINVAL(ls)                   !critical time (squared)
     IF( lmin < thick ) lmin = lmin* lmin/thick
     lmin = lmin/fac                     !critical time (squared)
     !IF(lmin < tcrit) tcrit = lmin       !compare with previous min
     IF(lmin < tcrit)THEN
       tcrit = lmin       !compare with previous min
     END IF 
   ELSE                         !anisotropic material
     DO i=1,4 !for each side
       ! large strain version (does not include thickness change)
       !g(:,1) = MATMUL(x,elem%cartd(:,1,i))     ! x(1)
       !g(:,2) = MATMUL(x,elem%cartd(:,2,i))     ! x(2)
       !CALL vecpro(g(1,1),g(1,2),g(1,3))
       !CALL vecpro(g(1,2),g(1,3),t1(1))
       !CALL vecpro(g(1,3),g(1,1),t2(1))
       !t1 = t1*elem%lb
       !t2 = t2*elem%lb
       ! small strain version
       t1 = MATMUL(x,elem%cartd(:,1,i))         ! x(1)
       t2 = MATMUL(x,elem%cartd(:,2,i))         ! x(2)
       a1 = DOT_PRODUCT(t1,ts(:,i))  !component along local coord 1
       a2 = DOT_PRODUCT(t2,ts(:,i))  !component along local coord 2
       b2 = ls(i)**3                 !length^6
       lmin  = b2/(d11*a1**4+d12*a1**2*a2**2+d22*a2**4)
       IF( b2 < fac ) lmin = lmin*b2**(0.3333d0)/thick
       IF(lmin < tcrit) tcrit = lmin       !compare with previous min
     END DO
   END IF
   elem => elem%next                         !point to next element
 END DO

 tcrit = SQRT(tcrit)                      !critical time increment
 IF( tcrit < deltc ) deltc = tcrit

 RETURN
 END SUBROUTINE delt25
