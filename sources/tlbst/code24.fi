 SUBROUTINE code24(nelem ,heade ,deltc , coord )

 !     Compute Critical time increment for a set of elements BST-CST (TLF)
 !     for thermal analysis
 USE ctrl_db, ONLY : tmscal
 IMPLICIT NONE
 INTEGER (kind=4), INTENT (IN) :: nelem  !number of elements in the set
 REAL (kind=8), INTENT(IN) :: coord(:,:) !present mesh coordinates
 REAL (kind=8), INTENT(IN OUT) :: deltc  !critical time increment
 TYPE (ele24), POINTER :: heade          !pointer to first element

 INTEGER :: ielem,in,jn,kn,isec,osec,secty
 REAL (kind=8) :: t1(3),densi,condu,capac,fac, &
                  a1,a2,a3,a4,ma,thick,tcrit
 TYPE (ele24), POINTER :: elem
 TYPE (section), POINTER :: sec

 osec = -1                             !previous material number
 elem => heade                         !point to first element
 tcrit = deltc                         !initializes
 DO ielem=1,nelem                      !for each element in the set
   isec = elem%matno                  !associated material
   IF( isec /= osec )THEN             !present and previous material diff
     sec => psecs(isec)%p
     secty  = sec%secty
     !        this must be modified
     SELECT CASE (ndoft)
     CASE (1)
       thick = 1
     CASE (2)
       thick  = sec%rprop(1)            !thickness
     CASE (3)
       thick  = sec%rprop(1)/2.5d0      !thickness
     END SELECT
     thick = thick*thick
     ! properties may be temperature dependant
     densi  = sec%mtbas%prope(5)*tmscal   !mass density (scaled)
     capac  = sec%mtbas%prope(7)      !heat capacity
     condu  = sec%mtbas%prope(8)      !heat conductivity
     fac    = capac*densi/condu       !factor
     osec  = isec                     !update previous material
   END IF

   in = elem%lnods(1)                  !first node of the element
   jn = elem%lnods(2)                  !second node of the element
   kn = elem%lnods(3)                  !third node of the element
   t1 = coord(1:3,jn) - coord(1:3,in)  !first side 1-2
   a1 = DOT_PRODUCT(t1,t1)             !first side length (squared)
   t1 = coord(1:3,kn) - coord(1:3,jn)  !second side 2-3
   a2 = DOT_PRODUCT(t1,t1)             !first side length (squared)
   t1 = coord(1:3,kn) - coord(1:3,in)  !third side 1-3
   a3 = DOT_PRODUCT(t1,t1)             !third side length (squared)
   a4 = thick*elem%lb*elem%lb
   ma = MIN(a1,a2,a3,a4)*fac
   tcrit = MIN(tcrit,ma)
   elem => elem%next                         !point to next element
 END DO

 IF( tcrit < deltc ) deltc = tcrit

 RETURN
 END SUBROUTINE code24
