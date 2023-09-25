 SUBROUTINE tdel19(nelem,heade,deltc,coord)

 !     Compute Critical time increment for a set of elements linear triangle (TLF)
 !     for thermal analysis
 USE ctrl_db, ONLY : tmscal
 IMPLICIT NONE
 INTEGER (kind=4), INTENT (IN) :: nelem  !number of elements in the set
 REAL (kind=8), INTENT(IN) :: coord(:,:) !present mesh coordinates
 REAL (kind=8), INTENT(IN OUT) :: deltc  !critical time increment
 TYPE (ele19), POINTER :: heade          !pointer to first element

 INTEGER :: i,ielem,in,jn,kn,isec,osec,secty,numct, & 
            ncurv(2)  !associated curves for therm dependant properties
 REAL (kind=8) :: t1(2),densi,condu,capac,fac, &
                  a1,a2,a3,crite,tcrit,val
 TYPE (ele19), POINTER :: elem
 TYPE (section), POINTER :: sec
 TYPE (mater), POINTER :: mat    !pointer to a section data 
 TYPE (curve), POINTER :: crvc,crvk !pointer to curves

 osec = -1                             !previous material number
 elem => heade                         !point to first element
 tcrit = deltc                         !initializes
 DO ielem=1,nelem                      !for each element in the set
   isec = elem%matno                  !associated material
   IF( isec /= osec )THEN             !present and previous material diff
     sec => psecs(isec)%p
     secty  = sec%secty
     mat => sec%mtbas     
     densi  = mat%prope(5)*tmscal   !mass density (scaled)
     capac  = mat%prope(7)          !heat capacity
     condu  = mat%prope(8)          !heat conductivity
     fac    = capac*densi/condu       !thermal diffusivity
     !associated curves for thermal dependant properties
     ncurv(1) = INT(mat%props(7))     !heat capacity
     ncurv(2) = INT(mat%props(8))     !heat conductivity
     CALL cur_point(ncurv(1),mat%chead,crvc) !pointer to young curve
     CALL cur_point(ncurv(2),mat%chead,crvk) !pointer to poisson curve
     osec  = isec                     !update previous material     
   END IF
   numct = SUM(ncurv)
   IF( numct > 0 )THEN !compute properties as a function of temperature
     val = SUM(tempe(1,elem%lnods(1:3)))/3d0 ! temperature at element center
     i = 1
     IF(ncurv(1) > 0) capac = inte_cr(crvc%val,crvc%np,val,i)
     IF(ncurv(2) > 0) condu = inte_cr(crvk%val,crvk%np,val,i)
     fac    = capac*densi/condu       !evaluates thermal diffusivity 
   END IF

   in = elem%lnods(1)                  !first node of the element
   jn = elem%lnods(2)                  !second node of the element
   kn = elem%lnods(3)                  !third node of the element   
   t1 = coord(1:2,jn) - coord(1:2,in)  !first side 1-2
   a1 = DOT_PRODUCT(t1,t1)             !first side length (squared)
   t1 = coord(1:2,kn) - coord(1:2,jn)  !second side 2-3
   a2 = DOT_PRODUCT(t1,t1)             !first side length (squared)
   t1 = coord(1:2,kn) - coord(1:2,in)  !third side 1-3
   a3 = DOT_PRODUCT(t1,t1)             !third side length (squared)
   crite = MIN(a1,a2,a3)*fac  !element thermal critical time step
   tcrit = MIN(tcrit,crite)   !minimum thermal critical time step

   elem => elem%next                         !point to next element
 END DO

 IF( tcrit < deltc ) deltc = tcrit !

 RETURN
 END SUBROUTINE tdel19
