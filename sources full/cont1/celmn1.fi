 SUBROUTINE celmn1(ttime,dtime,coord,emass,fcont,velnp)

 !.... perform contact searching & computes contact forces

 !USE cont1_db    !INTENT(IN OUT)
 IMPLICIT NONE
 !     Dummy arguments
 REAL (kind=8), INTENT(IN) :: ttime,dtime,emass(:,:),velnp(:,:)
 REAL (kind=8), POINTER :: coord(:,:)           !INTENT(IN)
 REAL (kind=8), INTENT(IN OUT) :: fcont(:,:)
 !     local variables
 INTEGER (kind=4) ipair,ims,iss,i
 TYPE (pair1_db), POINTER :: pair
 TYPE (surf1_db), POINTER :: master,slave


 pair => headp           !pointer to first pair
 DO ipair = 1, npair     !for each pair
   ! if pair is active for this time
   IF( pair%start <= ttime .AND. pair%end >= ttime)THEN
     !.... identify master and slave surfaces numbers
     ims = pair%imast    !master surface order in list
     master => shead     !point to first surface
     DO i=1,ims-1        !loop until pointer is correct
       master => master%next
     END DO
     IF( master%nsegm > 0 )THEN
       IF( .NOT.pair%prev )  pair%prev  = .TRUE.    !set Previous Step as TRUE

       iss = pair%islav    !slave surface order in list
       slave => shead      !point to first surface
       DO i=1,iss-1        !loop until pointer is correct
         slave => slave%next
       END DO
       IF( slave%ncnod > 0 )THEN
         !.... go into next routine
         CALL csrfc1(ipair,coord,pair,master,slave, &
                     dtime,emass,fcont,velnp,.TRUE.)
       END IF
     END IF
   END IF
   pair => pair%next     !point to next pair
 END DO
 RETURN
 END SUBROUTINE celmn1
