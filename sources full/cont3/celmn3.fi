 SUBROUTINE celmn3(ttime,dtime,td,disma,coord,emass,fcont,coorb,coort)

 !.... perform contact searching & computes contact forces

 IMPLICIT NONE
 !     Dummy arguments
 REAL (kind=8), INTENT(IN) :: ttime,dtime,td,disma,emass(:,:)
 REAL (kind=8), POINTER :: coord(:,:),coorb(:,:),coort(:,:) !INTENT(IN)
 REAL (kind=8), INTENT(IN OUT) :: fcont(:,:)
 !     local variables
 INTEGER (kind=4) ipair,ims,iss,i
 INTEGER (kind=4), SAVE :: istep = 0
 REAL (kind=8), SAVE :: ffdis = 1d0
 REAL (kind=8), POINTER :: coors(:,:),coorm(:,:)
 TYPE (pair3_db), POINTER :: pair
 TYPE (surf3_db), POINTER :: master,slave

 IF(oldis < disma)THEN               !modify factor ffdis to account
   ffdis = oldis/disma               !for velocity changes
   IF(ffdis > 0.975) ffdis = 1d0     !minor changes => set to 1
 ELSE
   ffdis = 1d0                       !no changes
 END IF
 oldis = disma           !keep actual maximum displacement for comparison

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
       IF( .NOT.pair%prev )THEN
         pair%prev  = .TRUE.    !set Previous Step as TRUE
         master%cxc = .TRUE.    !set to recompute coordinates
       END IF
       !If necessary set to recompute coordinates of triangles centers
       IF( MOD(istep,pair%freq) == 0 )master%cxc = .TRUE.

       SELECT CASE (pair%mtsur) !select coordinates to use for master surface
       CASE(-2)
         coorm => coorb               !bottom
       CASE(-1,1)
         coorm => coord               !central
       CASE( 2)
         coorm => coort               !top
       END SELECT

       IF( pair%auto )THEN
         CALL csrfc3a(ipair,coorm,oldis,ffdis,pair,master, &
                      dtime,emass,fcont,.TRUE.,td)

       ELSE
         iss = pair%islav    !slave surface order in list
         slave => shead      !point to first surface
         DO i=1,iss-1        !loop until pointer is correct
           slave => slave%next
         END DO
         IF( slave%ncnod > 0 )THEN
           SELECT CASE (pair%slsur) !select coordinates to use for slave surface
           CASE(-2)
             coors => coorb               !bottom
           CASE(-1,1)
             coors => coord               !central
           CASE( 2)
             coors => coort               !top
           END SELECT

           !.... go into next routine
           CALL csrfc3(ipair,coors,coorm,oldis,ffdis,pair,master,slave, &
                       dtime,emass,fcont,.TRUE.,td)
         END IF
       END IF
     END IF
   END IF
   pair => pair%next     !point to next pair
 END DO
 istep = istep + 1                   !increments counter
 RETURN
 END SUBROUTINE celmn3
