 SUBROUTINE nearst(lcseg,nsegm,nearn,xs,x)

 !.... look for the nearest master segment using a GLOBAL searching algorithm
 !
 !.... input
 !....   lcseg(2,nsegm) connectivities
 !....   x  = current spatial coordinates
 !....   xs = current spatial coordinates of the slave node
 !.... output
 !....   nearn = number of the nearest master segment

 IMPLICIT NONE
 !     arguments
 INTEGER (kind=4), INTENT (IN) :: nsegm,lcseg(:,:)
 INTEGER (kind=4), INTENT (IN OUT) :: nearn
 REAL (kind=8), INTENT (IN) :: x(:,:),xs(:)
 !     local variables
 INTEGER (kind=4) iseg
 REAL    (kind=8) vdist(3),d,dmin
 LOGICAL :: first,down

 vdist = xs - (x(:,lcseg(1,nearn)) + x(:,lcseg(2,nearn)))/2d0
 dmin  = DOT_PRODUCT(vdist,vdist)  !present distance
 first = .TRUE.                    !first time
 down  = .TRUE.
 !.... loop over all master segments
 iseg = nearn

 DO
   IF( down )THEN
     IF( iseg < nsegm )THEN         !look forward
       iseg = iseg + 1              !increase segment
       vdist = xs - (x(:,lcseg(1,iseg)) + x(:,lcseg(2,iseg)))/2d0 !distance vector
       d = DOT_PRODUCT(vdist,vdist)    ! distance (squared)
       IF(d < dmin) THEN               ! check for minimum distance
         dmin = d                      ! updates minimun distance
         nearn = iseg                  ! updates nearest node
         first = .FALSE.               ! search updated forward
       ELSE
         down = .FALSE.
         iseg = iseg - 1
       END IF
     END IF
   ELSE
     IF( first .AND. iseg > 1)THEN
       iseg = iseg - 1
       vdist = xs - (x(:,lcseg(1,iseg)) + x(:,lcseg(2,iseg)))/2d0
       d = DOT_PRODUCT(vdist,vdist)    ! distance (squared)
       IF(d < dmin) THEN               ! check for minimum distance
         dmin = d                      ! updates minimun distance
         nearn = iseg                  ! updates nearest node
       ELSE
         EXIT
       END IF
     ELSE
       EXIT
     END IF
   END IF
 END DO
 RETURN
 END SUBROUTINE nearst
