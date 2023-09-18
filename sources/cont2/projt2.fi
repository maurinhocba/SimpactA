 SUBROUTINE projt2(xs,nearn,x,lcseg,rssdb,   &
                   nhseg,issdb,prdat,disma,cutof,gapin,curms,cu)

 !.... project the slave node onto a 2-d master segment

 IMPLICIT NONE
 !     arguments
 LOGICAL, INTENT(IN) :: curms    !curved master surface?
 INTEGER (kind=4), INTENT(IN) :: lcseg(:,:),  & !connectivities of master surface
                                 nhseg(:,:)     !NeigHbour segments for each master segment
 INTEGER (kind=4), INTENT(IN OUT) :: nearn,   & !nearest master segment to node (candidate)
                                     issdb(:)   !slave node Integer data base
 REAL (kind=8), INTENT (IN) :: x(:,:),        & !actual nodal coordinates
                               disma,         & !maximum incremental displacement
                               cutof,         & !maximum allowed penetration
                               gapin,         & !maximum allowed incremental gap
                               xs(:),         & !slave node position
                               cu(:)            !curvatures of master surface
 REAL (kind=8), INTENT (IN OUT) :: rssdb(:)     !slave node Real data base
 REAL (kind=8), INTENT (OUT) :: prdat(:)        !penetration data

 !     local variables
 INTEGER (kind=4) :: imn,jmn,  & !first and second node of master segment
                     nchck,    & !number of segments with effective penetration
                     masts,    & !final master segment
                     loop,     & !number of candidates tested
                     nnear,    & !new element tested
                     i,        & !counter
                     o(0:2)      !tested segments
 REAL (kind=8) :: vd(2),       & !vector from slave node to first master node
                  vj(2),       & !segment vector from first to second master node
                  vt(2),       & !proyection vector over master segment
                  lengt,       & !segment length
                  gap,         & !penetration gap
                  dist,        & !distance from slave node to first master node
                  y(2),        & !equivalent position in averaged penetrations
                  mgap,        & !maximum gap found
                  xita,xita1,  & !local coordinate in master segment [0,1]
                  vn(2),vnn(2),& !segment normal and average nodal normal
                  vns(2),gaps, & !candidate segment normal and gap
                  bx,          & !curvature correction
                  pgap,        & !previous step gap
                  gapin0,      & !maximum incremental gap
                  cutof0         !maximum penetration
 LOGICAL ::  chck1,            & !projects over segment?
             chck2,            & !gap < 0? (penetration)
             sig                 !effective projected penetration ?
 REAL (kind=8), PARAMETER :: toler = 0.01d0, tole1 = 1.01d0  !tolerance for end projections

! REAL (kind=8), SAVE :: maxg = -1d-6  !variable for debug

 !.... initialize to zero some control parameters
 pgap = MIN( rssdb(2), 0d0)         !previos penetration
 gapin0 = MIN(gapin,pgap/2d1)       !maximum incremental gap
 cutof0 = cutof + pgap/2d0          !update maximum penetration allowed
 cutof0 = MAX(cutof0, 3d0*cutof )   !maximum penetration for this node
 nchck  = 0                         !initializes number of candidates to 0
 sig    = .FALSE.                   !initializes to NO full projection
 mgap   = 1d10                      !initializes maximum gap found to a huge number
 vnn    = 0d0                       !initializes weighted normal
 masts  = 0                         !initializes definite master segment to Nothing
 loop   = 0                         !initializes loop search
 o(loop) = nearn                    !store first candidate

 !....  First SEARCH for projection
 search : DO              !.... loop until projection
   imn = lcseg(1,nearn)               !1st master node
   jmn = lcseg(2,nearn)               !2nd master node
   vd = xs(1:2) - x(1:2,imn)          !distance vector from first master node to slave node
   dist = SQRT(DOT_PRODUCT(vd,vd))    !distance
   mgap = MIN(mgap,dist)              !for boundaries use distane as minimum gap
   IF(dist == 0d0) THEN               !IF nodes coincide
     nchck = 0                        !no penetration
     sig = .TRUE.                     !keep the master node and segment
     mgap = 0d0                       !update max gap
     masts = nearn                    !master segment is present candidate
     EXIT search                      !exit search
   END IF                             !close null distance check
   vj = x(1:2,jmn) - x(1:2,imn)       !form tangent vector to master segment, side i-->j
   lengt = SQRT(DOT_PRODUCT(vj,vj))   !segment length
   vj = vj/lengt                      !unit length side
   vn(1) = -vj(2)                     !normal vector (unit length) (x component)
   vn(2) =  vj(1)                     !                            (y component)
   gap = DOT_PRODUCT(vd,vn)           !project distance vector over the normal ==> gap
   vt = vd - gap*vn                   !the projection of the distance vector over tangent line
   xita1 = DOT_PRODUCT(vt,vj)/lengt   !local coordinate [0,1] means projects in
   chck1 = xita1 >= 0d0 .AND. xita1 <= 1d0 !check IF the projection is over the segment
   !modify gap and normal according to surface curvature
   IF( curms )THEN                    !if master surface is defined as curved
     gap = gap - xita1*(1d0-xita1) *cu(nearn)/2d0   !corrects penetration
     bx = (0.5d0-xita1)*cu(nearn)/2d0               !second derivative of position
     vn(1) = -(vj(2)*lengt+bx*vn(2))                !corrects normal (x component)
     vn(2) =  (vj(1)*lengt+bx*vn(1))                !                (y component)
     vn = vn/SQRT(DOT_PRODUCT(vn,vn))               !unit normal vector
   END IF                             !close curved master surface corrections
   !discard segment if GAP is too large
   IF( gap < cutof0 ) THEN            !compare actual penetration with maximum allowed
     !IF(chck1 .AND. gap/cutof0 < 3d0) &
     !  WRITE(58,"('ct',4i6,2e13.4)",ERR=9999)istep,isn,lcseg(:,nearn),gap,rssdb(2)  !debug print
     IF( loop == 2 )EXIT search    ! too many segments considered (3 is the limit in 2-D)
     !select a new candidate according to projection
     IF( xita1 <= 0.5d0 )THEN      !if local coordinate is less than 1/2
       nnear = nhseg(1,nearn)      !previous segment
     ELSE                          !if local coordinate is greater than 1/2
       nnear = nhseg(2,nearn)      !next segment
     END IF                        !close segment selection
     IF( nnear ==  0 )EXIT search  !No side element, EXIT
     DO i=0,loop                      !see if element not checked yet
       IF( nnear == o(i) )EXIT search !element already checked
     END DO                           !
     loop = loop + 1      !increase number of candidates considered
     IF( loop == 3 )EXIT search  !too many checks, EXIT search
     nearn = nnear        !new element to check
     o(loop) = nnear      !store element as checked
     CYCLE search         !new search
   END IF           !close check of maximum allowed penetration
   !continue with check if penetration is in the allowed range
   mgap = MIN(gap,mgap)             !Compare with previous gap
   chck2 = (gap < 0d0 )             !gap < 0  == effective penetration
   IF( chck1 )THEN                        !IF projects over the segment
     sig = .TRUE.                         !keep full projection
     masts = nearn                        !adopt segment candidate
     xita  = xita1                        !keep local coordinate
     IF ( chck2 )THEN                         !IF Penetration
       gaps  = gap                            !keep gap
       nchck = 1                              !number of efective penetrations
     END IF                                   !close penetration check
   ELSE                                   !projection is out of the element boundary
     IF( xita1 < 0d0 )THEN                    !if projects over previous segment
       nnear = nhseg(1,nearn)                 !previous segment
     ELSE   !(xita1 > 1d0)                    !if projects may be on next segment
       nnear = nhseg(2,nearn)                 !next segment
     END IF                                   !close selection
     IF ( chck2 ) THEN                        !if penetrated make another check
       IF( xita1 > -toler .AND. xita1 < tole1 )THEN !check if projects over the expanded segment
         ! next two lines could be commented)
         IF(xita1 < 0d0 ) xita1 = 0d0   !compute corrected coordinates
         IF(xita1 > 1d0 ) xita1 = 1d0
         nchck = nchck + 1   !increase number of projected elements
         IF(nchck == 1 .OR. gap > gaps) THEN  !if the first or nearer save as segment (ims) for projection
           masts = nearn      ! adopt segment candidate
           vns   = vn         ! adopt segment normal
           gaps  = gap        ! adopt gap candidate
           xita = xita1       ! keep local coordinate
         END IF  ! nchck = 1 or gap > gaps
         vnn = vnn - vn/gap   !weigthed normal (gap inverse weighting)
       END IF  ! projection over expanded segment
     END IF ! if gap is negative (effective penetration)
     IF( nnear > 0 )THEN       !if a new segment to search exists
       DO i=0,loop             !see if it has not been checked yet
         IF( nnear == o(i) ) EXIT search  !if checked EXIT search
       END DO
       loop = loop + 1         !not checked, new loop
       IF( loop == 3 )EXIT search  !too many check, EXIT search
       nearn = nnear           !new segment to check
       o(loop) = nnear         !store element as checked
       CYCLE search            !go for a new check
     END IF                    !close check if a new candidate exists
   END IF !close selection acording to segment projection
   EXIT  search     !segment OK, or no side element EXIT search
 END DO search   !projection loop seacrh

 !.... Second SELECT projection type

 IF(nchck >= 1) THEN         !if at least a candidate exists
   IF( issdb(3) > 0 )issdb(2) = issdb(1)   !store previous segment
   issdb(1) = nearn          !actual segment
   IF(nchck == 1) THEN       !CASE 1: project onto one master segment
     gap = gaps              !gap
     prdat(2:3) = vn         !normal
   ELSE                      !CASE 2: weigthed projection onto node
     gap = SQRT(DOT_PRODUCT(vnn,vnn))  !auxiliar value
     vnn = vnn/gap                     !unit weighted normal
     gap = gaps/ DOT_PRODUCT(vnn,vns)  !definite gap
     y = xs(1:2) - x(1:2,imn) - vnn * gap  !equivalent position over segment
     xita = DOT_PRODUCT(y,vj)/lengt    !recompute local coordinate
     prdat(2:3) = vnn                  !store normal
!     !WRITE(58,"(e12.4,3f9.6)",ERR=9999)gap,vnn,xita  !debug print
   END IF                    !end selection (number of nchck)
   IF(gap - pgap < gapin0)THEN  !check that incremental gap is less than the allowed
     !WRITE(58,"('mg',4i6,2e13.4)",ERR=9999)istep,isn,lcseg(:,nearn),gap,rssdb(2) !debug print
     gap = pgap + gapin0        !corrects gap to comply with maximum incremental gap
   END IF                       !end check
   prdat(1) = gap       !store penetration
   prdat(4) = xita      !store local coordinate
   rssdb(2) = gap       !keep previous penetration
!   IF( gap < maxg )THEN                     ! compare maximum gap for debug purpose
!     WRITE(58,"(' mgap ',2e15.3)",ERR=9999)gap,maxg  ! debug print
!     maxg = gap                             ! update maximum gap (debug)
!   END IF                                   ! close check        (debug)
 ELSE                       !no candidate found, modify data base for future checks
   issdb(3) = issdb(3) - 1  !increase number of steps without penetration
   ! Next line have sense if no Estrict penetration is enough for friction
   IF( -issdb(3) > 0 ) issdb(2) = 0  !forget previous proyection
   prdat(1) = 0d0           ! gap >= 0
   ! keep node candidate for future search
   IF(mgap > disma)THEN     !   **** IF distance is large enough
     issdb(1) = - nearn - 10000000      !no search for 10 steps
   ELSE IF(sig)THEN         !if an appropiate projection exists
     issdb(1) =  nearn      !keep segment with projection
   ELSE                     !if no projection was found
     issdb(1) =  o(0)       !keep original nearest segment
   END IF                   !end selection according to mgap
   IF(mgap > 0d0) THEN      !if a positive distance was computed
      rssdb(2) = mgap       !keep minimum distance
   ELSE                     !if a negative distance was computed
      rssdb(2) = MAX(0d0,rssdb(2))  !keep minimum positive distance
   END IF                   !close selection
 END IF                     !selection according to number of candidates found

 RETURN
 9999 STOP
 END SUBROUTINE projt2
