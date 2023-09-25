 SUBROUTINE projt1(xs,nearn,x,lcseg,rssdb,   &
                   issdb,prdat,dd,nsegm) !,curms,cu)

 !.... project the slave node onto a 2-d master segment

 IMPLICIT NONE
 !     arguments
 !LOGICAL, INTENT(IN) :: curms    !curved master surface?
 INTEGER (kind=4), INTENT(IN) :: nsegm,       & !number of segments in master surface
                                 lcseg(:,:)     !connectivities of master surface

 INTEGER (kind=4), INTENT(IN OUT) :: nearn,   & !nearest master segment to node (candidate)
                                     issdb(:)   !slave node Integer data base
 REAL (kind=8), INTENT (IN) :: x(:,:),        & !actual nodal coordinates
                               dd,            & !radius difference
                               xs(:)            !slave node position
                               !cu(:)            !curvatures of master surface
 REAL (kind=8), INTENT (IN OUT) :: rssdb(:)     !slave node Real data base
 REAL (kind=8), INTENT (OUT) :: prdat(:)        !penetration data (10)

 !     local variables
 INTEGER (kind=4) :: imn,jmn     !first and second node of master segment
 REAL (kind=8) :: vd(3),       & !vector from slave node to first master node
                  vj(3),       & !segment vector from first to second master node
                  vt(3),       & !proyection vector over master segment
                  lengt,       & !segment length
                  gap,         & !penetration gap
                  xita,proj,   & !local coordinate in master segment [0,1]
                  vn(3)          !segment normal and average nodal normal

 LOGICAL ::  sig                !effective projected penetration ?
 REAL (kind=8), PARAMETER :: toler = -0.01d0, tole1 = 1.01d0  !tolerance for end projections

! REAL (kind=8), SAVE :: maxg = -1d-6  !variable for debug

 !.... initialize to zero some control parameters
 sig    = .FALSE.                   !initializes to NO full projection

 !....  First SEARCH for projection
 search : DO              !.... loop until projection
   imn = lcseg(1,nearn)               !1st master node
   jmn = lcseg(2,nearn)               !2nd master node
   vd = xs(:) - x(:,imn)              !distance vector from first master node to slave node
   vt = x(:,jmn) - x(:,imn)           !form tangent vector to master segment, side i-->j
   lengt = SQRT(DOT_PRODUCT(vt,vt))   !segment length
   vt = vt/lengt                      !unit length side
   proj = DOT_PRODUCT(vd,vt)          !project distance vector over the tangent (proyection)
   xita = proj/lengt                  !project distance vector over the tangent (proyection)
   vj = vd - proj*vt                  !distance vector from node to axis
   gap = DOT_PRODUCT(vj,vj)           !distance squared
   ! according to proyection
   IF( xita < 0d0 )THEN               ! must go to previous segment
     IF( xita > toler )THEN           ! but near first node
       IF( sig )THEN                  ! if there is a previous proyection
         vj = ( vn + vj )/2d0           ! average distance vectors
         gap = DOT_PRODUCT(vn,vn)       ! distance squared
         xita = 0d0                     ! send to vertex
         EXIT search
       ELSE                           ! else keep data
         vn = vj
         sig = .TRUE.
       END IF
     ELSE IF(sig) THEN
       nearn = nearn - 1
       vj = vn
       gap = DOT_PRODUCT(vj,vj)           !distance squared
       xita = 1d0
       EXIT search
     END IF
     nearn = nearn - 1                ! point to previous segment
     IF( nearn == 0 )THEN
       nearn = 1
       IF( sig )  xita = 0d0
       EXIT search
     END IF
   ELSE IF( xita > 1d0 )THEN          ! must go to next segment
     IF( xita < tole1 )THEN           ! if the node is near second node
       IF( sig )THEN                  ! if there is a previous proyection
         vj = ( vn + vj )/2d0           ! average distance vector
         gap = DOT_PRODUCT(vn,vn)       ! distance squared
         xita = 1d0                     ! send to vertex
         EXIT search                    ! exit
       ELSE                           ! else keep data
         vn = vj
         sig = .TRUE.
       END IF
     ELSE IF(sig) THEN
       nearn = nearn + 1
       vj = vn
       gap = DOT_PRODUCT(vj,vj)           !distance squared
       xita = 0d0
       EXIT search
     END IF
     nearn = nearn + 1                ! point to next segment
     IF( nearn > nsegm ) THEN
       nearn = nsegm
       IF( sig )  xita = 1d0
     EXIT search                      !exit search
     END IF
   ELSE
     sig = .TRUE.
     EXIT search                      !exit search
   END IF
 END DO search   !projection loop seacrh
 !.... Second SELECT projection type
 ! keep vn, gap
 IF( sig )THEN
   gap = SQRT(gap)   !distance between axis
   vj = vj/gap       !unit vector normal
   gap =  dd - gap   !penetration
   prdat(1) = gap       !store penetration positive or negative
   IF( gap < 0d0 )THEN  !effective penetratin
     !modify gap and normal according to surface curvature
     !IF( curms )THEN                    !if master surface is defined as curved
     !  gap = gap - xita*(1d0-xita) *cu(nearn)/2d0   !corrects penetration
     !  bx = (0.5d0-xita)*cu(nearn)/2d0               !second derivative of position
     !  vj(1) = -(vj(2)*lengt+bx*vn(2))                !corrects normal (x component)
     !  vj(2) =  (vj(1)*lengt+bx*vn(1))                !                (y component)
     !  vj = vj/SQRT(DOT_PRODUCT(vj,vj))               !unit normal vector
     !END IF                             !close curved master surface corrections
     !continue with check if penetration is in the allowed range
     prdat(2)   = xita       !store local coordinate
     prdat(3:5) = -vj        !distance vector
     prdat(6:8) = vt         !tangent vector
     prdat(9) = lengt        !segment length
   ELSE                 !no contact
     issdb(2:3) = 0       !forget previous proyection
   END IF
   issdb(1) = nearn     !actual segment
   rssdb(2) = gap       !keep previous penetration ????
 END IF

 RETURN
 9999 STOP
 END SUBROUTINE projt1
