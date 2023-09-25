 SUBROUTINE nearsa(xc,nsegm,nearn,xs,tc,tn,maxdi,cutoff)

 !.... look for the nearest master segment using a GLOBAL searching algorithm
 !.... output:  nearn = number of the nearest master segment
 USE ctrl_db, ONLY : ndime
 IMPLICIT NONE
 !     arguments
 INTEGER (kind=4), INTENT (IN) :: nsegm   !number of segments in surface
 INTEGER (kind=4), INTENT (OUT) :: nearn  !segment candidate
 REAL (kind=8), INTENT (IN) :: xc(:,:), & ! xc = current coordinates of the segment center
                               xs(:),   & ! xs = current coordinates of the slave node
                               tc(:,:), & ! tc = current normals at segments
                               tn(:),   & ! tn = current normal at node
                               maxdi,   & ! maximum distance
                               cutoff     ! maximum penetration
 !     local variables
 INTEGER (kind=4) iseg
 REAL    (kind=8) vdist(ndime),d,dmin,cang

 !.... initialize values
 dmin  = 1d10        !a large value
 nearn = 0           !no candidate
 !.... loop over all master segments
 DO iseg = 1, nsegm
   cang = DOT_PRODUCT(tc(:,iseg),tn)     !angle between normals
   IF( cang > 0d0 )CYCLE                 !if |angle| < pi/2 discard segment
   vdist = xs(:) - xc(:,iseg)            !distance vector  xs - xc
   cang = DOT_PRODUCT(tc(:,iseg),vdist)  !penetration
   IF( cang < cutoff .OR. cang > maxdi )CYCLE !penetration is too large or point too distant
   d = DOT_PRODUCT(vdist,vdist)    ! distance (squared)
   IF(d < dmin) THEN               ! check for minimum distance
     dmin = d                      ! updates minimun distance
     nearn = iseg                  ! updates nearest segment
   END IF
 END DO
 RETURN
 END SUBROUTINE nearsa
