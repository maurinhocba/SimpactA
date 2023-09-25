 SUBROUTINE nears3(xc,nsegm,nearn,xs)

!.... look for the nearest master segment using a GLOBAL searching algorithm
!
!.... input
!....   xc = current spatial coordinates of the segment center
!....   xs = current spatial coordinates of the slave node
!.... output
!....   nearn = number of the nearest master segment

 IMPLICIT NONE
!     arguments
 INTEGER (kind=4), INTENT (IN) :: nsegm
 INTEGER (kind=4), INTENT (OUT) :: nearn
 REAL (kind=8), INTENT (IN) :: xc(:,:),xs(:)
!     local variables
 INTEGER (kind=4) iseg
 REAL    (kind=8) vdist(3),d,dmin

!.... initialize values
 dmin  = 1.0e+10
 nearn = 0
!.... loop over all master segments
 DO iseg = 1, nsegm
   vdist = xs(1:3) - xc(1:3,iseg)  !  distance vector  xs - xc
   d = DOT_PRODUCT(vdist,vdist)    ! distance (squared)
   IF(d < dmin) THEN               ! check for minimum distance
     dmin = d                      ! updates minimun distance
     nearn = iseg                  ! updates nearest node
   END IF
 END DO
 RETURN
 END SUBROUTINE nears3
