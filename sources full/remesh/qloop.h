SUBROUTINE qloop(nelbnd,cdbnd,nlbnd)
! delete nodes in the loop of contour lines
IMPLICIT NONE

  !--- Dummy variables
  INTEGER(kind=4),INTENT(INOUT):: nelbnd  !number of segments defining contour
  REAL(kind=8),POINTER:: cdbnd(:,:)       !coordinates of segments defining contour
  INTEGER(kind=4),POINTER:: nlbnd(:)      !segments in each contour lines

END SUBROUTINE  qloop
