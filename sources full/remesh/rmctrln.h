SUBROUTINE rmctrln(nelem,heade,nelbnd,cdbnd,lnbnd,nlbnd,fside)
!-----------------------------------------------------------------------------------
!  Create a conectivities and coordinates array for contour line (closed continuous line)
!  Routine works for 2-node segments only (2D problem)
!-----------------------------------------------------------------------------------
USE npo_db,ONLY: coora
USE surf_db,ONLY: srf_seg
IMPLICIT NONE

  !--- Dummy variables
  INTEGER(kind=4),INTENT(IN):: nelem     !number of elements
  INTEGER(kind=4),INTENT(OUT):: nelbnd   !Number of contour segments
  INTEGER(kind=4),POINTER:: lnbnd(:)     !Conectivities of nodes segments
  INTEGER(kind=4),POINTER:: nlbnd(:)     !Number of continuos lines
  INTEGER(kind=4),POINTER:: fside(:)     !Free side array
  REAL(kind=8),POINTER:: cdbnd(:,:)      !Coordinates of nodes segments
  TYPE(srf_seg),POINTER:: heade

END SUBROUTINE rmctrln
