SUBROUTINE plines(cdbnd,nlbnd,fside,ifx_b)
! Analyze contour fixity codes to find polylines
USE ifx_db,ONLY: nd1       !number of fixity codes per node
USE meshmo_db,ONLY: pline, headp, tailp, nstart, nplins, new_pline, add_pline
IMPLICIT NONE

  !--- Dummy variables
  INTEGER(kind=4):: ifx_b(:,:)   !(nd1,nelbnd) fixities at second nodes
  REAL(kind=8):: cdbnd(:,:)      !(ndime,nelbnd) coords at second nodes
  INTEGER(kind=4):: nlbnd(:)     !(nline) number of segments in each boundary line
  INTEGER(kind=4):: fside(:)     !(nelbnd) free-side segments

END SUBROUTINE plines
