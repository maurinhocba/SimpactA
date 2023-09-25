SUBROUTINE pregid(scalef,nnode,nelbnd,cdbnd,nlbnd,ifx_b)
! preparing a GiD batch file
!USE meshmo_db,ONLY: pline, headp, ltype, lstrt, nplins, r_elm_size
!USE c_input,ONLY: openfi
!USE param_db,ONLY: mlen, mnam, mstl
!USE name_db,ONLY: output
!USE dfport
IMPLICIT NONE

  !--- Dummy variables
  INTEGER(kind=4),INTENT(IN):: nnode, nelbnd
  REAL(kind=8),INTENT(IN):: scalef
  REAL(kind=8),POINTER:: cdbnd(:,:)
  INTEGER(kind=4),POINTER:: nlbnd(:)
  INTEGER(kind=4),POINTER:: ifx_b(:,:)

END SUBROUTINE pregid
