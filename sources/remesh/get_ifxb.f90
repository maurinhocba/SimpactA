SUBROUTINE get_ifxb(nelbnd,lnbnd,ifx_b)
! get ifx codes and coordinates for boundary nodes in a remeshed set
USE ctrl_db,ONLY: ndofn
USE ifx_db,ONLY: nd1, ifx_nod, ihead, itail, srch_ifx
USE npo_db,ONLY: label
IMPLICIT NONE

  !--- Dummy variables (INPUT)
  INTEGER(kind=4):: nelbnd               !number of segments defining contour
  INTEGER(kind=4),POINTER:: lnbnd(:)     !Conectivities of boundary segments
  !--- Dummy variables (OUTPUT)
  INTEGER(kind=4),POINTER:: ifx_b(:,:)   !(nd1,nelem)   boundary conditions
  !--- Local variables
  LOGICAL:: found
  INTEGER(kind=4):: i, n
  TYPE(ifx_nod),POINTER:: anter, posic

  !=============================

  ALLOCATE(ifx_b(ndofn+1,nelbnd))  !( NDOF+1(or not) , number of segments )
  ifx_b(1:ndofn+1,1:nelbnd) = 0    !Initializes
  ! loop over segments
  DO i=1,nelbnd                ! this seems to work for 2-node segments only (2-D problems only)
    n = label(lnbnd(i))      ! label of the second node
    CALL srch_ifx(ihead, anter, posic, n, found)  !search in the list of fixities
    IF (found) ifx_b(1:nd1,i)=posic%ifix(2:nd1+1) !assign nodal boundary cond.
  END DO

RETURN
END SUBROUTINE get_ifxb
