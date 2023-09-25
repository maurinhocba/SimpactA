SUBROUTINE dispcrv(npoin,nelem,lnods,lside,ldv2sd,nodcre,nwlnd,lindi,cd,cvdisp)
IMPLICIT NONE

  !--- Dummy variables
  INTEGER(kind=4),INTENT(IN):: npoin,          & !Number of nodes 
                               nelem             !Number of elements in the set
  INTEGER(kind=4),POINTER:: lnods(:,:),     & !Conectivities of the set
                            ldv2sd(:,:),    & !Element side subdivision
                            lside(:,:),     & !Neighbor elements of the set
                            nodcre(:,:),    & !Nodes between is generated a new node
                            nnodbd(:),      & !nnodbd = 0 If created node is not a boundary node
                                              !       = 1 If created node is a boundary node
                            nwlnd(:,:),     & !new conectivities of refined mesh
                            lindi(:)          !label of new element by old elements
  REAL(kind=8),INTENT(IN):: cd(:,:)     !Nodal coordinates
  REAL(kind=8),POINTER:: cvdisp(:,:)    !Displacements of culvature correction

END SUBROUTINE dispcrv
