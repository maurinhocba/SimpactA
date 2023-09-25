      SUBROUTINE fnod_fix(nnodbd,nodcre)
!--------------------------------------------------------------------
!     get ifx conditions for new boundary nodes in a refined set
!--------------------------------------------------------------------
      IMPLICIT NONE
      !--- Dummy variables
      INTEGER(kind=4),POINTER:: nodcre(:,:), & !Nodes between is generated a new node
                                nnodbd(:)      !nnodbd = 0 If created node is not a boundary node
                                               !       = 1 If created node is a boundary node

      END SUBROUTINE fnod_fix



