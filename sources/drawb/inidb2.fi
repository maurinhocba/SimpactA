SUBROUTINE inidb2 (irect, msr, iwork,  nseg, nmn, nrt, mst)
  !-----------------------------------------------------------------
  !                initialization of the sheet database
  !-----------------------------------------------------------------

  IMPLICIT NONE

  !   d u m m y   a r g u m e n t s

  INTEGER, INTENT(INOUT)  :: nmn  !number of nodes on elements
  INTEGER, INTENT(IN)  :: nrt  !number of elements to check
  INTEGER, INTENT (INOUT) :: mst !
  INTEGER, INTENT(IN) :: irect (:, :) !connectivities of element to check
  INTEGER :: msr (:)
  INTEGER :: iwork (:) !auxiliar array
  INTEGER :: nseg (:)  !nmn+1  pointers to an array

  INTERFACE
    INCLUDE 'sortnd.h'
  END INTERFACE

  !   l o c a l   v a r i a b l e s

  INTEGER :: j
  !-----------------------------------------------

  ! find nodes in sheet
  CALL countn (irect, msr, nmn, nrt, iwork) !=>msr(nmn)
  CALL sortnd (msr, nmn)    !  sort sheet nodes

  !  count sheet segments surrounding every node
  nseg = 0      !initializes pointers
  CALL counts (irect, msr, nseg(2:), nrt, nmn, iwork)
  ! generate of pointers
  nseg (1) = 1  !first pointer points to first position
  DO j = 2, nmn+1     !for each node in the mesh
     nseg (j) = nseg (j-1) + nseg (j)    !
  END DO
  mst = nseg(nmn+1) - 1
  RETURN
END SUBROUTINE inidb2
