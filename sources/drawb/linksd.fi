SUBROUTINE linksd (irect, lmsr,  msr, nseg, nrtm, &
 &                   nmn, iwork)
  !--------------------------------------------------------------
  !
  !  create the database of sheet segments surrounding every node
  !
  !--------------------------------------------------------------

  IMPLICIT NONE

  !   d u m m y   a r g u m e n t s

  INTEGER, INTENT (IN) :: nrtm         !number of elements in the mesh
  INTEGER, INTENT (IN) :: nmn          !number of nodes in the mesh
  INTEGER, INTENT (IN) :: irect (:,:)  !(4,:) connectivities 
  INTEGER, INTENT (INOUT) :: lmsr (:)  !IN 0
  INTEGER, INTENT (IN) :: msr (:)      !(nmn) nodes in the mesh
  INTEGER, INTENT (IN) :: nseg (:)     !
  INTEGER, INTENT (INOUT) :: iwork (:)
  !-----------------------------------------------
  !   l o c a l   v a r i a b l e s

  INTEGER :: i, j, il, ii, locn, nsgt
  !-----------------------------------------------

  iwork (msr(:nmn)) = 0  !initializes
  DO i = 1, nmn  !generates inverse relation
     iwork (msr(i)) = i
  END DO

  DO i = 1, nrtm  !loop over elements in the mesh
     DO j = 1, 4    !for each node
        IF (irect(j, i) /= 0) THEN   !unnecesary
           il = iwork (irect(j, i))  !node position in local array 
           locn = nseg (il)          !pointer to array LMSR
           nsgt = nseg (il+1) - locn !number of elements around IL node
           ii = COUNT (lmsr(locn:nsgt-1+locn) /= 0) !
           lmsr (locn+ii) = i
        END IF
     END DO
  END DO
  RETURN

END SUBROUTINE linksd
