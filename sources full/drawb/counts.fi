      SUBROUTINE counts (irect, mnn, nseg, nrt, nn, iscr)

      IMPLICIT NONE

      !   d u m m y   a r g u m e n t s

      INTEGER, INTENT (IN) :: nrt          !number of elements
      INTEGER, INTENT (IN) :: nn           !number of nodes in the mesh
      INTEGER, INTENT (IN) :: irect (:,:)  !(4,nrt) connectivities
      INTEGER, INTENT (IN) :: mnn (:)      !(nn) nodes in the mesh
      INTEGER, INTENT (IN OUT) :: nseg (:) !(nn) IN = 0 OUT number of connected elements
      INTEGER :: iscr (:)    !(npoin) auxiliar array IN=0  OUT=0
      !-----------------------------------------------
      !   l o c a l   v a r i a b l e s

      INTEGER :: i,j,k
      !-----------------------------------------------
      ! generate inverse relation MNN ==> ISCR 
      DO i = 1, nn
         iscr (mnn(i)) = i
      END DO
      ! loop over elements to be checked
      DO i = 1, nrt  !for each element
	    DO j=1,4 !why not NNODE defined in database 
		         !for 3-node elements last node appears twice 
          k = irect(j,i)
          IF ( k /= 0) nseg(iscr(k)) = nseg(iscr(k)) + 1
        END DO
      END DO

      iscr (mnn(:nn)) = 0 !back to 0

      RETURN
      END SUBROUTINE counts
