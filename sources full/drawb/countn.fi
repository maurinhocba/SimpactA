      SUBROUTINE countn (irect, mnn, n, nrt, iscr)

      !     Routine due to JR, used in CONT3 and DRAWB
      !     an identical routines exists in HYDF_DB.F90
      !     COUNTs Nodes in a connectivities array

      IMPLICIT NONE

      !   d u m m y   a r g u m e n t s

      INTEGER, INTENT (OUT) :: n            !number of nodes in the mesh
      INTEGER, INTENT (IN) :: nrt           !number of segments
      INTEGER, INTENT (IN) :: irect (:,:)   !(4,nrt) surface connectivities
      INTEGER, INTENT (OUT) :: mnn (:)      !(npoin) nodes included in the mesh
      INTEGER :: iscr (:)     !(npoin) auxiliar (0 at output)
      !-----------------------------------------------
      !   l o c a l   v a r i a b l e s

      INTEGER :: i,j
      !-----------------------------------------------


      IF (nrt == 0) RETURN !no segments
      n = 0                !initializes counter
      iscr = 0             !initializes array
      DO i = 1, nrt  !number of segments
        DO j = 1, 4    !number of nodes per segment
         IF (iscr(irect(j,i)) == 0) THEN !not considered yet
            n = n + 1                    !increase
            mnn (n) = irect (j,i)        !add node
            iscr (irect(j,i)) = n        !generate inverse relation
         END IF
        END DO
      END DO
      iscr (mnn(:n)) = 0 !back to 0

      RETURN
      END SUBROUTINE countn
