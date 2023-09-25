      SUBROUTINE velnpo ( npoin, nvelr, ifpre, nesdf, npsdf,      &
     &                   ftsdf, velor, veloc, velnp)
!***********************************************************************
!
!     rearrange velocity vector
!
!***********************************************************************
      IMPLICIT NONE

      INTEGER (kind=4), INTENT(IN) :: npoin, nvelr, ifpre(:,:),  &
     &                                nesdf(:),npsdf(:)
      REAL (kind=8), INTENT(IN)  :: ftsdf(:), veloc(:), velor(:,:)
      REAL (kind=8), INTENT(OUT) :: velnp(:,:)

      END SUBROUTINE velnpo
