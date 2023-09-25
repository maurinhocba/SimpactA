      SUBROUTINE ensmal(nvarl,neq,lm,locst,glost,npsdf,nesdf,ftsdf)
!*************************************************************************
!
!     assembles a local diagonal matrix into a global diagonal matrix
!     both matrices are stored as arrays
!
!*************************************************************************
      IMPLICIT NONE
!      INTEGER (kind=4),INTENT(IN) :: nvarl,neq,lm(nvarl),npsdf(:),      &
!     &                 nesdf(:)
!      REAL (kind=8),INTENT(IN) :: locst(nvarl),ftsdf(:)
!      REAL (kind=8),INTENT(IN OUT) :: glost(:)
      INTEGER (kind=4),INTENT(IN) :: nvarl,neq,lm(nvarl),npsdf(*),      &
     &                 nesdf(*)
      REAL (kind=8),INTENT(IN) :: locst(nvarl),ftsdf(*)
      REAL (kind=8),INTENT(IN OUT) :: glost(neq)

      END SUBROUTINE ensmal
