      SUBROUTINE ensvec(nvarl,lm,locvc,glovc,npsdf,nesdf,ftsdf)
      IMPLICIT NONE
      INTEGER (kind=4),INTENT(IN) :: nvarl,lm(nvarl),npsdf(:),nesdf(:)
      REAL (kind=8),INTENT(IN) :: locvc(nvarl),ftsdf(:)
      REAL (kind=8),INTENT(IN OUT) :: glovc(:)
      END SUBROUTINE ensvec
