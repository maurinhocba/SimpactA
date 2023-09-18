      SUBROUTINE reasdf(nescv,iwrit,npsdf,ifpre,nesdf,ftsdf,
     &                  label,npoin,ndofn)

!     READ restrictions on (translational) degrees of freedom

      IMPLICIT NONE
      INTEGER (kind=4),INTENT(IN) :: iwrit,npoin,label(npoin),ndofn
      INTEGER (kind=4),INTENT(IN OUT) :: ifpre(ndofn,*),nescv
      INTEGER (kind=4),INTENT(OUT) :: npsdf(*),nesdf(*)
      REAL (kind=8),INTENT(OUT) :: ftsdf(*)
      END SUBROUTINE reasdf
