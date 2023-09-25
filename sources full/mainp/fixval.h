      SUBROUTINE fixval(iwrit,ndime,ndofn,ifpre,nvfix,rot_free,iffix,label)

  !***  APPLY fixities

      IMPLICIT NONE
      LOGICAL, INTENT(IN) :: rot_free
      INTEGER (kind=4), INTENT(IN) :: iwrit,ndime,ndofn,label(:)
      INTEGER (kind=4), INTENT(IN OUT) :: ifpre(:,:),iffix(:)
      INTEGER (kind=4), INTENT(OUT) :: nvfix
      END SUBROUTINE fixval
