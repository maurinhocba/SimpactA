      SUBROUTINE rdvelr(nvelr,ndime,ndofn,npoin,iwrit,label,ifpre,nvfix, &
     &                  lcvel,velor)

  !***  Applies (?) rigid body velocities
      IMPLICIT NONE
      INTEGER (kind=4), INTENT(IN) :: ndime,ndofn,npoin,iwrit,label(:)
      INTEGER (kind=4), INTENT(IN OUT) :: ifpre(:,:),nvfix
      INTEGER (kind=4), INTENT(IN) :: nvelr
      INTEGER (kind=4), POINTER :: lcvel(:)
      REAL (kind=8), POINTER :: velor(:,:)
      END SUBROUTINE rdvelr
