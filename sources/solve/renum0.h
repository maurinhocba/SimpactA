      SUBROUTINE renum0(lnods,lpntn,nelem,nnode,npoin,krenu)
      IMPLICIT none
      INTEGER (kind=4),INTENT(IN) :: nelem,nnode,npoin,krenu
      INTEGER (kind=4),INTENT(IN OUT) :: lnods(:,:)
      INTEGER (kind=4),INTENT(OUT):: lpntn(:)
      END SUBROUTINE renum0
