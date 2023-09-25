      SUBROUTINE renum1(nnode,npoin,nelem,lnods,nodad,nposi)
      IMPLICIT none
      INTEGER (kind=4),INTENT(IN) :: nnode,npoin,nelem,lnods(:,:)
      INTEGER (kind=4),INTENT(OUT):: nodad(:),nposi
      END SUBROUTINE renum1
