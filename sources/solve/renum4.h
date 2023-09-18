       SUBROUTINE renum4(start,level,depth,adjnt,width,npoin,iroot,lhw)
       IMPLICIT none
       INTEGER (kind=4),INTENT(IN) :: npoin,iroot,adjnt(:)
       INTEGER (kind=4),INTENT(OUT):: depth,width,lhw,level(:),start(:)
       END SUBROUTINE renum4
