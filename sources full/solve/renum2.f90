      SUBROUTINE renum2(adjnt,start,level,lpntn,npoin,ns)
!***********************************************************************
!*
!****  computes a set of posibles starting npoin
!*
!***********************************************************************
      IMPLICIT none
!     routine parameters
      INTEGER (kind=4),INTENT(IN) :: adjnt(1:),npoin
      INTEGER (kind=4),INTENT(OUT):: lpntn(:),level(1:),start(1:),ns
!     local variables
      LOGICAL (kind=4) :: better
      INTEGER (kind=4) :: i,iroot,depth,width,nnoded,ndpth,nwdth,lhw

      INTERFACE
        INCLUDE 'renum4.h'
      END INTERFACE

!.... begin iteration
!.... select intial root node arbitrarily and generate its level
!.... structure

      iroot = 1
      DO

        CALL renum4(start,level,depth,adjnt,width,npoin,iroot,lhw)

        !..  create a list of npoin which are at maximum distance from root
        !..  node and store the root
        ns = iroot

        !..  loop over npoin at maximum distance from root node
        !..  generate level structure for each node
        !..  set switch IF a level structure of greater depth occurs
        better = .FALSE.

        DO i = 1,lhw
          nnoded = start(i)
          CALL renum4(lpntn,level,ndpth,adjnt,nwdth,npoin,nnoded,lhw)
          IF(ndpth >= depth) THEN
            IF((ndpth /= depth) .OR. (nwdth < width)) THEN
              iroot = nnoded
              depth = ndpth
              width = nwdth
              better = .TRUE.
            END IF
          END IF
        END DO
        IF (.NOT.better) EXIT
      END DO

      RETURN

      END SUBROUTINE renum2
