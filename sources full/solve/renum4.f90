      SUBROUTINE renum4(start,level,depth,adjnt,width,npoin,iroot,lhw)
!***********************************************************************
!*
!*****compute level structure rooted at iroot
!*
!***********************************************************************
      IMPLICIT none
!     routine parameters
      INTEGER (kind=4),INTENT(IN) :: npoin,iroot,adjnt(:)
      INTEGER (kind=4),INTENT(OUT):: depth,width,lhw,level(:),start(:)
!     local variables
      LOGICAL (kind=4) :: back
      INTEGER (kind=4) :: i,kount,lhwn,il,ip,nad,jp,nhalf,nn

!.... initialization

      level(1:npoin) = 0

      level(iroot) = 1
      start(1) = iroot
      back = .false.
      depth = 1
      kount = 1
      lhw = 1
      lhwn = 0
      width = 1

!.... assign levels

      DO
        IF(kount >= npoin)EXIT
        DO il = 1,lhw
          IF(.NOT.back) ip = start(il)
          IF(back) ip = start(npoin+1-il)
          nad = adjnt(ip)
          DO 
            IF(nad <= 0)EXIT
            jp = adjnt(nad)
            IF(level(jp) == 0) THEN
              lhwn = lhwn + 1
              level(jp) = depth + 1
              kount = kount + 1
              IF(back) start(lhwn) = jp
              IF(.NOT.back) start(npoin+1-lhwn) = jp
            END IF
            nad = adjnt(nad+1)
          END DO
        END DO
        lhw = lhwn
        lhwn = 0
        IF(lhw > width) width = lhw
        depth = depth + 1
        back = .NOT.back
      END DO

      IF(back) THEN
        nhalf = npoin/2
        DO i = 1,nhalf
          nn = start(i)
          start(i) = start(npoin+1-i)
          start(npoin+1-i) = nn
        END DO
      END IF

      RETURN

      END SUBROUTINE renum4
