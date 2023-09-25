      SUBROUTINE renum3(adjnt,start,level,lpntn,npoin,i)
!***********************************************************************
!*     This routine MUST not be compiled with Optimization beyond level(1)
!***** resequence npoin for minimum profile
!*
!***********************************************************************
      IMPLICIT none
!        routine parameters
      INTEGER (kind=4),INTENT(IN) :: npoin,i,adjnt(1:)
      INTEGER (kind=4),INTENT(OUT):: lpntn(:),start(1:),level(1:)
!        local variables
      INTEGER (kind=4) :: large,j,nac,maxfrt,nad,npj,k,minnew,lmin,iact,&
     &                    new,min,n,nif,npos,ilast,next

      large = 10000

!.... KING'S SCHEME
      lpntn(1:npoin) = 0
      level(1:npoin) = 0
      lpntn(i) = 1
      nac = 0

!.... negate all ndeg entries for npoin which are
!.... adjacent to starting node i

      maxfrt = 0
      nad = adjnt(i)
      DO WHILE (nad > 0)
        maxfrt = maxfrt + 1
        npj = adjnt(nad)
        IF(level(npj) == 0) THEN
          nac = nac + 1
          level(npj) = nac
          start(nac) = npj
        END IF
        nad = adjnt(nad+1)
      END DO
      level(i) = large

!.... loop over npoin to be renumbered

      DO k = 2,npoin
        minnew = large
        lmin = large

        !....   loop over active npoin
        !....   skip to next node IF old node is already renumbered

        DO iact = 1,nac
          j = start(iact)
          IF(lpntn(j) <= 0) THEN
            new = -1
            min = large

            !  compute the increment in active npoin for each node j
            !  compute when this node was first active by checking for
            !  renumbered neighbours with lowest numbers

            nad = adjnt(j)
            DO WHILE (nad > 0)
              n = adjnt(nad)
              IF(level(n) == 0) new = new + 1
              IF(lpntn(n) /= 0) THEN
                IF(lpntn(n) < min) min = lpntn(n)
              END IF
              nad = adjnt(nad+1)
            END DO

            !  select node with smallest increment in active npoin
            !  in the CASE of a tie, SELECT node which has been longest active

            IF(new <= minnew) THEN
              IF((new /= minnew) .OR. (min < lmin)) THEN
                minnew = new
                lmin = min
                next = j
              END IF
            END IF
          END IF
        END DO

        !....   renumber node and compute number of active npoin

        lpntn(next) = k
        nif = nif+minnew
        IF(nif > maxfrt) maxfrt = nif

        !....   set npoin which are adjacent to the node just renumbered
        !....   as actives npoin, deactivate next

        npos = level(next)
        ilast = start(nac)
        level(ilast) = npos
        start(npos) = ilast
        nac = nac - 1

        IF(minnew  /=  -1) THEN
          nad = ABS(adjnt(next))
          DO WHILE (nad > 0)
            n = adjnt(nad)
            IF(level(n) == 0) THEN
              nac = nac + 1
              level(n) = nac
              start(nac) = n
            END IF
            nad = adjnt(nad+1)
          END DO
        END IF
      END DO

      RETURN

      END SUBROUTINE renum3
