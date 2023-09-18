      SUBROUTINE wrtdb2 (idbfo, ndrawbead, del1, ttime, toutd, ihflag)
!-----------------------------------------------
!   m o d u l e s
!-----------------------------------------------
      USE c_input, ONLY : openfi
!-----------------------------------------------
!   i n t e r f a c e   b l o c k s
!-----------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------
!   d u m m y   a r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT (IN) :: ndrawbead
      INTEGER, INTENT (OUT) :: ihflag
      REAL (kind=8), INTENT (IN) :: del1
      REAL (kind=8), INTENT (IN) :: ttime
      REAL (kind=8), INTENT (IN) :: toutd
      INTEGER, INTENT (IN) :: idbfo (ndrawbead)
!-----------------------------------------------
!   l o c a l   v a r i a b l e s
!-----------------------------------------------
      INTEGER :: ndbfout, iprnums, i
      REAL (kind=8) :: time1, toler

      SAVE iprnums
!-----------------------------------------------
      DATA iprnums / 0 /

      ihflag = 0
      IF (iprnums == (-1)) RETURN

      IF (iprnums == 0) THEN
         ndbfout = 0

         ndbfout = ndbfout + count (idbfo(:) == 1)

         IF (ndbfout > 0) THEN
            CALL openfi (42)
            WRITE (42,ERR=9999) ndbfout, 1
            DO i = 1, ndrawbead
               IF (idbfo(i) == 1) write (42,ERR=9999) i
            END DO
            iprnums = 1
         ELSE
            iprnums = - 1
         END IF
      END IF

      IF (iprnums == 1) THEN

         !  checking history output

         toler = del1 * 1.0000000001d0 / 2d0
         time1 = ttime - Int (ttime/toutd) * toutd
         IF (time1+toler > toutd .OR. time1 < toler) THEN
            ihflag = 1
            WRITE (42,ERR=9999) ttime
         ELSE
            ihflag = 0
         END IF
      END IF

      RETURN
 9999 CALL runen2('')
      END SUBROUTINE wrtdb2
