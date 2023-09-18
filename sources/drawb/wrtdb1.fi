SUBROUTINE wrtdb1 (del1, ttime, toutd, numtdf, ihflag)

  IMPLICIT NONE

  !   d u m m y   a r g u m e n t s

  INTEGER :: numtdf, ihflag
  REAL (kind=8) :: del1, ttime, toutd

  !   l o c a l   v a r i a b l e s

  REAL (kind=8) :: time1, toler

     !  checking history output

  IF (numtdf == 0) THEN

     ihflag = 0

  ELSE

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
END SUBROUTINE wrtdb1
