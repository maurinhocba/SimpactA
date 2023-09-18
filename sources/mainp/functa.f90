      FUNCTION functa  (accer,afact,nacce,ttime)
!********************************************************************
!
!***  accelerogram interpolation
!
!********************************************************************
      IMPLICIT NONE

      INTEGER (kind=4),INTENT(IN) :: nacce
      REAL    (kind=8),INTENT(IN) :: accer(*),afact,ttime
      REAL (kind=8) :: functa

      INTEGER (kind=4) :: n,np1
      REAL    (kind=8) :: t

      IF(ttime <= 0 ) THEN
        functa = 0d0
      ELSE
        t = ttime/afact+1d0
        n   = INT(t)
        IF (n < nacce) THEN
          np1 = n+1
          t   = t - REAL(n,kind=8)
          functa = accer(n)*(1d0-t) + t*accer(np1)
        ELSE
          functa = 0d0
        END IF
      END IF
      RETURN
      END FUNCTION functa
