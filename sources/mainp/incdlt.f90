      SUBROUTINE incdlt( )
      !********************************************************************
      !
      !***   evaluation of time increment
      !
      !********************************************************************
      USE ctrl_db, ONLY : dtscal,dtuser,dtime,ttime,endtm,ncdlt,tdtuse,tdtime,tdtsca,therm
      USE lispa0, ONLY : lures
      IMPLICIT NONE

      REAL (kind=8) :: deltc,tdelt
      INTEGER(kind=4) :: n
      INTERFACE
        INCLUDE 'elemnt.h'
      END INTERFACE

      !***  loop over all the sets


      deltc = 1.0d+30        !set to a large value
      CALL elemnt ('INCDLT', deltc= deltc) !compute elements Maximum Value
      WRITE(lures,"(' Critical Time Increment',e15.5)")deltc
      IF(dtuser == 0d0)THEN        !check if the user gave a Maximum value
        dtime = deltc*dtscal       !use Maximum value
      ELSE
        dtime = MIN(deltc*dtscal,dtuser) !compare with user value
      END IF
      ! modify dtime to arrive exactly to ENDTM
      IF( ncdlt > 9 )THEN     !but only for ncdlt >= 10
        n = MAX(NINT((endtm-ttime)/dtime),1)      !steps to end of the strategy
        IF( n < ncdlt ) dtime = (endtm-ttime)/n   !modify dtime if end is near
      END IF
      WRITE(lures,"(' Time Increment',e15.5)")dtime

      IF( therm )THEN
        tdelt = 1.0d+30
        CALL elemnt ('TINCDT', deltc= tdelt)
        IF(tdtuse == 0d0)THEN
          tdtime = tdelt*tdtsca
        ELSE
          tdtime = MIN(tdelt*tdtsca,tdtuse)
        END IF
        IF(tdtime > dtime ) THEN
          tdtime = dtime  !O.K.
        ELSE
          dtime  = tdtime !BAD
        END IF
      END IF

      RETURN
      END SUBROUTINE incdlt
