SUBROUTINE db3udl ( )

  !                3d drawbead UpDate internal numbers

  USE npo_db, ONLY: oldlb
  IMPLICIT NONE

  INTEGER (KIND=4):: nrtm, nmn
  LOGICAL :: keep

  TYPE (db3pair_db), POINTER :: pair,prev

    ! updating internal node numbers in old pairs

    pair => headp
    NULLIFY(prev)
    DO
       IF ( .NOT.ASSOCIATED (pair) )EXIT
       nrtm = pair%iparm(2)
       nmn  = pair%iparm(3)
       CALL updnd3 (oldlb, nmn, nrtm,pair%msr, pair%irect, keep)
       IF( keep )THEN
         pair%iparm(2) = nrtm
         pair%iparm(3) = nmn
         prev => pair
         pair => pair%next
       ELSE
         ndbea = ndbea - 1
         DEALLOCATE(pair%msr,pair%irect)
         IF( ASSOCIATED(prev) )THEN
           prev%next => pair%next
           DEALLOCATE(pair)
           pair => prev%next
         ELSE
           headp => pair%next
           DEALLOCATE(pair)
           pair => headp
         END IF
       END IF
    END DO

  RETURN
END SUBROUTINE db3udl
