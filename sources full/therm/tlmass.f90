SUBROUTINE tlmass ()

  ! *** calculates lumped capacity matrix for each element
  !   & assembles global matrix

  USE ctrl_db, ONLY:  ndoft, npoin, tmscal
  USE c_input
  USE npo_db, ONLY : iftmp,label, tmass
  USE outp_db, ONLY: iwrit
  IMPLICIT NONE

  INTEGER (kind=4) :: i,j,ipoin
  REAL    (kind=8) :: xcmas(ndoft)

  INTERFACE
    INCLUDE 'elemnt.h'
  END INTERFACE

  tmass = 0d0  !initializes

  CALL elemnt ('TLMASS')

  IF(iwrit == 1) THEN
    IF(ndoft == 1)WRITE(lures,"(//,'  Lumped Capacity  Matrix' /,          &
                              &    ' Node',6X,'T')",ERR=9999)
    IF(ndoft == 2)WRITE(lures,"(//,'  Lumped Capacity  Matrix' /,' Node',  &
                & 6X,'TI',10X,'TS')",ERR=9999)
    IF(ndoft == 3)WRITE(lures,"(//,'  Lumped Capacity  Matrix' /,' Node',  &
                & 6X,'TN',10X,'TI',10X,'TS')",ERR=9999)
    DO ipoin=1,npoin
      xcmas = 0d0
      DO j=1,ndoft
        i = iftmp(j,ipoin)
        IF( i > 0 ) xcmas(j) = tmass(i)
      END DO
      WRITE(lures,"(i5,3e12.4)",ERR=9999)label(ipoin),xcmas(1:ndoft)
    END DO
  END IF
  tmass = tmass*tmscal  !scales CAPACITY matrix

RETURN
 9999 CALL runen2('')
END SUBROUTINE tlmass
