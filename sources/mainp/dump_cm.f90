      SUBROUTINE dump_cm (ndofn)

      !   dumps concentrated masses

      USE cms_db
      IMPLICIT NONE

      INTEGER (kind=4) :: ndofn

      INTEGER (kind=4) :: i

      WRITE (50,ERR=9999) nconm
      IF (nconm > 0) THEN
        DO i=1,nconm
          WRITE (50,ERR=9999) nodcms(i), cmass(1:ndofn,i)
        END DO
      END IF

      RETURN
 9999 CALL runen2('')
      END SUBROUTINE dump_cm

