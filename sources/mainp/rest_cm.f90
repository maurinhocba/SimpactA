      SUBROUTINE rest_cm (ndofn)

      !   restores concentrated masses at restart

      USE cms_db
      IMPLICIT NONE
      INTEGER (kind=4) :: ndofn

      INTEGER (kind=4) :: i

      READ (51) nconm
      IF (nconm > 0) THEN
        ALLOCATE ( nodcms(nconm), cmass(ndofn,nconm) )
        DO i=1,nconm
          READ (51) nodcms(i), cmass(1:ndofn,i)
        END DO
      END IF

      RETURN

      END SUBROUTINE rest_cm

