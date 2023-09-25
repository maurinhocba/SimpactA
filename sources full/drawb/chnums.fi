      SUBROUTINE chnums (nnode, heade, nelem, label)

      !changes the internal to external numbering

      USE surf_db
      IMPLICIT NONE
      INTEGER(kind=4) :: nelem,nnode,label(:)
      TYPE (srf_seg), POINTER  :: heade

      ! Local 
      INTEGER(kind=4) :: i
      TYPE (srf_seg), POINTER :: seg

      nelem = 0
      seg => heade
      DO 
        ! loop over segments

        IF (.NOT. ASSOCIATED(seg) ) EXIT
        nelem = nelem + 1
        DO i = 1,nnode
          IF (seg%nodes(i) /= 0) seg%nodes(i) = label(seg%nodes(i))
        END DO
        seg => seg%next

      END DO

      END SUBROUTINE chnums
