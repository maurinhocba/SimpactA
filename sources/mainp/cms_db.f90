      MODULE cms_db
        USE lispa0
        IMPLICIT NONE
        SAVE

        ! Derived type for a list containing concentrated masses

        ! list
        TYPE cms_nod
          INTEGER (kind=4) :: node
          REAL (kind=8) :: mass(6)
          TYPE (cms_nod), POINTER :: next
        END TYPE cms_nod

        ! array
        INTEGER (kind=4) :: nconm = 0
        INTEGER (kind=4), POINTER :: nodcms(:)
        REAL (kind=8), POINTER :: cmass(:,:)

      CONTAINS

        SUBROUTINE ini_cms (head, tail)
          !initialize a list

          !Dummy arguments
          TYPE (cms_nod), POINTER :: head, tail

          NULLIFY (head, tail)

        END SUBROUTINE ini_cms

        SUBROUTINE add_cms (new, head, tail)
          !This subroutine adds data to the end of the list
          !Dummy arguments
          TYPE (cms_nod), POINTER :: new, head, tail

          !Check if a list is empty
          IF (.NOT. ASSOCIATED (head)) THEN
            !list is empty, start it
            head => new
            tail => new
            NULLIFY (tail%next)

          ELSE
            !add a segment to the list
            tail%next => new
            NULLIFY (new%next)
            tail => new

          ENDIF
        END SUBROUTINE add_cms

        SUBROUTINE store_cms (head,nrotd)
          ! generates an array with all the data, and release list memory
          IMPLICIT NONE
          !Dummy argument
          TYPE (cms_nod), POINTER :: head
          INTEGER (kind=4) :: nrotd

          !Local variables
          TYPE (cms_nod), POINTER :: ptr,prev
          INTEGER (kind=4) :: n

          ptr => head

          n = 0
          DO
            IF (.NOT.ASSOCIATED(ptr)) EXIT
            n = n + 1
            nodcms(n) = ptr%node
            cmass(1:nrotd,n) = ptr%mass(1:nrotd)
            prev => ptr
            ptr => ptr%next
            DEALLOCATE( prev )
          END DO

          RETURN
          9999 CALL runen2('')
        END SUBROUTINE store_cms

      SUBROUTINE dump_cm (nrotd)

      !   dumps concentrated masses

      IMPLICIT NONE

      INTEGER (kind=4) :: nrotd

      INTEGER (kind=4) :: i

      WRITE (50,ERR=9999) nconm
      IF (nconm > 0) THEN
        DO i=1,nconm
          WRITE (50,ERR=9999) nodcms(i), cmass(1:nrotd,i)
        END DO
      END IF

      RETURN
      9999 CALL runen2('')
      END SUBROUTINE dump_cm

      SUBROUTINE rest_cm (nrotd)

      !   restores concentrated masses at restart

      IMPLICIT NONE
      INTEGER (kind=4) :: nrotd

      INTEGER (kind=4) :: i

      READ (51) nconm
      IF (nconm > 0) THEN
        ALLOCATE ( nodcms(nconm), cmass(nrotd,nconm) )
        DO i=1,nconm
          READ (51) nodcms(i), cmass(1:nrotd,i)
        END DO
      END IF

      RETURN

      END SUBROUTINE rest_cm


      SUBROUTINE updlon_cm (nrotd)

      !   update list of concentrated masses

      IMPLICIT NONE

      INTEGER (kind=4) :: nrotd
      !local
      INTEGER (kind=4) :: n,i,lab,chnode

      IF (nconm > 0) THEN
        n = 0
        DO i=1,nconm
          lab = chnode(nodcms(i))
          IF( lab > 0 )THEN
            n = n + 1
            IF( n == i) CYCLE
            nodcms(n) = nodcms(i)
            cmass(1:nrotd,n) = cmass(1:nrotd,i)
          END IF
        END DO
        nconm = n
      END IF

      RETURN
      9999 CALL runen2('')
      END SUBROUTINE updlon_cm

      END MODULE cms_db
