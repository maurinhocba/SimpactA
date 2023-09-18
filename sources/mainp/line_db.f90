      MODULE line_db
        IMPLICIT NONE
        SAVE

        ! Derived type for the list of 2-node segments defining a line
        TYPE lsg
          INTEGER (kind=4) :: nods(2)
          TYPE (lsg), POINTER :: next
        END TYPE lsg
        TYPE (lsg), POINTER :: headli, tailli

      CONTAINS

        SUBROUTINE ini_lsg (head, tail)
          !initialize the list of segments

          !Dummy arguments
          TYPE (lsg), POINTER :: head, tail

          NULLIFY (head, tail)
        END SUBROUTINE ini_lsg

        SUBROUTINE add_lsg (new, head, tail)
          !add a segment to the end of the list
          !Dummy arguments
          TYPE (lsg), POINTER :: new, head, tail

          !Check if the list is empty
          IF (.NOT. ASSOCIATED (head)) THEN
            !list is empty, start it
            head => new
            tail => new
            NULLIFY (tail%next)

          ELSE
            !add a node to the list
            tail%next => new
            NULLIFY (new%next)
            tail => new

          ENDIF
        END SUBROUTINE add_lsg

        SUBROUTINE srch_left (head, number, found)
          !This subroutine searches for a seg with the 2nd node = number
          !Dummy arguments
          LOGICAL :: found
          INTEGER (kind=4) :: number
          TYPE (lsg), POINTER :: head, posic

          found = .FALSE.
          NULLIFY (posic)
          !Check if the list is empty
          IF (ASSOCIATED (head)) THEN
            posic => head
            DO
              IF(posic%nods(2) == number) THEN
                found = .TRUE.
                EXIT
              END IF
              IF (ASSOCIATED(posic%next) ) THEN
                posic => posic%next
              ELSE
                EXIT
              END IF
            END DO
          ENDIF
          IF (.NOT.found) NULLIFY (posic)
        END SUBROUTINE srch_left

        SUBROUTINE srch_lsg (head, left, right, number, found)
          !This subroutine searches for a seg with the 1st node = number
          !Dummy arguments
          LOGICAL :: found
          INTEGER (kind=4) :: number
          TYPE (lsg), POINTER :: head, left, right

          found = .FALSE.
          NULLIFY (right,left)
          !Check if the list is empty
          IF (ASSOCIATED (head)) THEN
            right => head
            DO
              IF(right%nods(1) == number) THEN
                found = .TRUE.
                EXIT
              END IF
              IF (ASSOCIATED(right%next) ) THEN
                left => right
                right => right%next
              ELSE
                EXIT
              END IF
            END DO
          ENDIF
          IF (.NOT.found) NULLIFY (right,left)
        END SUBROUTINE srch_lsg

        SUBROUTINE mov2h_lsg (head, tail, anter, posic)
          !moves a lsg pointed with posic to the head
          TYPE (lsg), POINTER :: head, tail, anter, posic

          IF (.NOT.ASSOCIATED (anter)) RETURN
          anter%next => posic%next
          ! if posic == tail
          IF (.NOT.ASSOCIATED (posic%next) ) tail => anter
          posic%next => head
          head => posic
        END SUBROUTINE mov2h_lsg

        SUBROUTINE mov_lsg (tail, anter, posic, left, right)
          !moves a lsg pointed with posic after left
          TYPE (lsg), POINTER :: tail, anter, posic, left, right

          anter%next => posic%next
          IF (.NOT.ASSOCIATED (posic%next) ) tail => anter! if posic == tail
          left%next => posic
          posic%next => right
        END SUBROUTINE mov_lsg

        SUBROUTINE sort_lsg (head, tail, sorted, closed)
          ! sorts a list of lsg in a chain lsg%nods(2) = lsg%next%nods(1)
          USE lispa0 ,ONLY:  lures
          USE outp_db ,ONLY:  iwrit
          !Dummy arguments
          LOGICAL :: sorted, closed
          TYPE (lsg), POINTER :: head, tail

          LOGICAL :: found
          TYPE (lsg), POINTER :: anter, posic, left, right, aux, seg

          closed = .TRUE.  ! closed line assumed initially
          NULLIFY (anter)
          seg => head
          DO
            IF (.NOT.ASSOCIATED (seg) ) EXIT ! IF the line is closed
            ! search left extreme
            CALL srch_left (head, seg%nods(1), found)
            IF (found) THEN
              anter => seg
              seg => seg%next
            ELSE
              CALL mov2h_lsg (head, tail, anter, seg)
              closed = .FALSE.  ! open line found
              EXIT
            ENDIF
          END DO

          sorted = .FALSE.
          NULLIFY (posic,anter,left,right)
          !Check if the list is empty
          IF (ASSOCIATED (head)) THEN
            left => head
            DO
              right => left%next
              IF (ASSOCIATED (right)) THEN
                IF(left%nods(2) == right%nods(1)) THEN
                  left => right
                ELSE
                  found = .FALSE.
                  CALL srch_lsg (head, anter, posic, left%nods(2),found)
                  IF (found) THEN
                    aux => posic
                    CALL mov_lsg (tail, anter, posic, left, right)
                    left => aux
                  ELSE
                    CALL runend ('SORT_LSG: Line not continous.      ')
                  END IF
                ENDIF
              ELSE
                sorted = .TRUE.
                EXIT
              END IF
            END DO
          ENDIF
          ! control printout
          seg => head
          DO
            IF (.NOT.ASSOCIATED (seg) ) EXIT ! IF the line is closed
            ! search left extreme
            IF(iwrit > 0 )WRITE (lures,*,ERR=9999) seg%nods
            seg => seg%next
          END DO
        RETURN
   9999 CALL runen2('')
        END SUBROUTINE sort_lsg

      END MODULE line_db
