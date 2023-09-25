 MODULE ift_db
   ! this MODULE contains a data base to manage Temperature boundary conditions

   USE ctrl_db, ONLY : ndoft  !number of temperature DOFs per node
   USE param_db,ONLY: mnam    !lenght of labels

   IMPLICIT NONE
   INTEGER(kind=4), PARAMETER :: nmf=3  !(>= NDOFT)
   SAVE

   !******* for fixed values

   ! Derived type for fixity data
   ! list
   TYPE ift_nod
     INTEGER (kind=4) :: ifix(1+nmf)      !1:label 2-4: fix codes
     REAL (kind=8) :: val(nmf)            !fixed values
     TYPE (ift_nod), POINTER :: next      !pointer to next node in the list
   END TYPE ift_nod

   TYPE (ift_nod), POINTER :: head,tail    ! pointers to first and last values

   INTEGER (kind=4) :: nift=0, &  !Number of nodes in the list
                       ntfix      !number of total fixed values

  REAL (kind=8),POINTER :: fxtem(:)       ! (NTFIX) fixed temperatures

   !******* for prescribed values (changing in time with a known curve)

  ! node in a list
  TYPE rpt_nod
    INTEGER (kind=4) :: node        !node label
    REAL (kind=8) :: v(nmf)         !fixed temperatures
    TYPE (rpt_nod), POINTER :: next !next position
  END TYPE rpt_nod

  ! set (list of nodes)
  TYPE rpt_set
    CHARACTER (len=mnam) :: lc                !associated curpt label
    INTEGER (kind=4) :: nrv                   !number of nodes in the set
    REAL (kind=8) :: factor                   !curpt factor
    TYPE (rpt_nod), POINTER :: head,tail      !pointers to first and last
    TYPE (rpt_set), POINTER :: next           !pointer to next set
  END TYPE rpt_set

  ! global pointers
  TYPE (rpt_set), POINTER :: headv,tailv   !pointers to first & last sets

  INTEGER (kind=4) :: &
    nprev,            &  !number of prescribed values
    npret                !number of constrained temperature sets

  INTEGER (kind=4),POINTER :: &
    lctmp(:)      ! associated curves to prescribed temp. sets

  REAL (kind=8),POINTER :: &
    prtmp(:,:)     ! prescribed temperatures in time


 CONTAINS

   SUBROUTINE ini_ift (head, tail)
     !initialize a dependent nodes list

     !Dummy arguments
     TYPE (ift_nod), POINTER :: head, tail

     NULLIFY (head, tail)

   END SUBROUTINE ini_ift

   SUBROUTINE add_ift (new, head, tail)
     !This subroutine adds data to the end of the list
     !Dummy arguments
     TYPE (ift_nod), POINTER :: new, head, tail

     !Check if a list is empty
     IF (.NOT. ASSOCIATED (head)) THEN
       !list is empty, start it
       head => new
       tail => new
       NULLIFY (tail%next)

     ELSE
       !add data to the list
       tail%next => new
       NULLIFY (new%next)
       tail => new

     END IF
   END SUBROUTINE add_ift

   SUBROUTINE srch_ift (head, anter, posic, node, found)
     !This subroutine searches the pointer for a "node" ==> POSIC
     !Dummy arguments
     LOGICAL :: found               !flag
     INTEGER (kind=4) :: node       !node label
     TYPE (ift_nod), POINTER :: head, anter, posic !head of list ,previous and search result

     found = .FALSE.                         !initializes flag and pointers
     NULLIFY (posic,anter)                   !
     !Check if a list is empty
     IF (ASSOCIATED (head)) THEN
       posic => head                         !point to first
       DO
         IF(posic%ifix(1) == node) THEN      !check label
           found = .TRUE.                    ! O.K.
           EXIT
         END IF
         IF (ASSOCIATED(posic%next) ) THEN   !if there are more points
           anter => posic                    !keep previous position
           posic => posic%next               !point to next
         ELSE
           EXIT                              !no more points EXIT
         END IF
       END DO
     END IF
     IF (.NOT.found) NULLIFY (posic,anter)   !not found, nullify pointers
     RETURN
   END SUBROUTINE srch_ift

   SUBROUTINE del_ift (head, tail, anter, posic)

     !This subroutine deletes data pointed with posic and nullifies posic & anter

     TYPE (ift_nod), POINTER :: head, tail, anter, posic

     IF (.NOT.ASSOCIATED (anter)) THEN   !see if we are at the beginning
       head => posic%next             !then modifies head pointer
     ELSE
       anter%next => posic%next       !else skip posic
     END IF
              ! if posic == tail
     IF (.NOT.ASSOCIATED (posic%next) ) tail => anter
     DEALLOCATE (posic)               !release memory
     NULLIFY (anter)                  !both posic and anter are null now
     RETURN
   END SUBROUTINE del_ift

   SUBROUTINE dump_ift
   IMPLICIT NONE
     !dumps data

     !Local variables
     TYPE (ift_nod), POINTER :: ift
     !INTEGER :: i,n,node

     WRITE (50,ERR=9999) nift       !nift = Nodes with resctrictions
     IF (nift > 0) THEN        !if restrictions exist

       ift => head             !point to first

       DO
         WRITE (50,ERR=9999) ift%ifix(1:ndoft+1),ift%val(1:ndoft)     !writes label, codes and values
         ift => ift%next                  !point to next
         IF (.NOT.ASSOCIATED(ift)) EXIT
       END DO
     END IF
     RETURN
 9999 CALL runen2('')
   END SUBROUTINE dump_ift

   SUBROUTINE rest_ift
     !restores variables from disk
     !Local variables
     TYPE (ift_nod), POINTER :: ift
     INTEGER :: i

     !Initialize empty list
     CALL ini_ift(head,tail)

     READ (51) nift    !nd1 = NDOFN+(1) , nift = Nodes with resctrictions
                           !read list
     DO i = 1,nift
       ALLOCATE (ift)                   !get memory
       READ (51) ift%ifix(1:ndoft+1),ift%val(1:ndoft)      !read data
       CALL add_ift( ift, head, tail )  !add to list
     END DO

     RETURN
   END SUBROUTINE rest_ift

        ! functions for rpt nodes

        SUBROUTINE ini_rptn (head, tail)
          !initialize a list of nodes

          !Dummy arguments
          TYPE (rpt_nod), POINTER :: head, tail

          NULLIFY (head, tail)

        END SUBROUTINE ini_rptn

        SUBROUTINE add_rptn (new, head, tail)
          !This subroutine adds data to the end of the list of nodes
          !Dummy arguments
          TYPE (rpt_nod), POINTER :: new, head, tail

          !Check if a list is empty
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

          END IF
          RETURN
        END SUBROUTINE add_rptn

        ! functions for rpt sets

        SUBROUTINE ini_rpts (head, tail)
          !initialize a list

          !Dummy arguments
          TYPE (rpt_set), POINTER :: head, tail

          NULLIFY (head, tail)

        END SUBROUTINE ini_rpts

        SUBROUTINE add_rpts (new, head, tail)
          !This subroutine adds a set to the end of the list of sets
          !Dummy arguments
          TYPE (rpt_set), POINTER :: new, head, tail

          !Check if a list is empty
          IF (.NOT. ASSOCIATED (head)) THEN
            !list is empty, start it
            head => new
            tail => new
            NULLIFY (tail%next)

          ELSE
            !add a set to the list
            tail%next => new
            NULLIFY (new%next)
            tail => new

          END IF
          RETURN
        END SUBROUTINE add_rpts

        SUBROUTINE srch_rpts (head, anter, posic, lc, found)
          !This subroutine searches for a set associated to curpt LC
          !Dummy arguments
          LOGICAL :: found            ! Answer
          CHARACTER (len=mnam) :: lc  ! curpt label
          TYPE (rpt_set), POINTER :: head, anter, posic

          found = .FALSE.                   !initializes
          NULLIFY (posic,anter)
          !Check if a list is empty
          IF (ASSOCIATED (head)) THEN
            posic => head                   !point to head of the list
            DO
              IF(posic%lc == lc) THEN       !compare with associated curpt
                found = .TRUE.              !found
                EXIT                        !leave search
              END IF
              IF (ASSOCIATED(posic%next) ) THEN     !more sets to check
                anter => posic                      !keep previous
                posic => posic%next                 !point to next
              ELSE
                EXIT                        !no more sets => exit
              END IF
            END DO
          ENDIF
          IF (.NOT.found) NULLIFY (posic,anter)   !not found => nullify pointers
          RETURN
        END SUBROUTINE srch_rpts

        SUBROUTINE del_rpts (head, tail, anter, posic)
          !deletes a set pointed with posic
          TYPE (rpt_set), POINTER :: head, tail, anter, posic

          IF (.NOT.ASSOCIATED (anter)) THEN  !if first set
            head => posic%next               !new head
          ELSE
            anter%next => posic%next         !skip posic
          END IF
          ! if posic == tail    (last set)
          IF (.NOT.ASSOCIATED (posic%next) ) tail => anter  !new last set
          CALL dalloc_rpts (posic)           !release memory
          NULLIFY (anter)                    !both anter & posic are null now
        END SUBROUTINE del_rpts

        SUBROUTINE dalloc_rpts (rpts)
          ! deallocates a rpt set (release memory)
          TYPE (rpt_set), POINTER :: rpts
          TYPE (rpt_nod), POINTER :: rptn, rptnaux

          rptn => rpts%head    !point to first node
          DO
            IF (.NOT.ASSOCIATED (rptn) ) EXIT
            rptnaux => rptn%next    !keep next pointer
            DEALLOCATE (rptn)       !release memory of the node
            rptn => rptnaux         !point to next
          END DO

          DEALLOCATE (rpts)         !release rest of the vars. lc nrv & factor
          RETURN
        END SUBROUTINE dalloc_rpts

        SUBROUTINE store_rpt (head,vel,nods,ndoft)
          !This subroutine stores the data

          !Dummy argument
          TYPE (rpt_nod), POINTER :: head
          INTEGER (kind=4) :: ndoft, nods(:)
          REAL (kind=8) :: vel(:,:)

          !Local variables
          TYPE (rpt_nod), POINTER :: ptr
          INTEGER :: n

          IF (ASSOCIATED(head))THEN
            ptr => head

            n = 0
            DO
              n = n + 1
              nods(n) = ptr%node
              vel(1:ndoft,n) = ptr%v(1:ndoft)
              ptr => ptr%next
              IF (.NOT.ASSOCIATED(ptr)) EXIT
            END DO
          ENDIF
        END SUBROUTINE store_rpt

 END MODULE ift_db
