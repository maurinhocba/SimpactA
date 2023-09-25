 MODULE ele02_db
   USE param_db,ONLY: mnam,midn
   USE mat_dba, ONLY: section,sect_search,mater,psecs
   USE c_input, ONLY : openfi,listen,exists,param,lures,getint,getrea,rdreqs
   IMPLICIT NONE
   INTEGER, PARAMETER :: nnode = 2

   SAVE

   ! Derived type for a  TRUSS element
   TYPE ele02
     INTEGER (kind=4) :: numel  ! label of element
     INTEGER (kind=4) :: matno  ! Material number
     INTEGER (kind=4) :: lnods(nnode)  ! Conectivities
     REAL (kind=8), POINTER    :: gausv(:)      !gaussian variables
     ! for each standard element  (1) l0   : original length
                                 !(2) S0   : initial Force
                                 !(3) S    : present Stress
                                 !(4) N    : present Force
                                 !(5) Ep   : plastic strain (Present)
                                 !(6) BS   : back stress (P)
                                 !(7) Eps  : Effective plastic strain(P)
                                 !(8) Nefpst :
     TYPE (ele02), POINTER :: next              !pointer to next element
   END TYPE ele02
   ! Derived type for a set of TRUSS elements
   TYPE ele02_set
     CHARACTER (len=mnam) :: sname ! set name
     INTEGER (kind=4) :: nelem, &  ! number of elements
                         nreqs, &  ! number of GP for hist. output
                         narch     ! number of output unit
     INTEGER (kind=4), POINTER :: ngrqs(:)    !elements for output

     TYPE (ele02), POINTER :: head,tail
     TYPE (ele02_set), POINTER :: next
   END TYPE ele02_set

   TYPE (ele02_set), POINTER :: head, tail

 CONTAINS

   SUBROUTINE ini_ele02 (head, tail)
     !initialize a list of ELE02 sets

     TYPE (ele02_set), POINTER :: head, tail

     NULLIFY (head, tail)

   END SUBROUTINE ini_ele02

   SUBROUTINE add_ele02 (new, head, tail)
     !This subroutine adds data to the end of the list
     !Dummy arguments
     TYPE (ele02_set), POINTER :: new, head, tail

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
   END SUBROUTINE add_ele02

   SUBROUTINE srch_ele02 (head, anter, posic, name, found)
     !This subroutine searches for a set named "name"
     !Dummy arguments
     LOGICAL :: found
     CHARACTER (len=*) :: name ! set name
     TYPE (ele02_set), POINTER :: head, anter, posic

     found = .FALSE.                     !initializes flag
     NULLIFY (posic,anter)               !initializes pointers
     !Check if a list is empty
     IF (ASSOCIATED (head)) THEN         !if there are sets
       posic => head                     !point to first set
       DO
         IF (TRIM(posic%sname) == TRIM(name)) THEN    !check name
           found = .TRUE.                !set flag to .TRUE.
           EXIT                          !O.K. exit search
         END IF
         IF (ASSOCIATED(posic%next) ) THEN   !there is a next set
           anter => posic                    !point anter to present
           posic => posic%next               !point present to next
         ELSE
           EXIT                              !list exhausted, Exit search
         END IF
       END DO
     END IF
     IF (.NOT.found) NULLIFY (posic,anter)   !set not found, null pointers
   END SUBROUTINE srch_ele02

   SUBROUTINE del_ele02 (head, anter, posic)

     !This subroutine deletes a set pointed with posic

     TYPE (ele02_set), POINTER :: head, anter, posic

     TYPE (ele02), POINTER :: ea,ep
     INTEGER (kind=4) :: iel

     IF (.NOT.ASSOCIATED (anter)) THEN  !if anter pointer is null => head
       head => posic%next               !point first to next
     ELSE
       anter%next => posic%next         !skip posic
     END IF

     ! deallocation of the set memory is next done
     NULLIFY( ea )                  !nullify previous element in list
     ep => posic%head               !point present element to first
     DO iel = 1,posic%nelem         !for each element in the set
       CALL del_ele02e (posic%head,posic%tail, ea, ep )  !deletes element
     END DO

     NULLIFY (posic,anter)          !point to nothing
   END SUBROUTINE del_ele02

   ! ******* functions for a list of elements in a set ********

   SUBROUTINE ini_ele02e (head, tail)
     !initialize a list of ELE02 elements

     TYPE (ele02), POINTER :: head, tail

     NULLIFY (head, tail)       !initializes first and last pointer

   END SUBROUTINE ini_ele02e

   SUBROUTINE add_ele02e (new, head, tail)
     !This subroutine adds data to the end of the list
     !Dummy arguments
     TYPE (ele02), POINTER :: new, head, tail

     !Check if a list is empty
     IF (.NOT. ASSOCIATED (head)) THEN       !list is empty, start it
       head => new                           !first element
       tail => new                           !last element
       NULLIFY (tail%next)                   !last poit to nothing

     ELSE                                    !add a segment to the list
       tail%next => new                      !point to the new element
       NULLIFY (new%next)                    !nothing beyond the last
       tail => new                           !new last element

     END IF
   END SUBROUTINE add_ele02e

   SUBROUTINE srch_ele02e (head, anter, posic, kelem, found)
     !This subroutine searches for an element labeled "kelem"
     !Dummy arguments
     LOGICAL :: found
     INTEGER (kind=4) :: kelem
     TYPE (ele02), POINTER :: head, anter, posic

     found = .FALSE.             !initializes flag and pointers
     NULLIFY (posic,anter)

     IF (ASSOCIATED (head)) THEN          !Check if a list is empty
       posic => head                      !begin at top
       DO
         IF(posic%numel == kelem) THEN    !check if label found
           found = .TRUE.                 !Found
           EXIT                           !element found, EXIT
         END IF
         IF (ASSOCIATED(posic%next) ) THEN    !check if there are more els
           anter => posic                     !remember previous element
           posic => posic%next                !new present element
         ELSE
           EXIT                               !no more elements EXIT
         END IF
       END DO
     END IF
     IF (.NOT.found) NULLIFY (posic,anter)    !point to nothing
     RETURN
   END SUBROUTINE srch_ele02e

   SUBROUTINE del_ele02e (head, tail, anter, posic)

     !This subroutine deletes element pointed with posic

     TYPE (ele02), POINTER :: head, tail, anter, posic
     TYPE (ele02), POINTER :: e

     IF (.NOT.ASSOCIATED (anter)) THEN    !
       head => posic%next
     ELSE
       anter%next => posic%next
     END IF
     e => posic%next                       !keep pointer to next element
     IF( .NOT.ASSOCIATED(e) )tail => anter !last element in list
     DEALLOCATE (posic%gausv)              !deallocate variable arrays
     DEALLOCATE (posic)                    !deallocate fixed space
     posic => e                            !point to next element
     ! NULLIFY (posic,anter)
     RETURN
   END SUBROUTINE del_ele02e

   SUBROUTINE cut_ele02e (head, anter, posic)
     !This subroutine deletes a set pointed with posic
     ! without nullifying anter    ???? what for ????
     TYPE (ele02), POINTER :: head, anter, posic

     IF (.NOT.ASSOCIATED (anter)) THEN
       head => posic%next
     ELSE
       anter%next => posic%next
     ENDIF
     NULLIFY (posic)
   END SUBROUTINE cut_ele02e

   INCLUDE 'acvdf2.fi'
   INCLUDE 'commv2.fi'
   INCLUDE 'deltc2.fi'
   INCLUDE 'dumpi2.fi'
   INCLUDE 'elmda2.fi'
   INCLUDE 'gauss2.fi'
   !INCLUDE 'loadp2.fi'
   INCLUDE 'lumas2.fi'
   INCLUDE 'masel2.fi'
   INCLUDE 'outdy2.fi'
   INCLUDE 'resta2.fi'
   INCLUDE 'resvp2.fi'
   INCLUDE 'stra02.fi'
   INCLUDE 'stre02.fi'
   INCLUDE 'updlo2.fi'

 END MODULE ele02_db
