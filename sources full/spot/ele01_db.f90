 MODULE ele01_db
   USE param_db,ONLY : mnam, midn
   USE mat_dba, ONLY : sect_search,section,psecs,curve,cur_point,inte_cr
   USE c_input !, ONLY : lures,param,listen,exists,runend,openfi,getint,rdreqs,nnpar,get_name
   USE ele13_db, ONLY : search_ele13
   IMPLICIT NONE

   SAVE

   ! Derived type for a SPOT element
   TYPE ele01
     INTEGER (kind=4) :: numel  ! label of element
     INTEGER (kind=4) :: matno  ! Section number
     INTEGER (kind=4), POINTER :: lnods(:)  ! Conectivities
     LOGICAL :: rot             ! rotation stiffness considered
     REAL (kind=8), POINTER    :: gausv(:)      !gaussian variables
     !                            1-  : internal force
     ! for 1-node element  (ndime) x0 : original coordinates
     !                                : Internal moment
     !                     (ndime x ndime) L0 : original local system
     TYPE (ele01), POINTER :: next              !pointer to next element
   END TYPE ele01
   ! Derived type for a set of SPOT elements
   TYPE ele01_set
     CHARACTER (len=mnam) :: sname ! set name
     INTEGER (kind=4) :: nnode, &  ! number of nodes
                         nelem, &  ! number of elements
                         nreqs, &  ! number of GP for hist. output
                         narch     ! number of output unit
     LOGICAL :: gauss
     INTEGER (kind=4), POINTER :: ngrqs(:)    !elements for output
     CHARACTER (len=mnam) :: suname,slname ! upper and lower surface names of sets
     TYPE (ele01), POINTER :: head,tail
     TYPE (ele01_set), POINTER :: next
   END TYPE ele01_set

   TYPE (ele01_set), POINTER :: head, tail

 CONTAINS

   SUBROUTINE ini_ele01 (head, tail)
     !initialize a list of ELE01 sets

     TYPE (ele01_set), POINTER :: head, tail

     NULLIFY (head, tail)

   END SUBROUTINE ini_ele01

   SUBROUTINE add_ele01 (new, head, tail)
     !This subroutine adds data to the end of the list
     !Dummy arguments
     TYPE (ele01_set), POINTER :: new, head, tail

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
   END SUBROUTINE add_ele01

   SUBROUTINE srch_ele01 (head, anter, posic, name, found)
     !This subroutine searches for a set named "name"
     !Dummy arguments
     LOGICAL :: found
     CHARACTER (len=*) :: name ! set name
     TYPE (ele01_set), POINTER :: head, anter, posic

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
   END SUBROUTINE srch_ele01

   SUBROUTINE del_ele01 (head, anter, posic)

     !This subroutine deletes a set pointed with posic

     TYPE (ele01_set), POINTER :: head, anter, posic

     TYPE (ele01), POINTER :: ea,ep
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
       CALL del_ele01e (posic%head,posic%tail, ea, ep )  !deletes element
     END DO

     NULLIFY (posic,anter)          !point to nothing
   END SUBROUTINE del_ele01

   ! ******* functions for a list of elements in a set ********

   SUBROUTINE ini_ele01e (head, tail)
     !initialize a list of ELE01 elements

     TYPE (ele01), POINTER :: head, tail

     NULLIFY (head, tail)       !initializes first and last pointer

   END SUBROUTINE ini_ele01e

   SUBROUTINE add_ele01e (new, head, tail)
     !This subroutine adds data to the end of the list
     !Dummy arguments
     TYPE (ele01), POINTER :: new, head, tail

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
   END SUBROUTINE add_ele01e

   SUBROUTINE srch_ele01e (head, anter, posic, kelem, found)
     !This subroutine searches for an element labeled "kelem"
     !Dummy arguments
     LOGICAL :: found
     INTEGER (kind=4) :: kelem
     TYPE (ele01), POINTER :: head, anter, posic

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
   END SUBROUTINE srch_ele01e

   SUBROUTINE del_ele01e (head, tail, anter, posic)

     !This subroutine deletes element pointed with posic

     TYPE (ele01), POINTER :: head, tail, anter, posic
     TYPE (ele01), POINTER :: e

     IF (.NOT.ASSOCIATED (anter)) THEN    !
       head => posic%next
     ELSE
       anter%next => posic%next
     END IF
     e => posic%next                       !keep pointer to next element
     IF( .NOT.ASSOCIATED(e) )tail => anter !last element in list
     DEALLOCATE (posic%gausv)              !deallocate variable arrays
     DEALLOCATE (posic%lnods)              !deallocate array of connectivities
     DEALLOCATE (posic)                    !deallocate fixed space
     posic => e                            !point to next element
     ! NULLIFY (posic,anter)
     RETURN
   END SUBROUTINE del_ele01e

   SUBROUTINE cut_ele01e (head, anter, posic)
     !This subroutine deletes a set pointed with posic
     ! without nullifying anter    ???? what for ????
     TYPE (ele01), POINTER :: head, anter, posic

     IF (.NOT.ASSOCIATED (anter)) THEN
       head => posic%next
     ELSE
       anter%next => posic%next
     ENDIF
     NULLIFY (posic)
   END SUBROUTINE cut_ele01e

   INCLUDE 'acvdf1.fi'
   INCLUDE 'commv1.fi'
   INCLUDE 'deltc1.fi'
   INCLUDE 'dumpi1.fi'
   INCLUDE 'elmda1.fi'
   INCLUDE 'gauss1.fi'
   !INCLUDE 'loadp1.fi'
   INCLUDE 'lumas1.fi'
   INCLUDE 'masel1.fi'
   INCLUDE 'outdy1.fi'
   INCLUDE 'resta1.fi'
   INCLUDE 'resvp1.fi'
   INCLUDE 'updlo1.fi'

 END MODULE ele01_db
