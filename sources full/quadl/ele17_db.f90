 MODULE ele17_db
   USE param_db,ONLY: mnam,midn
   USE mat_dba, ONLY: section,sect_search,mater,psecs,pmats
   USE c_input
   USE meshmo_db,ONLY: numpo, nodset
   IMPLICIT NONE

   ! Derived type for an ELE17 element
   ! Total Lagrangial (Log strain) 2-D quad element
   INTEGER, PARAMETER :: nnode = 4    !number of nodes per element
   TYPE ele17
     INTEGER (kind=4) :: numel  ! label of element
     INTEGER (kind=4) :: matno  ! Material number
     INTEGER (kind=4) :: lnods(nnode)  ! Conectivities
     REAL (kind=8) :: angle     ! angle between dir 1 and orthotropic dir 1
     REAL (kind=8), POINTER :: dvol(:),      & ! (NGAUS) Initial area
                               cartd(:,:,:), & !(nnode,2,ngaus) cartesian derivatives of SF
                               stint(:,:),   & !(4,ngaus) Actual Integrated forces
                               gausv(:,:)      !Gauss-point values
     TYPE (ele17), POINTER :: next              !pointer to next element
   END TYPE ele17

   ! Derived type for a set of ELE17 elements
   TYPE ele17_set
     CHARACTER (len=mnam) :: sname ! set name
     INTEGER (kind=4) :: nelem, &  ! number of elements
                         nreqs, &  ! number of GP for hist. output
                         narch, &  ! number of output unit
                         ngaus     !number of integration points
     LOGICAL :: gauss           ! .FALSE. -> Initial constants not
                                !  defined or not updated
     INTEGER :: plstr           ! compute Plastic Strain Flag
         ! -1 from Cauchy stress  0 - do not   1 from 2nd P-K
     REAL (kind=8) :: angdf, & ! angle between X_1 and orthotropic dir 1
                      stabs    ! stabilization coefficient

     TYPE (ele17), POINTER    :: head, tail !pointer to first and last elm.
     INTEGER (kind=4), POINTER :: ngrqs(:)  !gauss points for output
     TYPE (ele17_set), POINTER :: next      !pointer to next set
   END TYPE ele17_set
   TYPE (ele17_set), POINTER, SAVE :: head,tail !first and last elements sets

 CONTAINS
   SUBROUTINE ini_ele17 (head, tail)
     !initialize a list of ELE17 sets

     TYPE (ele17_set), POINTER :: head, tail

     NULLIFY (head, tail)

   END SUBROUTINE ini_ele17

   SUBROUTINE add_ele17 (new, head, tail)
     !This subroutine adds data to the end of the list
     !Dummy arguments
     TYPE (ele17_set), POINTER :: new, head, tail

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
   END SUBROUTINE add_ele17

   SUBROUTINE srch_ele17 (head, anter, posic, name, found)
     !This subroutine searches for a set named "name"
     !Dummy arguments
     LOGICAL :: found
     CHARACTER (len=*) :: name ! set name
     TYPE (ele17_set), POINTER :: head, anter, posic

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
   END SUBROUTINE srch_ele17

   SUBROUTINE del_ele17 (head, anter, posic)

     !This subroutine deletes a set pointed with posic

     TYPE (ele17_set), POINTER :: head, anter, posic

     TYPE (ele17), POINTER :: ea,ep
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
       CALL del_ele17e (posic%head,posic%tail, ea, ep )  !deletes element
     END DO

     NULLIFY (posic) !,anter)          !point to nothing
   END SUBROUTINE del_ele17

   ! ******* functions for a list of elements in a set ********

   SUBROUTINE ini_ele17e (head, tail)
     !initialize a list of ELE17 elements

     TYPE (ele17), POINTER :: head, tail

     NULLIFY (head, tail)       !initializes first and last pointer

   END SUBROUTINE ini_ele17e

   SUBROUTINE add_ele17e (new, head, tail)
     !This subroutine adds data to the end of the list
     !Dummy arguments
     TYPE (ele17), POINTER :: new, head, tail

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
   END SUBROUTINE add_ele17e

   SUBROUTINE srch_ele17e (head, anter, posic, kelem, found)
     !This subroutine searches for an element labeled "kelem"
     !Dummy arguments
     LOGICAL :: found
     INTEGER (kind=4) :: kelem
     TYPE (ele17), POINTER :: head, anter, posic

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
   END SUBROUTINE srch_ele17e

   SUBROUTINE del_ele17e (head, tail, anter, posic)

     !This subroutine deletes element pointed with posic

     TYPE (ele17), POINTER :: head, tail, anter, posic
     TYPE (ele17), POINTER :: e

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
   END SUBROUTINE del_ele17e

   SUBROUTINE cut_ele17e (head, anter, posic)
     !This subroutine deletes a set pointed with posic
     ! without nullifying anter    ???? what for ????
     TYPE (ele17), POINTER :: head, anter, posic

     IF (.NOT.ASSOCIATED (anter)) THEN
       head => posic%next
     ELSE
       anter%next => posic%next
     ENDIF
     NULLIFY (posic)
   END SUBROUTINE cut_ele17e

   INCLUDE 'acvd17.fi'
   INCLUDE 'bmat17.fi'
   INCLUDE 'comm17.fi'
   INCLUDE 'delt17.fi'
   INCLUDE 'dump17.fi'
   INCLUDE 'elmd17.fi'
   INCLUDE 'gaus17.fi'
   INCLUDE 'jaco17.fi'
   !INCLUDE 'load17.fi'
   !INCLUDE 'luma17.fi'
   INCLUDE 'mase17.fi'
   INCLUDE 'masm17.fi'
   INCLUDE 'nods17.fi'
   INCLUDE 'outd17.fi'
   INCLUDE 'rest17.fi'
   INCLUDE 'resv17.fi'
   INCLUDE 'resv17r.fi'
   INCLUDE 'slnods17.fi'
   INCLUDE 'surf17.fi'
   INCLUDE 'updl17.fi'


 END MODULE ele17_db
