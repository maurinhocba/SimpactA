 MODULE ele08_db
   USE param_db,ONLY: mnam,midn,mlin
   USE mat_dba, ONLY: section,mater,sect_search,psecs !,pmats,postv
   USE c_input
   !USE kinc_db, ONLY : nndpd
   IMPLICIT NONE
   SAVE

   ! Derived type for a  BEAM  element
   TYPE ele08
     INTEGER (kind=4) :: numel  ! label of element
     INTEGER (kind=4) :: matno  ! Material number
     INTEGER (kind=4), POINTER :: lnods(:)  ! Conectivities
     REAL (kind=8), POINTER    :: llbd(:,:,:)     !(3,3,axesc) local system
     REAL (kind=8), POINTER    :: jac(:)        !(ngaus) jacobian at gauss points
     REAL (kind=8), POINTER    :: stran(:,:)    !(nstre,ngaus)
     REAL (kind=8), POINTER    :: stres(:,:)    !(nstre,ngaus) gaussian variables
     REAL (kind=8), POINTER    :: epdef(:)      !(ngaus) gaussian variables
     REAL (kind=8), POINTER    :: sedef(:)      !(ngaus)
     TYPE (ele08), POINTER :: next              !pointer to next element
   END TYPE ele08
   ! Derived type for a set of BEAM elements
   TYPE ele08_set
     CHARACTER(len=mnam):: sname   ! set name
     INTEGER (kind=4) :: nelem, & ! number of elements
                         nnode, & ! number of nodes per element
                         ngaus, & ! number of gauss points per element
                         axesc, & ! axes code for local systems
                         mtype, & ! material model type
                         nreqs, & ! number of GP for hist. output
                         narch, & ! number of output unit
                         nvare    ! length of gausv
     INTEGER (kind=4), POINTER ::  ngrqs(:)     !Gauss points for OUTPUT
     REAL (kind=8), POINTER :: posgp(:),   & !Gauss point positions
                               shape(:,:), & !Nodal shape function at Gauss-P
                               deriv(:,:), & !Nodal function deriv at Gauss-P
                               weigh(:)      !Gauss point weigths
     TYPE (ele08), POINTER :: head,tail
     TYPE (ele08_set), POINTER :: next
   END TYPE ele08_set
   TYPE (ele08_set), POINTER :: head
   TYPE (ele08_set), POINTER :: tail

 CONTAINS
   SUBROUTINE ini_ele08 (head, tail)
     !initialize a list of BEAM sets

     TYPE (ele08_set), POINTER :: head, tail

     NULLIFY (head, tail)

   END SUBROUTINE ini_ele08

   SUBROUTINE add_ele08 (new, head, tail)
     !This subroutine adds data to the end of the list
     !Dummy arguments
     TYPE (ele08_set), POINTER :: new, head, tail

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
   END SUBROUTINE add_ele08

   SUBROUTINE srch_ele08 (head, anter, posic, name, found)
     !This subroutine searches for a set named "name"
     !Dummy arguments
     LOGICAL :: found
     CHARACTER(len=*):: name ! set name
     TYPE (ele08_set), POINTER :: head, anter, posic

     found = .FALSE.
     NULLIFY (posic,anter)
     !Check if a list is empty
     IF (ASSOCIATED (head)) THEN
       posic => head
       DO
         IF(TRIM(posic%sname) == TRIM(name)) THEN
           found = .TRUE.
           EXIT
         END IF
         IF (ASSOCIATED(posic%next) ) THEN
           anter => posic
           posic => posic%next
         ELSE
           EXIT
         END IF
       END DO
     ENDIF
     IF (.NOT.found) NULLIFY (posic,anter)
   END SUBROUTINE srch_ele08

   SUBROUTINE del_ele08 (head, tail, anter, posic)
     !This subroutine deletes a set pointed with posic
     TYPE (ele08_set), POINTER :: head, tail, anter, posic

     IF (.NOT.ASSOCIATED (anter)) THEN
       head => posic%next
     ELSE
       anter%next => posic%next
     ENDIF
     ! if posic == tail
     IF (.NOT.ASSOCIATED (posic%next) ) tail => anter
     CALL dalloc08 (posic)
     NULLIFY (anter)
   END SUBROUTINE del_ele08

   SUBROUTINE dalloc08 (elset)
     ! deallocates a set
     TYPE (ele08_set), POINTER :: elset
     DEALLOCATE ( elset%ngrqs )
     DEALLOCATE ( elset%posgp, elset%shape, elset%deriv, &
                  elset%weigh )
     DEALLOCATE (elset)
   END SUBROUTINE dalloc08

   ! ******* functions for a list of elements in a set ********

   SUBROUTINE ini_ele08e (head, tail)
     !initialize a list of ELE08 elements

     TYPE (ele08), POINTER :: head, tail

     NULLIFY (head, tail)       !initializes first and last pointer

   END SUBROUTINE ini_ele08e

   SUBROUTINE add_ele08e (new, head, tail)
     !This subroutine adds data to the end of the list
     !Dummy arguments
     TYPE (ele08), POINTER :: new, head, tail

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
   END SUBROUTINE add_ele08e

   SUBROUTINE srch_ele08e (head, anter, posic, kelem, found)
     !This subroutine searches for an element labeled "kelem"
     !Dummy arguments
     LOGICAL :: found
     INTEGER (kind=4) :: kelem
     TYPE (ele08), POINTER :: head, anter, posic

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
   END SUBROUTINE srch_ele08e

   SUBROUTINE del_ele08e (head, tail, anter, posic)

     !This subroutine deletes element pointed with posic

     TYPE (ele08), POINTER :: head, tail, anter, posic
     TYPE (ele08), POINTER :: e

     IF (.NOT.ASSOCIATED (anter)) THEN    !
       head => posic%next
     ELSE
       anter%next => posic%next
     END IF
     e => posic%next                       !keep pointer to next element
     IF( .NOT.ASSOCIATED(e) )tail => anter !last element in list
     DEALLOCATE (posic%llbd,posic%jac,posic%stran,posic%stres,posic%epdef,posic%sedef)   !deallocate variable arrays
     DEALLOCATE (posic)                    !deallocate fixed space
     posic => e                            !point to next element
     ! NULLIFY (posic,anter)
     RETURN
   END SUBROUTINE del_ele08e

   SUBROUTINE cut_ele08e (head, anter, posic)
     !This subroutine deletes a set pointed with posic
     ! without nullifying anter    ???? what for ????
     TYPE (ele08), POINTER :: head, anter, posic

     IF (.NOT.ASSOCIATED (anter)) THEN
       head => posic%next
     ELSE
       anter%next => posic%next
     ENDIF
     NULLIFY (posic)
   END SUBROUTINE cut_ele08e

   INCLUDE 'acvdf8.fi'
   INCLUDE 'commv8.fi'
   INCLUDE 'deltc8.fi'
   INCLUDE 'dumpi8.fi'
   INCLUDE 'elmda8.fi'
   INCLUDE 'gauss8.fi'
   INCLUDE 'intrf8.fi'
   INCLUDE 'istg18.fi'
   INCLUDE 'istg28.fi'
   !INCLUDE 'loadp8.fi'
   INCLUDE 'locla8.fi'
   !INCLUDE 'lumas8.fi'
   INCLUDE 'masel8.fi'
   INCLUDE 'masmt8.fi'
   INCLUDE 'nodxy8.fi'
   INCLUDE 'outdy8.fi'
   INCLUDE 'rare18.fi'
   INCLUDE 'resta8.fi'
   INCLUDE 'resvp8.fi'
   INCLUDE 'setga8.fi'
   INCLUDE 'stran8.fi'
   INCLUDE 'surf08.fi'
   INCLUDE 'updlo8.fi'

 END MODULE ele08_db
