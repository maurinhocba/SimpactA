 MODULE ele05_db
   USE param_db,ONLY: mnam,midn,mlin
   USE mat_dba, ONLY: section,psecs,pmats,sect_search,mater,postv,pmater,snn
   USE c_input  !All?
  IMPLICIT NONE

   ! Derived type for an ELE05 element
   ! Total Lagrangial (Log strain) 3-D prism Solid-Shell element
   INTEGER, PARAMETER :: nstre =  6, &! number of stress measures
                         nnb   =  6, &! number of basic nodes per element
                         ngaud =  3   ! number of Gauss point to define volume
   TYPE ele05
     INTEGER (kind=4) :: numel  ! label of element
     INTEGER (kind=4) :: matno  ! Material number
     INTEGER, POINTER  :: lnods(:)  ! Conectivities
     REAL (kind=8) :: angle,            &  ! Euler angle between local t1-t2 and orthotropic
                      dvol(3),          &  ! Initial volume
                      cartd(nnb),       &  ! cartesian der (y3) of Shape Functions at center (EAS)
                      nfdas(nnb,3,2),   &  ! Nodal Function Derivatives at the Assumed Strain points
                      jacin(2,2,2),     &  ! in-plane inverse jacobian for shear ANS
                      alpha                ! EAS variable
     REAL (kind=8), POINTER :: cdq(:,:,:,:),  & ! cartesian derivatives for QUAD approach (AS)
                               stint(:,:),    & !Actual Kirchhoff stresses
                               gausv(:,:)       !Gauss-point internal variables
     TYPE (ele05), POINTER :: next              !pointer to next element
   END TYPE ele05

   ! Derived type for a set of ELE05 elements
   TYPE ele05_set
     CHARACTER (len=mnam) :: sname ! set name
     INTEGER (kind=4) :: nelem, &  ! number of elements
                         nnode, &  ! number of nodes per element
                         nreqs, &  ! number of GP for History output
                         narch, &  ! number of output unit
                         ngaus     !number of integration points
     LOGICAL :: gauss           ! .FALSE. -> Initial constants not
                                !  defined or not updated
     LOGICAL :: small           ! .TRUE. -> use Green strain tensor
                                ! .FALSE. -> Use log strains (Default)
     LOGICAL :: lface           ! .TRUE. -> if extended connectivities have been computed
                                ! .FALSE.
     INTEGER :: plstr           ! compute Plastic Strain Flag
         ! -1 from Cauchy stress  0 - do not   1 from 2nd P-K
     REAL (kind=8) :: angdf            ! Default Euler angles between
                                  ! Global X-Y-Z and orthotropic system
     REAL (kind=8) :: psg(2,2)  ! factors for output
     INTEGER(kind=4) :: isg(2,2)! points for output
     LOGICAL :: quad            ! .TRUE. -> use in-plane quadratic approach
                                ! .FALSE. -> use standard formulation
     INTEGER :: locax           ! local x definition option
     REAL (kind=8) :: btscal    ! critical time scale factor

     TYPE (ele05), POINTER    :: head, tail !pointer to first and last elm.
     INTEGER (kind=4), POINTER :: ngrqs(:)  !gauss points for output
     TYPE (ele05_set), POINTER :: next      !pointer to next set
   END TYPE ele05_set
   TYPE (ele05_set), POINTER, SAVE :: head,tail !first and last elements sets

 CONTAINS

   !----------- Set managment routines

   SUBROUTINE ini_ele05 (head, tail)
     !initialize a list of ELE05 sets

     !Dummy arguments
     TYPE (ele05_set), POINTER :: head, tail

     NULLIFY (head, tail)

   END SUBROUTINE ini_ele05

   SUBROUTINE add_ele05 (new, head, tail)
     !This subroutine adds a SET to the end of the list

     !Dummy arguments
     TYPE (ele05_set), POINTER :: new, head, tail

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

     END IF
   END SUBROUTINE add_ele05

   SUBROUTINE srch_ele05 (head, anter, posic, name, found)
     !This subroutine searches for a set named "name"

     !Dummy arguments
     LOGICAL :: found
     CHARACTER (len=*) :: name ! set name
     TYPE (ele05_set), POINTER :: head, anter, posic

     found = .FALSE.                     !initializes flag
     NULLIFY (posic,anter)               !initializes pointers
     !Check if a list is empty
     IF (ASSOCIATED (head)) THEN         !if there are sets
       posic => head                     !point to first set
       DO
         IF(TRIM(posic%sname) == TRIM(name)) THEN    !check name
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
   END SUBROUTINE srch_ele05

   SUBROUTINE del_ele05 (head, anter, posic)

     !This subroutine deletes a set pointed with posic

     !Dummy arguments
     TYPE (ele05_set), POINTER :: head, anter, posic

     !local variables
     TYPE (ele05), POINTER :: ea,ep
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
       CALL del_ele05e (posic%head,posic%tail, ea, ep )  !deletes element
     END DO

     NULLIFY (posic,anter)          !point to nothing
   END SUBROUTINE del_ele05

   !----------- Element management routines

   SUBROUTINE ini_ele05e (head, tail)
     !initialize a list of ELE05 elements

     ! dummy arguments
     TYPE (ele05), POINTER :: head, tail

     NULLIFY (head, tail)       !initializes first and last pointer

   END SUBROUTINE ini_ele05e

   SUBROUTINE new_ele05e(elm)
   !Create a new element of ELE05 sets

     TYPE(ele05),POINTER:: elm

     ALLOCATE(elm)
     elm%numel = 0        !Initialize label of element
     elm%matno = 0        !     "     material number
     elm%angle = 0d0      !     "     angle between dir 1 and orthotropic dir 1
     elm%alpha = 0d0      !Initializes EAS parameter
     NULLIFY(elm%lnods,elm%cdq,elm%stint,elm%gausv)
     NULLIFY(elm%next)

   RETURN
   END SUBROUTINE new_ele05e

   SUBROUTINE add_ele05e (new, head, tail)
     !This subroutine adds data to the end of the list

     !Dummy arguments
     TYPE (ele05), POINTER :: new, head, tail

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
   END SUBROUTINE add_ele05e

   SUBROUTINE srch_ele05e (head, anter, posic, kelem, found)
     !This subroutine searches for an element labeled "kelem"

     !Dummy arguments
     LOGICAL :: found
     INTEGER (kind=4) :: kelem
     TYPE (ele05), POINTER :: head, anter, posic

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
   END SUBROUTINE srch_ele05e

   SUBROUTINE del_ele05e (head, tail, anter, posic)

     !This subroutine deletes element pointed with posic

     ! dummy arguments
     TYPE (ele05), POINTER :: head, tail, anter, posic
     ! local variables
     TYPE (ele05), POINTER :: e

     IF (.NOT.ASSOCIATED (anter)) THEN    !
       head => posic%next
     ELSE
       anter%next => posic%next
     END IF
     e => posic%next                       !keep pointer to next element
     IF( .NOT.ASSOCIATED(e) )tail => anter !last element in list
     DEALLOCATE (posic%lnods)              !deallocate connectivities array
     DEALLOCATE (posic%cdq)                !deallocate in-plane cartesyan derivatives array
     DEALLOCATE (posic%stint)              !deallocate stress array
     DEALLOCATE (posic%gausv)              !deallocate variable arrays
     DEALLOCATE (posic)                    !deallocate fixed space
     posic => e                            !point to next element
     ! NULLIFY (posic,anter)
     RETURN
   END SUBROUTINE del_ele05e

   SUBROUTINE cut_ele05e (head, anter, posic)
     !This subroutine deletes a element pointed with posic
     ! without nullifying anter, DOES NOT deallocate memory

     ! dummy arguments
     TYPE (ele05), POINTER :: head, anter, posic

     IF (.NOT.ASSOCIATED (anter)) THEN
       head => posic%next
     ELSE
       anter%next => posic%next
     ENDIF
     NULLIFY (posic)
   END SUBROUTINE cut_ele05e

   INCLUDE 'acvd05.fi'
   INCLUDE 'bmat05.fi'
   INCLUDE 'bmat05q.fi'
   INCLUDE 'bsma05.fi'
   INCLUDE 'comm05.fi'
   INCLUDE 'delt05.fi'
   INCLUDE 'dump05.fi'
   INCLUDE 'elmd05.fi'
   INCLUDE 'expo05.fi'
   INCLUDE 'gaus05.fi'
   INCLUDE 'impo05.fi'
   INCLUDE 'jacob05.fi'
   INCLUDE 'lcsy05.fi'
   !INCLUDE 'load05.fi'
   INCLUDE 'luma05.fi'
   INCLUDE 'mase05.fi'
   INCLUDE 'outd05.fi'
   INCLUDE 'rest05.fi'
   INCLUDE 'resv05.fi'
   INCLUDE 'slno05.fi'
   INCLUDE 'slum05.fi'
   INCLUDE 'surf05.fi'
   INCLUDE 'toar05.fi'
   INCLUDE 'updl05.fi'

 END MODULE ele05_db
