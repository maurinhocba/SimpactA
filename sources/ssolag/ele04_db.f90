 MODULE ele04_db
   USE param_db,ONLY: mnam,midn,mlin
   USE mat_dba, ONLY: section,psecs,pmats,sect_search,mater,postv,pmater,snn
   USE c_input  !All?
  IMPLICIT NONE

   ! Derived type for an ELE04 element
   ! Total Lagrangial (Log strain) 3-D hexahedrañ Solid-Shell element
   INTEGER, PARAMETER :: nvare = 14, &!number of internal variables per integration point
                         nstre =  6, &! number of stress measures
                         nnode =  8, &! number of nodes per element
                         ngaud =  3   ! number of Gauss point to define volume
   TYPE ele04
     INTEGER (kind=4) :: numel  ! label of element
     INTEGER (kind=4) :: matno  ! Material number
     INTEGER (kind=4) :: lnods(nnode)      ! Conectivities
     REAL (kind=8) :: angle,            &  ! Euler angle between local t1-t2 and orthotropic
                      dvol(3),          &  ! Initial volume
                      cartd(nnode),     &  ! cartesian der (y3) of Shape Functions at center (EAS)
                      cdq(4,2,4,2),     &  ! cartesian derivatives at mid side points
                      nfdas(nnode,4,2), &  ! Nodal Function y3 Derivatives at the Assumed Strain points
                      area(4,2),        &  ! jacobians at mid side points and its sum
                      jacin(2,2,2),     &  ! in-plane inverse jacobian for shear ANS at face centers
                      alpha,ehb(2),fac(2)  ! EAS variable
     REAL (kind=8), POINTER :: stint(:,:),    & !Actual Kirchhoff stresses
                               gausv(:,:)       !Gauss-point internal variables
     TYPE (ele04), POINTER :: next              !pointer to next element
   END TYPE ele04

   ! Derived type for a set of ELE04 elements
   TYPE ele04_set
     CHARACTER (len=mnam) :: sname ! set name
     INTEGER (kind=4) :: nelem, &  ! number of elements
                         nreqs, &  ! number of GP for History output
                         narch, &  ! number of output unit
                         ngaus     !number of integration points
     LOGICAL :: gauss           ! .FALSE. -> Initial constants not
                                !  defined or not updated
     LOGICAL :: small           ! .TRUE. -> use Green strain tensor
                                ! .FALSE. -> Use log strains (Default)
     INTEGER :: plstr           ! compute Plastic Strain Flag
         ! -1 from Cauchy stress  0 - do not   1 from 2nd P-K
     REAL (kind=8) :: angdf     ! Default Euler angles between
                                ! Global X-Y-Z and orthotropic system
     REAL (kind=8) :: beta(3)   ! stabilization factors
     REAL (kind=8) :: psg(2,2)  ! factors for output
     INTEGER(kind=4) :: isg(2,2)! points for output
     INTEGER :: locax           ! local x definition option
     REAL (kind=8) :: btscal    ! critical time scale factor

     TYPE (ele04), POINTER    :: head, tail !pointer to first and last elm.
     INTEGER (kind=4), POINTER :: ngrqs(:)  !gauss points for output
     TYPE (ele04_set), POINTER :: next      !pointer to next set
   END TYPE ele04_set
   TYPE (ele04_set), POINTER, SAVE :: head,tail !first and last elements sets

 CONTAINS

   !----------- Set managment routines

   SUBROUTINE ini_ele04 (head, tail)
     !initialize a list of ELE04 sets

     !Dummy arguments
     TYPE (ele04_set), POINTER :: head, tail

     NULLIFY (head, tail)

   END SUBROUTINE ini_ele04

   SUBROUTINE add_ele04 (new, head, tail)
     !This subroutine adds a SET to the end of the list

     !Dummy arguments
     TYPE (ele04_set), POINTER :: new, head, tail

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
   END SUBROUTINE add_ele04

   SUBROUTINE srch_ele04 (head, anter, posic, name, found)
     !This subroutine searches for a set named "name"

     !Dummy arguments
     LOGICAL :: found
     CHARACTER (len=*) :: name ! set name
     TYPE (ele04_set), POINTER :: head, anter, posic

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
   END SUBROUTINE srch_ele04

   SUBROUTINE del_ele04 (head, anter, posic)

     !This subroutine deletes a set pointed with posic

     !Dummy arguments
     TYPE (ele04_set), POINTER :: head, anter, posic

     !local variables
     TYPE (ele04), POINTER :: ea,ep
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
       CALL del_ele04e (posic%head,posic%tail, ea, ep )  !deletes element
     END DO

     NULLIFY (posic,anter)          !point to nothing
   END SUBROUTINE del_ele04

   !----------- Element management routines

   SUBROUTINE ini_ele04e (head, tail)
     !initialize a list of ELE04 elements

     ! dummy arguments
     TYPE (ele04), POINTER :: head, tail

     NULLIFY (head, tail)       !initializes first and last pointer

   END SUBROUTINE ini_ele04e

   SUBROUTINE new_ele04e(elm)
   !Create a new element of ELE04 sets

     TYPE(ele04),POINTER:: elm

     ALLOCATE(elm)
     elm%numel = 0        !Initialize label of element
     elm%matno = 0        !     "     material number
     elm%angle = 0d0      !     "     angle between dir 1 and orthotropic dir 1
     elm%alpha = 0d0      !Initializes EAS parameter
     elm%lnods = 0
     elm%cdq   = 0d0
     elm%ehb   = 0d0
     elm%fac   = 1d0
     NULLIFY(elm%stint,elm%gausv)
     NULLIFY(elm%next)

   RETURN
   END SUBROUTINE new_ele04e

   SUBROUTINE add_ele04e (new, head, tail)
     !This subroutine adds data to the end of the list

     !Dummy arguments
     TYPE (ele04), POINTER :: new, head, tail

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
   END SUBROUTINE add_ele04e

   SUBROUTINE srch_ele04e (head, anter, posic, kelem, found)
     !This subroutine searches for an element labeled "kelem"

     !Dummy arguments
     LOGICAL :: found
     INTEGER (kind=4) :: kelem
     TYPE (ele04), POINTER :: head, anter, posic

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
   END SUBROUTINE srch_ele04e

   SUBROUTINE del_ele04e (head, tail, anter, posic)

     !This subroutine deletes element pointed with posic

     ! dummy arguments
     TYPE (ele04), POINTER :: head, tail, anter, posic
     ! local variables
     TYPE (ele04), POINTER :: e

     IF (.NOT.ASSOCIATED (anter)) THEN    !
       head => posic%next
     ELSE
       anter%next => posic%next
     END IF
     e => posic%next                       !keep pointer to next element
     IF( .NOT.ASSOCIATED(e) )tail => anter !last element in list
     DEALLOCATE (posic%stint)              !deallocate stress array
     DEALLOCATE (posic%gausv)              !deallocate variable arrays
     DEALLOCATE (posic)                    !deallocate fixed space
     posic => e                            !point to next element
     ! NULLIFY (posic,anter)
     RETURN
   END SUBROUTINE del_ele04e

   SUBROUTINE cut_ele04e (head, anter, posic)
     !This subroutine deletes a element pointed with posic
     ! without nullifying anter, DOES NOT deallocate memory

     ! dummy arguments
     TYPE (ele04), POINTER :: head, anter, posic

     IF (.NOT.ASSOCIATED (anter)) THEN
       head => posic%next
     ELSE
       anter%next => posic%next
     ENDIF
     NULLIFY (posic)
   END SUBROUTINE cut_ele04e

   INCLUDE 'acvd04.fi'
   INCLUDE 'bmat04.fi'
   INCLUDE 'bmat04q.fi'
   INCLUDE 'bsma04.fi'
   INCLUDE 'comm04.fi'
   INCLUDE 'delt04.fi'
   INCLUDE 'dump04.fi'
   INCLUDE 'elmd04.fi'
   INCLUDE 'expo04.fi'
   INCLUDE 'gaus04.fi'
   INCLUDE 'impo04.fi'
   INCLUDE 'jacob04.fi'
   INCLUDE 'lcsy04.fi'
   !INCLUDE 'load04.fi'
   INCLUDE 'luma04.fi'
   INCLUDE 'mase04.fi'
   INCLUDE 'outd04.fi'
   INCLUDE 'rest04.fi'
   INCLUDE 'resv04.fi'
   INCLUDE 'slno04.fi'
   INCLUDE 'slum04.fi'
   INCLUDE 'surf04.fi'
   INCLUDE 'updl04.fi'

 END MODULE ele04_db
