 MODULE ele18_db
   USE param_db,ONLY: mnam,midn,mlin
   USE mat_dba, ONLY: section,psecs,pmats,sect_search,mater,postv,pmater,snn
   USE c_input
   IMPLICIT NONE
     REAL(kind=8), PARAMETER :: sq3i =  0.577350269189626D+00

   ! Derived type for an ELE18 element
   ! Total Lagrangial (Log strain) 3-D brick element
   INTEGER(kind=4), PARAMETER :: nnode = 8    !number of nodes per element

   TYPE ele18
     INTEGER (kind=4) :: numel  ! label of element
     INTEGER (kind=4) :: matno  ! Material number
     INTEGER (kind=4) :: lnods(nnode)  ! Conectivities
     REAL (kind=8) :: angle(3)                  ! Euler angle between XYZ and orthotropic
     REAL (kind=8), POINTER :: dvol(:)          ! Initial volume
     REAL (kind=8), POINTER :: cartd(:,:,:)     !cartesian derivatives of SF
     REAL (kind=8), POINTER :: stint(:,:)       ! Actual Kirchhoff stresses
     REAL (kind=8), POINTER :: gausv(:,:)       !Gauss-point internal variables
     REAL (kind=8), POINTER :: nfdas(:,:),   &  ! Nodal Function Derivatives at the Assumed Strain points
                               jacin(:,:,:)     ! in-plane inverse jacobian
     TYPE (ele18), POINTER :: next              !pointer to next element
   END TYPE ele18

   TYPE ele18p
     TYPE(ele18), POINTER :: p
   END TYPE ele18p

   ! Derived type for a set of ELE18 elements
   TYPE ele18_set
     CHARACTER (len=mnam) :: sname ! set name
     INTEGER (kind=4) :: nelem, &  ! number of elements
                         nreqs, &  ! number of GP for hist. output
                         narch, &  ! number of output unit
                         ngaus     !number of integration points
     LOGICAL :: gauss           ! .FALSE. -> Initial constants not
                                !  defined or not updated
     LOGICAL :: small           ! .TRUE. -> use Green strain tensor
                                ! .FALSE. -> Use log strains (Default)
     INTEGER :: plstr           ! compute Plastic Strain Flag
         ! -1 from Cauchy stress  0 - do not   1 from 2nd P-K
     REAL (kind=8) :: angdf(3)    ! Default Euler angles between
                                  ! Global X-Y-Z and orthotropic system
     REAL(kind=8), POINTER :: gpc(:,:) !Gauss-points coordinates
     LOGICAL :: shell           ! .TRUE. -> modify transverse shear
                                ! .FALSE. -> use standard formulation
     LOGICAL :: bbar            ! .TRUE. -> use 1 GP for volume straian
                                ! .FALSE. -> use standard formulation
     LOGICAL :: check           ! .TRUE. -> modify transverse shear
     INTEGER :: locax           ! local x definition option
     REAL (kind=8) :: btscal      ! critical time scale factor

     TYPE (ele18), POINTER    :: head, tail !pointer to first and last elm.
     INTEGER (kind=4), POINTER :: ngrqs(:)  !gauss points for output
     TYPE (ele18_set), POINTER :: next      !pointer to next set
   END TYPE ele18_set
   TYPE (ele18_set), POINTER, SAVE :: head,tail !first and last elements sets

 CONTAINS
   SUBROUTINE ini_ele18 (head, tail)
     !initialize a list of ELE18 sets

     TYPE (ele18_set), POINTER :: head, tail

     NULLIFY (head, tail)

   END SUBROUTINE ini_ele18

   SUBROUTINE add_ele18 (new, head, tail)
     !This subroutine adds data to the end of the list
     !Dummy arguments
     TYPE (ele18_set), POINTER :: new, head, tail

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
   END SUBROUTINE add_ele18

   SUBROUTINE srch_ele18 (head, anter, posic, name, found)
     !This subroutine searches for a set named "name"
     !Dummy arguments
     LOGICAL :: found
     CHARACTER (len=*) :: name ! set name
     TYPE (ele18_set), POINTER :: head, anter, posic

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
   END SUBROUTINE srch_ele18

   SUBROUTINE del_ele18 (head, anter, posic)

     !This subroutine deletes a set pointed with posic

     TYPE (ele18_set), POINTER :: head, anter, posic

     TYPE (ele18), POINTER :: ea,ep
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
       CALL del_ele18e (posic%head,posic%tail, ea, ep )  !deletes element
     END DO

     NULLIFY (posic,anter)          !point to nothing
   END SUBROUTINE del_ele18

   ! ******* functions for a list of elements in a set ********

   SUBROUTINE ini_ele18e (head, tail)
     !initialize a list of ELE18 elements

     TYPE (ele18), POINTER :: head, tail

     NULLIFY (head, tail)       !initializes first and last pointer

   END SUBROUTINE ini_ele18e

   SUBROUTINE new_ele18e(elm)
   !Create a new element of ELE18 sets

     TYPE(ele18),POINTER:: elm

     ALLOCATE(elm)
     elm%numel = 0        !Initialize label of element
     elm%matno = 0        !     "     material number
     elm%lnods(1:8) = 0   !     "     conectivities
     elm%angle = 0d0      !     "     angle between dir 1 and orthotropic dir 1
     NULLIFY(elm%dvol,elm%cartd,elm%stint,elm%gausv,elm%nfdas,elm%jacin)
     NULLIFY(elm%next)

   RETURN
   END SUBROUTINE new_ele18e

   SUBROUTINE add_ele18e (new, head, tail)
     !This subroutine adds data to the end of the list
     !Dummy arguments
     TYPE (ele18), POINTER :: new, head, tail

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
   END SUBROUTINE add_ele18e

   SUBROUTINE srch_ele18e (head, anter, posic, kelem, found)
     !This subroutine searches for an element labeled "kelem"
     !Dummy arguments
     LOGICAL :: found
     INTEGER (kind=4) :: kelem
     TYPE (ele18), POINTER :: head, anter, posic

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
   END SUBROUTINE srch_ele18e

   SUBROUTINE del_ele18e (head, tail, anter, posic)

     !This subroutine deletes element pointed with posic

     TYPE (ele18), POINTER :: head, tail, anter, posic
     TYPE (ele18), POINTER :: e

     IF (.NOT.ASSOCIATED (anter)) THEN    !
       head => posic%next
     ELSE
       anter%next => posic%next
     END IF
     e => posic%next                       !keep pointer to next element
     IF( .NOT.ASSOCIATED(e) )tail => anter !last element in list
     DEALLOCATE (posic%dvol)               !deallocate vol arrays
     DEALLOCATE (posic%stint)              !deallocate stress array
     IF(ASSOCIATED(posic%gausv))DEALLOCATE (posic%gausv)  !deallocate variable arrays
     DEALLOCATE (posic%cartd)              !deallocate fixed space
     IF(ASSOCIATED(posic%nfdas))DEALLOCATE (posic%nfdas,posic%jacin)  !deallocate variable arrays
     DEALLOCATE (posic)                    !deallocate fixed space
     posic => e                            !point to next element
     ! NULLIFY (posic,anter)
     RETURN
   END SUBROUTINE del_ele18e

   SUBROUTINE cut_ele18e (head, anter, posic)
     !This subroutine deletes a set pointed with posic
     ! without nullifying anter    ???? what for ????
     TYPE (ele18), POINTER :: head, anter, posic

     IF (.NOT.ASSOCIATED (anter)) THEN
       head => posic%next
     ELSE
       anter%next => posic%next
     ENDIF
     NULLIFY (posic)
   END SUBROUTINE cut_ele18e

   INCLUDE 'acvd18.fi'
   INCLUDE 'bmat18.fi'
   INCLUDE 'bsma18.fi'
   INCLUDE 'comm18.fi'
   INCLUDE 'delt18.fi'
   INCLUDE 'dump18.fi'
   INCLUDE 'eige18.fi'
   INCLUDE 'elmd18.fi'
   INCLUDE 'expo18.fi'
   INCLUDE 'gaus18.fi'
   INCLUDE 'impo18.fi'
   !INCLUDE 'jacobi.fi'
   INCLUDE 'lgst18.fi'
   INCLUDE 'lcas18.fi'
   !INCLUDE 'luma18.fi'
   INCLUDE 'masm18.fi'
   INCLUDE 'mase18.fi'
   INCLUDE 'nods18.fi'
   INCLUDE 'outd18.fi'
   INCLUDE 'rest18.fi'
   INCLUDE 'resv18.fi'
   INCLUDE 'rubber3d.fi'
   INCLUDE 'rubber3dn.fi'
   INCLUDE 'slum18.fi'
   INCLUDE 'surf18.fi'
   INCLUDE 'secd18.fi'
   INCLUDE 'slno18.fi'
   INCLUDE 'updl18.fi'

 END MODULE ele18_db
