 MODULE ele25_db
   USE param_db,ONLY: mnam,midn,mlin
   USE mat_dba, ONLY: section,sect_search,psecs,pmats,mater,postv,pmater
   USE kinc_db, ONLY : nndpd
   USE c_input
   IMPLICIT NONE

   ! Derived type for an ELE25 element
   ! Basic thin Shell Quadrilateral (BSQ)
   ! including Non-smooth surfaces and branching

   ! Reference: C. Estrada & F. Flores
   !            Unpublished

   SAVE

   TYPE sideb
     INTEGER (kind=4) :: nn
     INTEGER (kind=4), POINTER :: lnods(:)  !-1:2*nn = connectivities
     REAL (kind=8), POINTER  :: alph0(:)    !nn   = initial angles + side length
     REAL (kind=8), POINTER  :: gamma(:)    !nn   = distorsions
     REAL (kind=8), POINTER  :: fc(:,:)     !nn   = factors
     REAL (kind=8), POINTER  :: t0(:,:)     !nn   = factors
     REAL (kind=8), POINTER  :: tb(:,:)     !nn   = factors
     REAL (kind=8), POINTER  :: c(:,:,:)    !nn   = normal derivatives
     REAL (kind=8), POINTER  :: bbr(:,:)    !ndof*(nn+2) = nodal B matrix
     TYPE (sideb), POINTER :: next
   END TYPE sideb

   TYPE pside
      TYPE (sideb), POINTER :: p
   END TYPE pside

   INTEGER (kind=4), PARAMETER ::  &
      ln(4,4) = RESHAPE(           &     !element connectivities
               (/ 2,1,5,6, 3,2,7,8, 4,3,9,10, 1,4,11,12 /), (/4,4/) )
   !   hh= side element connectivities
   INTEGER(kind=4), PARAMETER :: hh(3,4) = RESHAPE((/ 4,2,1, 5,3,2, 6,4,3, 6,1,4 /), (/3,4/) ), &
      fn(4) = (/5,7,9,11/),        &     !first node of each side
      nextn(7) = (/2,3,4,1,2,3,4/)       !cycling list

   TYPE ele25
     INTEGER (kind=4) :: numel  ! label of element
     INTEGER (kind=4) :: matno  ! Material number
     INTEGER (kind=4) :: lnods(12) ! Conectivities
     INTEGER (kind=4) :: lside(4)  ! side elements
     REAL (kind=8) :: area(0:4), & ! Initial area
                      lb,        & ! Thickness ratio
                      angle,     & ! angle between dir 1 and orthotropic dir 1
                      cartd(4,2,0:8),  & !cartesian derivatives
                      normd(4,2,1:8),  & !normal derivatives
                      ns(2,4),   & ! side normals
                      a0(4),     & ! initial angles with side elements
                      ci(4),     & ! coefficients for angle change
                      gamma(4),  & ! distorsion at each side
                      stra0(3),  & ! Original curvature tensor
                      stra1(6)     ! Actual mid-surface metric tensor
     REAL (kind=8) :: sf (1) !stabilization forces
     REAL (kind=8), POINTER :: gausv(:,:)       !layer values
     TYPE (pside), POINTER :: si(:)
     TYPE (ele25), POINTER :: next              !pointer to next element
   END TYPE ele25

   ! Derived type for a set of ELE25 elements
   TYPE ele25_set
     CHARACTER (len=mnam) :: sname ! set name
     INTEGER (kind=4) :: nelem  ! number of elements
     INTEGER (kind=4) :: nreqs  ! number of GP for hist. output
     INTEGER (kind=4) :: narch  ! number of output unit
     INTEGER (kind=4) :: nbs    ! number of branching sides
     LOGICAL :: logst           ! use logarithmic strain
     LOGICAL :: lside           ! .FALSE. -> topological arrays not
                                !  defined or not updated
     LOGICAL :: gauss           ! .FALSE. -> Initial constants not
                                !  defined or not updated
     INTEGER :: plstr           ! compute Plastic Strain Flag
         ! -1 from Cauchy stress  0 - do not   1 from 2nd P-K
     INTEGER :: locax           ! local x definition option
     REAL (kind=8) ::  angdf     ! angle between X_1 and orthotropic dir 1
     REAL (kind=8) ::  stabs, &  ! stabilization factor for membrane
                       stabb     ! stabilization factor for bending
     REAL (kind=8), POINTER :: stint(:,:)   !moments, forces and shears
     TYPE (ele25), POINTER    :: head, tail !pointer to first and last elm.
     TYPE (sideb), POINTER :: bhead , btail !pointers to branching data base
     INTEGER (kind=4), POINTER :: ngrqs(:)  !gauss points for output
     INTEGER (kind=4), POINTER :: sidel(:,:)  !side elements (temporary, no dumping)
     TYPE (ele25_set), POINTER :: next      !pointer to next set
   END TYPE ele25_set
   TYPE (ele25_set), POINTER, SAVE :: head,tail !first and last elements sets

 CONTAINS

   FUNCTION atan4(a,b,c,d)
    IMPLICIT NONE
    REAL(kind=8) :: atan4
    REAL(kind=8), INTENT(IN) :: a,b,c,d
    REAL (kind=8),PARAMETER :: twopi=6.283185307179586d0, &
                                  pi=3.141592553589793d0
      atan4 = ATAN2(a,b) - c
      !  limit angle change to Pi (180 degrees)
      IF( atan4-d > pi  )THEN
        atan4 = atan4 - twopi
      ELSE IF( atan4-d < -pi )THEN
        atan4 = atan4 + twopi
      END IF
    END FUNCTION

   SUBROUTINE new_ele25(elset)
   !Create a new element of ELE25 sets

     TYPE(ele25_set),POINTER:: elset

     ALLOCATE(elset)
     elset%sname = ''       !Initialize set name
     elset%nelem = 0        !     "     number of elements
     elset%nreqs = 0        !     "     number of GP for hist. output
     elset%narch = 0        !     "     number of output unit
     elset%nbs   = 0        !     "     number of branching sides
     elset%logst = .FALSE.  !     "     use logarithmic strain
     elset%lside = .FALSE.  !     "     flag to compute LSIDE
     elset%gauss = .FALSE.  !     "     flag to compute Gauss constants
     elset%plstr = 0        !     "     compute Plastic Strain Flag
     elset%locax = 3        !     "     local x definition
     elset%angdf = 0d0      !     "     angle between X_1 and orthotropic dir 1
     elset%stabs = 0d0      !     "     membrane stabilization factor
     elset%stabb = 0d0      !     "     bending stabilization factor
     NULLIFY(elset%head,elset%tail,elset%ngrqs,elset%stint,elset%bhead,elset%btail)
     NULLIFY(elset%next)

   RETURN
   END SUBROUTINE new_ele25

   SUBROUTINE ini_ele25 (head, tail)
     !initialize a list of ELE25 sets

     TYPE (ele25_set), POINTER :: head, tail

     NULLIFY (head, tail)

   END SUBROUTINE ini_ele25

   SUBROUTINE add_ele25 (new, head, tail)
     !This subroutine adds data to the end of the list
     !Dummy arguments
     TYPE (ele25_set), POINTER :: new, head, tail

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
   END SUBROUTINE add_ele25

   SUBROUTINE srch_ele25 (head, anter, posic, name, found)
     !This subroutine searches for a set named "name"
     !Dummy arguments
     LOGICAL :: found
     CHARACTER (len=*) :: name ! set name
     TYPE (ele25_set), POINTER :: head, anter, posic

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
   END SUBROUTINE srch_ele25

   SUBROUTINE del_ele25 (head, anter, posic)

     !This subroutine deletes a set pointed with posic

     TYPE (ele25_set), POINTER :: head, anter, posic

     TYPE (ele25), POINTER :: ea,ep
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
       CALL del_ele25e (posic%head,posic%tail, ea, ep )  !deletes element
     END DO

     NULLIFY (posic,anter)          !point to nothing
   END SUBROUTINE del_ele25

   ! ******* functions for a list of elements in a set ********

   SUBROUTINE ini_ele25e (head, tail)
     !initialize a list of ELE13 elements

     TYPE (ele25), POINTER :: head, tail

     NULLIFY (head, tail)       !initializes first and last pointer

   END SUBROUTINE ini_ele25e

   SUBROUTINE new_ele25e(elm)
   !Create a new element of ELE25 sets

     TYPE(ele25),POINTER:: elm

     ALLOCATE(elm)
     elm%numel = 0        !Initialize label of element
     elm%matno = 0        !     "     material number
     elm%lnods = 0        !     "     conectivities
     elm%lside(1:4) = 0   !     "     neighbour elements
     !elm%ldv2sd(1:4) = 1  !     "     status of refined side
     elm%Area  = 0d0      !     "     initial area
     elm%angle = 0d0      !     "     angle between dir 1 and orthotropic dir 1
     elm%lb = 0d0         !     "     initial thickness ratio
     elm%cartd = 0d0      !     "     cartesian derivatives
     elm%normd = 0d0
     elm%ns = 0d0         !     "     side normals
     elm%a0 = 0d0         !     "     initial angles
     elm%ci = 0d0
     elm%gamma = 0d0
     elm%stra0 = 0d0      !     "     original curvature tensor
     elm%stra1 = 0d0      !     "     actual mid-surface metric tensor
     elm%sf = 0d0         !     "     stabilization forces
     ALLOCATE( elm%si(4))
     NULLIFY(elm%gausv,elm%si(1)%p,elm%si(2)%p,elm%si(3)%p,elm%si(4)%p)
     NULLIFY(elm%next)

   RETURN
   END SUBROUTINE new_ele25e

   SUBROUTINE add_ele25e (new, head, tail)
     !This subroutine adds data to the end of the list
     !Dummy arguments
     TYPE (ele25), POINTER :: new, head, tail

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
   END SUBROUTINE add_ele25e

   SUBROUTINE srch_ele25e (head, anter, posic, kelem, found)
     !This subroutine searches for an element labeled "kelem"
     !Dummy arguments
     LOGICAL :: found
     INTEGER (kind=4) :: kelem
     TYPE (ele25), POINTER :: head, anter, posic

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
   END SUBROUTINE srch_ele25e

   SUBROUTINE del_ele25e (head, tail, anter, posic)

     !This subroutine deletes element pointed with posic

     TYPE (ele25), POINTER :: head, tail, anter, posic

     IF (.NOT.ASSOCIATED (anter)) THEN    !
       head => posic%next
     ELSE
       anter%next => posic%next
     END IF
     IF( .NOT.ASSOCIATED(posic%next) )tail => anter !last element in list
     IF(ASSOCIATED(posic%gausv)) DEALLOCATE (posic%gausv) !deallocate variable arrays
     DEALLOCATE (posic)                    !deallocate fixed space
     IF( ASSOCIATED( anter) )posic => anter%next                   !point to next element
     RETURN
   END SUBROUTINE del_ele25e

   SUBROUTINE cut_ele25e (head, anter, posic)
     !This subroutine deletes a set pointed with posic
     ! without nullifying anter    ???? what for ????
     TYPE (ele25), POINTER :: head, anter, posic

     IF (.NOT.ASSOCIATED (anter)) THEN
       head => posic%next
     ELSE
       anter%next => posic%next
     ENDIF
     NULLIFY (posic)
   END SUBROUTINE cut_ele25e

   SUBROUTINE dalloc25e(elem)
     ! deallocates an element pointed with elem (releasing memory)
     TYPE (ele25), POINTER :: elem
     IF(ASSOCIATED(elem%gausv))DEALLOCATE (elem%gausv)              !deallocate variable arrays
     NULLIFY(elem%next)
     DEALLOCATE(elem)
   END SUBROUTINE dalloc25e

   INCLUDE 'acvd25.fi'
   INCLUDE 'axep25.fi'
   INCLUDE 'bbra25.fi'
   INCLUDE 'bfle25.fi'
   INCLUDE 'bmem25.fi'
   INCLUDE 'bran25.fi'
   INCLUDE 'comm25.fi'
   INCLUDE 'delt25.fi'
   INCLUDE 'dump25.fi'
   INCLUDE 'elmd25.fi'
   INCLUDE 'gaus25.fi'
   INCLUDE 'luma25.fi'
   INCLUDE 'mase25.fi'
   INCLUDE 'nghb25.fi'
   INCLUDE 'outd25.fi'
   INCLUDE 'rest25.fi'
   INCLUDE 'resv25.fi'
   INCLUDE 'slnods25.fi'
   INCLUDE 'stra25.fi'
   INCLUDE 'surf25.fi'
   INCLUDE 'updl25.fi'

 END MODULE ele25_db

