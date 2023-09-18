 MODULE ele15_db
   USE param_db,ONLY: mnam,midn,mlin
   USE mat_dba, ONLY: section,sect_search,mater,postv,psecs
   USE c_input
   IMPLICIT NONE

   ! Derived type for an ELE15 element
   ! Basic thin Shell Triangle (BST) + reinforcement

   ! Reference: F.Flores Umpublished


   SAVE

   TYPE reinf
     INTEGER (kind=4) :: numel              !element label
     INTEGER (kind=4) :: lnods(4)           !1:2  = nodes  3:4 elements
     INTEGER (kind=4) :: secno              !section number
     REAL (kind=8)    :: l0                 !nn   = initial length
     REAL (kind=8)    :: curv(2,2)          !actual curvatures
     REAL (kind=8)    :: fc(2)              !factors
     REAL (kind=8)    :: dc(2)              !direction cosines
     REAL (kind=8)    :: stint(3)           !force and moments
     TYPE (reinf), POINTER :: next
   END TYPE reinf

   !Derived type for an array of pointers to reinforcements
   TYPE preinf
      TYPE (reinf), POINTER :: p
   END TYPE preinf

   !   hh= side element connectivities
   INTEGER(kind=4), PARAMETER :: hh(3,3) = RESHAPE((/ 4,3,2, 5,1,3, 6,2,1 /), (/3,3/) )
   REAL(kind=8), PARAMETER :: alp1 = 0.70d0, & !Max angle to apply membrane cuadratic approach
                              alp2 = 1.00d0, & !angle to use standard CST for membrane
                              alp3 = 0.3d0     !alp2-alp1 auxiliar value

   TYPE ele15
     INTEGER (kind=4) :: numel  ! label of element
     INTEGER (kind=4) :: matno  ! Material number
     INTEGER (kind=4) :: lnods(6)  ! Conectivities
     INTEGER (kind=4) :: lside(3)  ! Neighbour elements
     REAL (kind=8) :: Area1,     & ! Initial area
                      lb,        & ! Thickness ratio
                      angle,     & ! angle between dir 1 and orthotropic dir 1
                      a(3),      & ! side projections in dir 1
                      b(3),      & ! side projections in dir 2
                      c(3,3,2),  & ! side projections in dir n
                      cd(4,2,3), & ! cartesian derivatives at mid-side points
                      a0(3),     & ! initial angles with side elements
                      ci(3),     & ! coefficients for angle change
                      gamma(3),  & ! distorsion at each side
                      stra1(6)     ! Actual mid-surface metric tensors

     REAL (kind=8), POINTER :: gausv(:,:)       !layer values
     LOGICAL :: mems(3)                         !use quadratic approach for that side
     TYPE (preinf), POINTER :: si(:)
     TYPE (ele15), POINTER :: next              !pointer to next element
   END TYPE ele15

   ! Derived type for a set of ELE15 elements
   TYPE ele15_set
     CHARACTER (len=mnam) :: sname ! set name
     INTEGER (kind=4) :: nelem  ! number of elements
     INTEGER (kind=4) :: nreqs  ! number of GP for hist. output
     INTEGER (kind=4) :: narch  ! number of output unit
     INTEGER (kind=4) :: nrf    ! number of reinforced sides
     LOGICAL :: logst           ! use logarithmic strain
     LOGICAL :: lside           ! .FALSE. -> topological arrays not
                                !  defined or not updated
     LOGICAL :: gauss           ! .FALSE. -> Initial constants not
                                !  defined or not updated
     INTEGER :: plstr           ! compute Plastic Strain Flag
         ! -1 from Cauchy stress  0 - do not   1 from 2nd P-K
     REAL (kind=8) ::  angdf     ! angle between X_1 and orthotropic dir 1
     REAL (kind=8), POINTER :: stint(:,:)   !forces, moments, shears & Mises stresses
     TYPE (ele15), POINTER    :: head, tail !pointer to first and last elm.
     TYPE (reinf), POINTER :: rhead , rtail !pointers to branching data base
     INTEGER (kind=4), POINTER :: ngrqs(:)  !gauss points for output
     ! for shear evaluation
     INTEGER (kind=4) :: shear              !compute shear forces
     INTEGER (kind=4), POINTER :: ninv(:)   !inverse nodal relation (temporary, no dumping)
     REAL (kind=8), POINTER :: shears(:,:) !smoothed nodal moments (temporary, no dumping)
     REAL (kind=8), POINTER :: factors(:)   !nodal factors for smoothing (temporary, no dumping)
     ! end
     TYPE (ele15_set), POINTER :: next      !pointer to next set
   END TYPE ele15_set
   TYPE (ele15_set), POINTER, SAVE :: head,tail !first and last elements sets

 CONTAINS

    FUNCTION atan4(a,b,c,d)
    IMPLICIT NONE
    REAL(kind=8) :: atan4
    REAL(kind=8), INTENT(IN) :: a,b,c,d
    REAL (kind=8),PARAMETER :: twopi=6.283185307179586d0, &
                                  pi=3.141592653589793d0
      atan4 = ATAN2(a,b) - c
      !  limit angle change to Pi (180 degrees)
      IF( atan4-d > pi  )THEN
        atan4 = atan4 - twopi
      ELSE IF( atan4-d < -pi )THEN
        atan4 = atan4 + twopi
      END IF
    END FUNCTION

   SUBROUTINE new_ele15(elset)
   !Create a new element of ELE15 sets

     TYPE(ele15_set),POINTER:: elset

     ALLOCATE(elset)
     elset%sname = ''       !Initialize set name
     elset%nelem = 0        !     "     number of elements
     elset%nreqs = 0        !     "     number of GP for hist. output
     elset%narch = 0        !     "     number of output unit
     elset%nrf   = 0        !     "     number of reinforcement sides
     elset%logst = .FALSE.  !     "     use logarithmic strain
     elset%lside = .FALSE.  !     "     flag to compute LSIDE
     elset%gauss = .FALSE.  !     "     flag to compute Gauss constants
     elset%shear = 0        !     "     flag to compute shear forces
     elset%plstr = 0        !     "     compute Plastic Strain Flag
     elset%angdf = 0d0      !     "     angle between X_1 and orthotropic dir 1
     NULLIFY(elset%head,elset%tail,elset%ngrqs,elset%stint,elset%rhead,elset%rtail)
     NULLIFY(elset%shears,elset%factors,elset%ninv)
     NULLIFY(elset%next)

   RETURN
   END SUBROUTINE new_ele15

   SUBROUTINE ini_ele15 (head, tail)
     !initialize a list of ELE15 sets

     TYPE (ele15_set), POINTER :: head, tail

     NULLIFY (head, tail)

   END SUBROUTINE ini_ele15

   SUBROUTINE add_ele15 (new, head, tail)
     !This subroutine adds data to the end of the list
     !Dummy arguments
     TYPE (ele15_set), POINTER :: new, head, tail

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
   END SUBROUTINE add_ele15

   SUBROUTINE srch_ele15 (head, anter, posic, name, found)
     !This subroutine searches for a set named "name"
     !Dummy arguments
     LOGICAL :: found
     CHARACTER (len=*) :: name ! set name
     TYPE (ele15_set), POINTER :: head, anter, posic

     found = .FALSE.                     !initializes flag
     NULLIFY (posic,anter)               !initializes pointers
     !Check if a list is empty
     IF (ASSOCIATED (head)) THEN         !if there are sets
       posic => head                     !point to first set
       DO
         IF(posic%sname == name) THEN    !check name
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
   END SUBROUTINE srch_ele15

   SUBROUTINE del_ele15 (head, anter, posic)

     !This subroutine deletes a set pointed with posic

     TYPE (ele15_set), POINTER :: head, anter, posic

     TYPE (ele15), POINTER :: ea,ep
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
       CALL del_ele15e (posic%head,posic%tail, ea, ep )  !deletes element
     END DO

     NULLIFY (posic,anter)          !point to nothing
   END SUBROUTINE del_ele15

   ! ******* functions for a list of elements in a set ********

   SUBROUTINE ini_ele15e (head, tail)
     !initialize a list of ELE15 elements

     TYPE (ele15), POINTER :: head, tail

     NULLIFY (head, tail)       !initializes first and last pointer

   END SUBROUTINE ini_ele15e

   SUBROUTINE new_ele15e(elm)
   !Create a new element of ELE15 sets

     TYPE(ele15),POINTER:: elm
     INTEGER(kind=4) :: i

     ALLOCATE(elm)
     elm%numel = 0        !Initialize label of element
     elm%matno = 0        !     "     material number
     elm%lnods(1:6) = 0   !     "     conectivities
     elm%lside(1:3) = 0   !     "     neighbour elements
     !elm%ldv2sd(1:3) = 1  !     "     status of refined side
     elm%Area1 = 0d0      !     "     initial area
     elm%angle = 0d0      !     "     angle between dir 1 and orthotropic dir 1
     elm%lb = 0d0         !     "     initial area
     elm%a = 0d0          !     "     side projections in dir 1
     elm%b = 0d0          !     "     side projections in dir 2
     elm%c = 0d0          !     "     side projections in dir n
     elm%cd = 0d0         !     "     cartesian derivatives
     elm%a0 = 0d0         !     "     initial angles
     elm%gamma = 0d0      !     "     distorsions
     elm%stra1(1:3) = 0d0 !     "     actual mid-surface metric tensor
     elm%mems = .TRUE.    !     "
     ALLOCATE(elm%si(3))
     DO i=1,3
       NULLIFY(elm%si(i)%p)
     END DO
     NULLIFY(elm%gausv)
     NULLIFY(elm%next)

   RETURN
   END SUBROUTINE new_ele15e

   SUBROUTINE add_ele15e (new, head, tail)
     !This subroutine adds data to the end of the list
     !Dummy arguments
     TYPE (ele15), POINTER :: new, head, tail

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
   END SUBROUTINE add_ele15e

   SUBROUTINE srch_ele15e (head, anter, posic, kelem, found)
     !This subroutine searches for an element labeled "kelem"
     !Dummy arguments
     LOGICAL :: found
     INTEGER (kind=4) :: kelem
     TYPE (ele15), POINTER :: head, anter, posic

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
   END SUBROUTINE srch_ele15e

   SUBROUTINE del_ele15e (head, tail, anter, posic)

     !This subroutine deletes element pointed with posic

     TYPE (ele15), POINTER :: head, tail, anter, posic
     TYPE (ele15), POINTER :: e

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
   END SUBROUTINE del_ele15e

   SUBROUTINE cut_ele15e (head, anter, posic)
     !This subroutine deletes a set pointed with posic
     ! without nullifying anter    ???? what for ????
     TYPE (ele15), POINTER :: head, anter, posic

     IF (.NOT.ASSOCIATED (anter)) THEN
       head => posic%next
     ELSE
       anter%next => posic%next
     ENDIF
     NULLIFY (posic)
   END SUBROUTINE cut_ele15e

   SUBROUTINE new_reinf(side)
   IMPLICIT NONE
   TYPE (reinf), POINTER :: side

     ALLOCATE(side)

      side%numel  = 0         !element label
      side%lnods  = 0         !1:2  = nodes  3:4 elements
      side%secno  = 0         !section number
      side%l0     = 0d0       !nn   = initial length
      side%curv   = 0d0       !initial curvatures
      side%fc     = 0d0       !factors
      side%stint  = 0d0       !force and moments
      NULLIFY(side%next)      !pointer

     RETURN
   END SUBROUTINE new_reinf

   INCLUDE 'acvd15.fi'
   INCLUDE 'axep15.fi'
   INCLUDE 'bfle15.fi'
   INCLUDE 'bmem15.fi'
   INCLUDE 'comm15.fi'
   INCLUDE 'delt15.fi'
   INCLUDE 'dump15.fi'
   INCLUDE 'elmd15.fi'
   INCLUDE 'elmd15r.fi'
   INCLUDE 'gaus15.fi'
   !INCLUDE 'load15.fi'
   INCLUDE 'mase15.fi'
   INCLUDE 'luma15.fi'
   INCLUDE 'outd15.fi'
   INCLUDE 'rest15.fi'
   INCLUDE 'resv15.fi'
   INCLUDE 'stra15.fi'
   !INCLUDE 'stst15.fi'
   INCLUDE 'toar15.fi'

 END MODULE ele15_db
