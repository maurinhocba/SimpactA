 MODULE ele14_db
   USE param_db,ONLY: mnam,midn,mlin
   USE mat_dba, ONLY: section,sect_search,psecs,pmats,mater,postv,pmater,snn
   USE c_input
   USE kinc_db, ONLY : nndpd
   IMPLICIT NONE

   ! Derived type for an ELE14 element
   ! Basic thin Shell Triangle (BST)

   ! Reference: F.Flores & E.Oñate
   !            A Basic thin shell triangle for large strain plasticity
   !            Int. J. Num. Methods in Engng. (2001)

   INTEGER (kind=4), PARAMETER :: hh(2,3) = RESHAPE( &
                     [ 3,2, 1,3, 2,1 ], [ 2,3 ] ), &    !nodes of each side
                                  kk(3,3) = RESHAPE( &
                     [ 4,3,2, 5,1,3, 6,2,1 ], [3,3] ) !side element connectivities



   TYPE ele14
     INTEGER (kind=4) :: numel  ! label of element
     INTEGER (kind=4) :: matno  ! Material number
     INTEGER (kind=4) :: lnods(6)  ! Conectivities
     INTEGER (kind=4) :: lside(3)  ! Neighbour elements
     REAL (kind=8) :: Area1,     & ! Initial area
                      lb,        & ! Thickness ratio
                      angle,     & ! angle between dir 1 and orthotropic dir 1
                      a(3,0:3),  & ! side projections in dir 1
                      b(3,0:3),  & ! side projections in dir 2
                      ci(3),     & ! coefficients for curvature computations
                      stra0(3),  & ! Original curvature tensor
                      stra1(6)     ! Actual mid-surface metric tensor and curvatures
     REAL (kind=8), POINTER :: gausv(:,:)       !layer values
     TYPE (ele14), POINTER :: next              !pointer to next element
   END TYPE ele14

   ! Derived type for a set of ELE14 elements
   TYPE ele14_set
     CHARACTER (len=mnam) :: sname ! set name
     INTEGER (kind=4) :: nelem  ! number of elements
     INTEGER (kind=4) :: nreqs  ! number of GP for hist. output
     INTEGER (kind=4) :: narch  ! number of output unit
     LOGICAL :: logst           ! use logarithmic strain
     LOGICAL :: lside           ! .FALSE. -> topological arrays not
                                !  defined or not updated
     LOGICAL :: gauss           ! .FALSE. -> Initial constants not
                                !  defined or not updated
     LOGICAL :: nonrg           ! .TRUE. -> Non-regular
     LOGICAL :: origl           ! .TRUE. -> Original labels stored in connectivities
     INTEGER :: plstr           ! compute Plastic Strain Flag
         ! -1 from Cauchy stress  0 - do not   1 from 2nd P-K
     INTEGER :: locax           ! local x definition option
     REAL (kind=8) ::  angdf     ! angle between X_1 and orthotropic dir 1
     REAL (kind=8), POINTER :: stint(:,:)   !moments, forces and shears
     TYPE (ele14), POINTER    :: head, tail !pointer to first and last elm.
     INTEGER (kind=4), POINTER :: ngrqs(:)  !gauss points for output
     TYPE (ele14_set), POINTER :: next      !pointer to next set
   END TYPE ele14_set
   TYPE (ele14_set), POINTER, SAVE :: head, tail !first and last elements sets

 CONTAINS

   SUBROUTINE new_ele14(elset)
   !Create a new element of ELE14 sets

     TYPE(ele14_set),POINTER:: elset

     ALLOCATE(elset)
     elset%sname = ''       !Initialize set name
     elset%nelem = 0        !     "     number of elements
     elset%nreqs = 0        !     "     number of GP for hist. output
     elset%narch = 0        !     "     number of output unit
     elset%logst = .FALSE.  !     "     use logarithmic strain
     elset%lside = .FALSE.  !     "     flag to compute LSIDE
     elset%gauss = .FALSE.  !     "     flag to compute Gauss constants
     elset%nonrg = .FALSE.  !     "     flag to consider non-regular meshes
     elset%origl = .FALSE.  !     "     flag to consider original labels
     elset%plstr = 0        !     "     compute Plastic Strain Flag
     elset%locax = 3        !     "     local x definition
     elset%angdf = 0d0      !     "     angle between X_1 and orthotropic dir 1
     NULLIFY(elset%head,elset%tail,elset%ngrqs,elset%stint)
     NULLIFY(elset%next)

   RETURN
   END SUBROUTINE new_ele14

   SUBROUTINE add_ele14 (new, head, tail)
     !This subroutine adds data to the end of the list
     !Dummy arguments
     TYPE (ele14_set), POINTER :: new, head, tail

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
   END SUBROUTINE add_ele14

   SUBROUTINE srch_ele14 (head, anter, posic, name, found)
     !This subroutine searches for a set named "name"
     !Dummy arguments
     LOGICAL :: found
     CHARACTER (len=*) :: name ! set name
     TYPE (ele14_set), POINTER :: head, anter, posic

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
   END SUBROUTINE srch_ele14

   SUBROUTINE del_ele14 (head, anter, posic)

     !This subroutine deletes a set pointed with posic

     TYPE (ele14_set), POINTER :: head, anter, posic

     TYPE (ele14), POINTER :: ea,ep
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
       CALL del_ele14e (posic%head,posic%tail, ea, ep )  !deletes element
     END DO

     NULLIFY (posic,anter)          !point to nothing
   END SUBROUTINE del_ele14

   ! ******* functions for a list of elements in a set ********

   SUBROUTINE new_ele14e(elm)
   !Create a new element of ELE14 sets

     TYPE(ele14),POINTER:: elm

     ALLOCATE(elm)
     elm%numel = 0        !Initialize label of element
     elm%matno = 0        !     "     material number
     elm%lnods(1:6) = 0   !     "     conectivities
     elm%lside(1:3) = 0   !     "     neighbour elements
     elm%area1 = 0d0      !     "     initial area
     elm%angle = 0d0      !     "     angle between dir 1 and orthotropic dir 1
     elm%lb = 0d0         !     "     initial area
     elm%a = 0d0          !     "     side projections in dir 1
     elm%b = 0d0          !     "       "       "       "  "  2
     elm%stra0(1:3) = 0d0 !     "     original curvature tensor
     elm%stra1(1:6) = (/ 1d0,1d0,0d0,0d0,0d0,0d0 /) !     "     actual mid-surface metric tensor
     elm%ci(3)      = 1d0 !reserve space for coefficients
     NULLIFY(elm%gausv)
     NULLIFY(elm%next)

   RETURN
   END SUBROUTINE new_ele14e

   SUBROUTINE add_ele14e (new, head, tail)
     !This subroutine adds data to the end of the list
     !Dummy arguments
     TYPE (ele14), POINTER :: new, head, tail

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
   END SUBROUTINE add_ele14e

   SUBROUTINE ins_ele14e(new, top, tail)
   !This subroutine is simmilar to add_ele, but insert a member list between two pointers

     !--- Dummy variables
     TYPE(ele14),POINTER:: new, top, tail

     IF (ASSOCIATED(top%next)) THEN
       new%next => top%next
       top%next => new
     ELSE
       NULLIFY(new%next)
       top%next => new
       !It should update the tail of the list if necessary
       tail => new
     ENDIF

   RETURN
   END SUBROUTINE ins_ele14e

   SUBROUTINE srch_ele14e (head, anter, posic, kelem, found)
     !This subroutine searches for an element labeled "kelem"
     !Dummy arguments
     LOGICAL :: found
     INTEGER (kind=4) :: kelem
     TYPE (ele14), POINTER :: head, anter, posic

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
   END SUBROUTINE srch_ele14e

   SUBROUTINE del_ele14e (head, tail, anter, posic)

     !This subroutine deletes element pointed with posic

     TYPE (ele14), POINTER :: head, tail, anter, posic

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
   END SUBROUTINE del_ele14e

   SUBROUTINE cut_ele14e(head, tail, anter, posic)
     !This subroutine deletes a set pointed with posic
     ! without nullifying anter    ???? what for ????
     TYPE(ele14),POINTER:: head, tail, anter, posic

     IF (.NOT.ASSOCIATED (anter)) THEN
       head => posic%next
     ELSE
       anter%next => posic%next
     ENDIF
     ! if posic == tail
     IF (.NOT.ASSOCIATED(posic%next)) tail=>anter
     CALL dalloc14e(posic)

   RETURN
   END SUBROUTINE cut_ele14e

   SUBROUTINE dalloc14e(elem)
     ! deallocates an element pointed with elem (releasing memory)
     TYPE (ele14), POINTER :: elem
     IF(ASSOCIATED(elem%gausv))DEALLOCATE (elem%gausv)              !deallocate variable arrays
     NULLIFY(elem%next)
     DEALLOCATE(elem)
   END SUBROUTINE dalloc14e

   SUBROUTINE copy_ele14e(elold,elnew)
   !---------------------------------------------------
   !Copy a member of the list ELE14
   !---------------------------------------------------
   IMPLICIT NONE

     !Dummy arguments
     TYPE(ele14),POINTER:: elold, elnew
     !Local variables
     INTEGER(kind=4)::  s1, s2 !,sz

     elnew%numel = elold%numel    ! label of element
     elnew%matno = elold%matno    ! Material number
     elnew%Area1 = elold%Area1    ! Initial area
     elnew%lb = elold%lb          ! Thickness ratio
     elnew%angle = elold%angle    ! Angle between dir 1 and orthotropic dir 1

     elnew%lnods  = elold%lnods   ! Conectivities
     elnew%lside  = elold%lside   ! Neighbour elements
     elnew%a      = elold%a       ! side projections in dir 1
     elnew%b      = elold%b       ! side projections in dir 2
     elnew%ci     = elold%ci      ! coefficients for curvature computation
     elnew%stra0  = elold%stra0   ! Original curvature tensor
     elnew%stra1  = elold%stra1   ! Actual mid-surface metric tensor

     IF (ASSOCIATED(elold%gausv)) THEN
       s1 = SIZE(elold%gausv,DIM=1) ; s2 = SIZE(elold%gausv,DIM=2)
       ALLOCATE(elnew%gausv(s1,s2))
       elnew%gausv(1:s1,1:s2) = elold%gausv(1:s1,1:s2)   !layer values
     END IF

     elnew%next => elold%next

   RETURN
   END SUBROUTINE copy_ele14e

   INCLUDE 'acvd14.fi'
   INCLUDE 'axep14.fi'
   INCLUDE 'bfle14.fi'
   !INCLUDE 'boun14.fi' !not used presently
   INCLUDE 'comm14.fi'
   INCLUDE 'delt14.fi'
   INCLUDE 'dump14.fi'
   INCLUDE 'elmd14.fi'
   INCLUDE 'expo14.fi'
   INCLUDE 'gaus14.fi'
   INCLUDE 'impo14.fi'
   INCLUDE 'inig14.fi'
   INCLUDE 'luma14.fi'
   INCLUDE 'mase14.fi'
   INCLUDE 'nods14.fi'
   INCLUDE 'outd14.fi'
   INCLUDE 'rest14.fi'
   INCLUDE 'resv14.fi'
   INCLUDE 'secd14.fi'
   INCLUDE 'slum14.fi'
   INCLUDE 'stra14.fi'
   INCLUDE 'stst14.fi'
   INCLUDE 'surf14.fi'
   INCLUDE 'toar14.fi'
   INCLUDE 'updl14.fi'

 END MODULE ele14_db
