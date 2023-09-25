MODULE ele19_db
  USE param_db,ONLY: mnam,midn
  USE c_input
  USE mat_dba, ONLY: section,sect_search,psecs,pmats,mater,curve,inte_cr,cur_point  !USE POINTERS
  USE meshmo_db,ONLY: strnd, numpo, nodset, r_last, itransf, sg_old, &
                      ne_new, maxnn, l_new, sg_new, r_elm_zone,      &
                      r_crit, r_l_dist, r_b_dist, r_z_actv, nodlb
  USE ctrl_db, ONLY : therm, ndoft, itemp, tref, tdtime
  USE npo_db, ONLY : iftmp, tempe, dtemp
  IMPLICIT NONE

  ! Derived type for an ELE19 element
  ! Total & Updated Lagrangian 2-D triangular element for coupled thermo-mechanical analysis

  !  Reference: ??
  !
  !

  INTEGER, PARAMETER :: nnode = 3, & !number of nodes per element
                        ngaus = 1    !number of integration points

  INTEGER(kind=4), PARAMETER  :: kk(2,3) = RESHAPE((/ 3,2, 1,3, 2,1 /), &
                                                    (/2,3/) )

  TYPE ele19
    INTEGER (kind=4) :: numel  ! label of element
    INTEGER (kind=4) :: matno  ! Material number
    INTEGER (kind=4) :: lnods(6)  ! Conectivities
    REAL (kind=8) :: Area1,       & ! Initial area
                     angle,       & ! angle between dir 1 and orthotropic dir 1
                     cd(6,3),     & ! cartesian derivatives
                     stint(4),    & ! Actual stresses
                     beta(3),     & ! elemental distortion value for remeshing
                     pwork(3)       ! coupled thermo-mech factors (plastic power, elastic strain change,
                                    ! and plastic dissipation for temperature change)
    LOGICAL       :: repla   ! .TRUE. -> element must be replaced (distorted)
    REAL (kind=8), POINTER :: gausv(:), &    !Gauss-point internal variables
                              strpl(:)       !plastic strain for thermal coupled problem
    TYPE (ele19), POINTER :: next            !pointer to next element
  END TYPE ele19

  ! Derived type for a set of ELE19 elements
  TYPE ele19_set
    CHARACTER (len=mnam) :: sname ! set name
    INTEGER (kind=4) :: nelem, &  ! number of elements
                        nreqs, &  ! number of GP for hist. output
                        narch     ! number of output unit

    LOGICAL :: gauss, &        ! .FALSE. -> Initial constants not defined or not updated
               lside, &        ! .FALSE. -> patch node not computed yet
               origl, &        ! .TRUE. if original labels appears in connectivities
               eulrf, &        ! .TRUE. for Eulerian Formulation (Garino J2 algorithm)
               swapc           ! .FALSE. not swap the corner alone elements (default is TRUE)

    INTEGER :: plstr           ! compute Plastic Strain Flag
        ! -1 from Cauchy stress  0 - do not   1 from 2nd P-K
    REAL (kind=8) :: angdf     ! angle between X_1 and orthotropic dir 1

    TYPE (ele19), POINTER    :: head, tail !pointer to first and last elm.
    INTEGER (kind=4), POINTER :: ngrqs(:)  !gauss points for output
    INTEGER (kind=4), POINTER :: eside(:,:)!neigbor elements in patch
    TYPE (ele19_set), POINTER :: next      !pointer to next set
  END TYPE ele19_set

  TYPE(ele19_set),POINTER,SAVE:: head=>NULL(), tail=>NULL()  !first and last elements sets

 CONTAINS

  SUBROUTINE new_ele19(elms)
  !Create a new ELE19 sets of the list

    !--- Dummy variables
    TYPE(ele19_set),POINTER:: elms

    ALLOCATE(elms)   !Create ele19_set element
    elms%sname = ''        !Initialise set name
    elms%nelem = 0         !     "     number of elements
    elms%nreqs = 0         !     "     number of GP for hist. output
    elms%narch = 0         !     "     number of output unit
    elms%gauss = .FALSE.   !     "     flag for initial constants defined or updated
    elms%lside = .FALSE.   !     "     flag to compute Gauss constants
    elms%eulrf = .TRUE.    !     "     flag to consider Finger tensor
    elms%origl = .FALSE.   !     "     flag indicates original labels
    elms%swapc = .TRUE.    !     "     flag to swap triangles with one neighbor
    elms%plstr = 0         !     "     compute Plastic Strain Flag
    elms%angdf = 0d0       !     "     angle between X_1 and orthotropic dir 1
    NULLIFY(elms%head, elms%tail, elms%ngrqs, elms%eside)  !Initialises pointer
    NULLIFY(elms%next)   !Initialises pointer to next element of the list

  RETURN
  END SUBROUTINE new_ele19

  SUBROUTINE add_ele19(new, head, tail)
    !This subroutine adds data to the end of the list
    !Dummy arguments
    TYPE (ele19_set), POINTER :: new, head, tail

    !Check if a list is empty
    IF (.NOT. ASSOCIATED (head)) THEN
      !list is empty, start it
      head => new
      tail => new
      NULLIFY (tail%next)

    ELSE
      !add an element to the list
      tail%next => new
      NULLIFY (new%next)
      tail => new

    ENDIF
  END SUBROUTINE add_ele19

  SUBROUTINE srch_ele19(head, anter, posic, name, found)
    !This subroutine searches for a set named "name"
    !Dummy arguments
    LOGICAL :: found
    CHARACTER (len=*) :: name ! set name
    TYPE (ele19_set), POINTER :: head, anter, posic

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
  END SUBROUTINE srch_ele19

  SUBROUTINE del_ele19(head, anter, posic)
    !This subroutine deletes a set pointed with posic

    TYPE (ele19_set), POINTER :: head, anter, posic

    TYPE (ele19), POINTER :: ea,ep
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
      CALL del_ele19e (posic%head,posic%tail, ea, ep )  !deletes element
    END DO

    NULLIFY (posic,anter)          !point to nothing
  END SUBROUTINE del_ele19

  ! ******* functions for a list of elements in a set ********

  SUBROUTINE new_ele19e(elm)
  !Create a new ELE19 sets of the list

    !--- Dummy variables
    TYPE(ele19),POINTER:: elm

    ALLOCATE(elm)   !Create ele19 element
    elm%numel = 0     !Initialise label of element
    elm%matno = 0     !     "     material number
    elm%lnods = 0     !     "     conectivities
    elm%Area1 = 0d0   !     "     initial area
    elm%angle = 0d0   !     "     angle between dir 1 and orthotropic dir 1
    elm%cd = 0d0      !     "     cartesian derivatives
    elm%stint = 0d0   !     "     actual stresses
    elm%beta = 1d0    !     "     element distortion
    elm%pwork = 0d0   !     "     coupled thermo-mechanical dissipation factors
    elm%repla = .FALSE. !Initializes the flag for replace distorted element
    NULLIFY(elm%gausv)  !Initialises pointer
    NULLIFY(elm%next)   !Initialises pointer to next element of the list

  RETURN
  END SUBROUTINE new_ele19e

  SUBROUTINE add_ele19e(new, head, tail)
    !This subroutine adds data to the end of the list
    !Dummy arguments
    TYPE (ele19), POINTER :: new, head, tail

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
  END SUBROUTINE add_ele19e

  SUBROUTINE srch_ele19e(head, anter, posic, kelem, found)
    !This subroutine searches for an element labeled "kelem"
    !Dummy arguments
    LOGICAL :: found
    INTEGER (kind=4) :: kelem
    TYPE (ele19), POINTER :: head, anter, posic

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
  END SUBROUTINE srch_ele19e

  SUBROUTINE del_ele19e(head, tail, anter, posic)
    !This subroutine deletes element pointed with posic

    TYPE (ele19), POINTER :: head, tail, anter, posic
    TYPE (ele19), POINTER :: e

    IF (.NOT.ASSOCIATED (anter)) THEN    !
      head => posic%next
    ELSE
      anter%next => posic%next
    END IF
    e => posic%next                       !keep pointer to next element
    IF( .NOT.ASSOCIATED(e) )tail => anter !last element in list
    IF( ASSOCIATED( posic%gausv ) )DEALLOCATE (posic%gausv)          !deallocate internal variables
    !IF( ASSOCIATED( posic%naxis ) )DEALLOCATE( posic%naxis )
    DEALLOCATE (posic)                    !deallocate fixed space
    posic => e                            !point to next element
    ! NULLIFY (posic,anter)
    RETURN
  END SUBROUTINE del_ele19e

  SUBROUTINE cut_ele19e(head, anter, posic)
    !This subroutine deletes a set pointed with posic
    ! without nullifying anter    ???? what for ????
    TYPE (ele19), POINTER :: head, anter, posic

    IF (.NOT.ASSOCIATED (anter)) THEN
      head => posic%next
    ELSE
      anter%next => posic%next
    ENDIF
    NULLIFY (posic)
  END SUBROUTINE cut_ele19e

  SUBROUTINE delete_ele19e(head, tail)
  !This subroutine deletes element pointed with posic
  IMPLICIT NONE
    !--- Dummy variables
    TYPE(ele19),POINTER:: head, tail
    !--- Local variables
    TYPE(ele19),POINTER:: elm

    NULLIFY(tail)
    DO
      IF (.NOT.ASSOCIATED(head)) EXIT
      elm => head%next                  !keep pointer to next element
      DEALLOCATE(head%gausv)            !deallocate internal variables
      DEALLOCATE(head)                  !deallocate element
      head => elm                       !point to next element
    END DO
  RETURN
  END SUBROUTINE delete_ele19e

  INCLUDE 'acvd19.fi'
  INCLUDE 'axep19.fi'
  INCLUDE 'bmat19.fi'
  INCLUDE 'boun19.fi'
  INCLUDE 'cdac19.fi'
  INCLUDE 'comm19.fi'
  INCLUDE 'delt19.fi'
  INCLUDE 'dist19.fi'
  INCLUDE 'dump19.fi'
  INCLUDE 'elmd19.fi'
  INCLUDE 'gaus19.fi'
  INCLUDE 'getv19.fi'
  INCLUDE 'luma19.fi'
  INCLUDE 'mase19.fi'
  INCLUDE 'nods19.fi'
  INCLUDE 'outd19.fi'
  INCLUDE 'real19.fi'
  INCLUDE 'rest19.fi'
  INCLUDE 'resv19.fi'
  INCLUDE 'slnd19.fi'
  INCLUDE 'smog19.fi'
  INCLUDE 'smom19.fi'
  INCLUDE 'smth19.fi'
  INCLUDE 'stra19.fi'
  INCLUDE 'surf19.fi'
  INCLUDE 'tdel19.fi'
  INCLUDE 'tlma19.fi'
  INCLUDE 'toar19.fi'
  INCLUDE 'tres19.fi'
  INCLUDE 'updl19.fi'

END MODULE ele19_db
