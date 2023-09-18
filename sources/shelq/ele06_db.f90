  MODULE ele06_db
   USE param_db,ONLY: mnam,midn,mlin
   USE mat_dba, ONLY: section,sect_search,psecs,pmats,mater,postv
   USE c_input
    IMPLICIT NONE
    INTEGER (kind=4), PARAMETER :: &
      nnode = 4, &  ! number of nodes per element
      ngaus = 4, &  ! number of GP
      ngamm = 4     ! number of shear strain values

    SAVE

   REAL (kind=8) :: posgp(2,ngaus),    & !Gauss point position
                    shape(nnode,ngaus)   !shape functions

   TYPE ele06
     INTEGER (kind=4) :: numel  ! label of element
     INTEGER (kind=4) :: matno  ! Material number
     INTEGER (kind=4) :: lnods(nnode)  ! Conectivities
     REAL (kind=8) :: angle                  ! angle of local system
     REAL (kind=8) :: dvolu(ngaus),        & ! Gauss associated area
                      cartd(nnode,2,ngaus),& ! Cartesian derivatives
                      jacin(2,2,ngaus),    & ! Jacobian inverse
                      stra0(6,ngaus),      & ! initial strains
                      gamm0(ngamm),        & ! initial shear strains
                      ambda(2,ngaus)         ! thickness ratio
     REAL (kind=8), POINTER :: & !output & plasticity variables
                      jas(:,:),        & !(2,4)       side vector at mid point
                      nangl(:,:),      & !(2,nnode)   Cos - Sin of nodal angles
                      strsg(:,:),      & !(nstre,ngaus)   present forces and moments
                      ehist(:,:),      & !(5,ngaus)       first set
                      strap(:,:),      & !(nstre,ngaus)   second set
                      stres(:,:,:)       !(5,nlayr,ngaus) third set
     !ehist = (/ epstr-actual,eqstr,epstr-lastconverged,consist-param(2) /)
     TYPE (ele06), POINTER :: next              !pointer to next element
   END TYPE ele06

    ! Derived type for a set of SHELQ elements
    TYPE ele06_set
      CHARACTER (len=mnam) :: sname   ! set name
      INTEGER (kind=4) :: nstre = 8, &  ! number of stresses/GP
                          nelem, & ! number of elements
                          nreqs, & ! number of GP for hist. output
                          narch    ! number of output unit
      LOGICAL :: gauss           ! .FALSE. -> Initial constants not
                                 !  defined or not updated
      LOGICAL :: zigzag = .FALSE.  ! .TRUE. use superimposed zig-zag function
      INTEGER :: plstr           ! compute Plastic Strain Flag
         ! -1 from Cauchy stress  0 - do not   1 from 2nd P-K
      INTEGER :: locax           ! local x definition option
      REAL (kind=8) ::  angdf    ! angle between X_1 and orthotropic dir 1
      TYPE (ele06), POINTER    :: head, tail !pointer to first and last elm.
      INTEGER (kind=4), POINTER :: ngrqs(:)  !(nreqs) Points for history output
      TYPE (ele06_set), POINTER :: next
    END TYPE ele06_set
    TYPE (ele06_set), POINTER :: head
    TYPE (ele06_set), POINTER :: tail

  CONTAINS
    SUBROUTINE ini_ele06 (head, tail)
      !initialize a list of SHELQ sets

      TYPE (ele06_set), POINTER :: head, tail

      NULLIFY (head, tail)

    END SUBROUTINE ini_ele06

    SUBROUTINE add_ele06 (new, head, tail)
      !This subroutine adds data to the end of the list
      !Dummy arguments
      TYPE (ele06_set), POINTER :: new, head, tail

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
    END SUBROUTINE add_ele06

    SUBROUTINE srch_ele06 (head, anter, posic, name, found)
      !This subroutine searches for a set named "name"
      !Dummy arguments
      LOGICAL :: found
      CHARACTER (len=*) :: name ! set name
      TYPE (ele06_set), POINTER :: head, anter, posic

      found = .FALSE.
      NULLIFY (posic,anter)
      !Check if a list is empty
      IF (ASSOCIATED (head)) THEN
        posic => head
        DO
          IF (TRIM(posic%sname) == TRIM(name)) THEN
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
    END SUBROUTINE srch_ele06

    SUBROUTINE del_ele06 (head, anter, posic)
      !This subroutine deletes a set pointed with posic
      TYPE (ele06_set), POINTER :: head, anter, posic

      TYPE (ele06), POINTER :: ea,ep
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
        CALL del_ele06e (posic%head,posic%tail, ea, ep )  !deletes element
      END DO

     NULLIFY (posic,anter)          !point to nothing
    END SUBROUTINE del_ele06

   SUBROUTINE new_ele06(elset)
   !Create a new element of ELE06 sets

     TYPE(ele06_set),POINTER:: elset

     ALLOCATE(elset)
     elset%sname = ''       !Initialize set name
     elset%nelem = 0        !     "     number of elements
     elset%nreqs = 0        !     "     number of GP for hist. output
     elset%narch = 0        !     "     number of output unit
     elset%gauss = .FALSE.  !     "     flag to compute Gauss constants
     elset%plstr = 0        !     "     compute Plastic Strain Flag
     elset%locax = 3        !     "     local x definition
     elset%angdf = 0d0      !     "     angle between X_1 and orthotropic dir 1
     NULLIFY(elset%head,elset%tail,elset%ngrqs)
     NULLIFY(elset%next)

   RETURN
   END SUBROUTINE new_ele06

   SUBROUTINE ini_ele06e (head, tail)
     !initialize a list of ELE06 elements

     TYPE (ele06), POINTER :: head, tail

     NULLIFY (head, tail)       !initializes first and last pointer

   END SUBROUTINE ini_ele06e

   SUBROUTINE add_ele06e (new, head, tail)
     !This subroutine adds data to the end of the list
     !Dummy arguments
     TYPE (ele06), POINTER :: new, head, tail

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
   END SUBROUTINE add_ele06e

   SUBROUTINE srch_ele06e (head, anter, posic, kelem, found)
     !This subroutine searches for an element labeled "kelem"
     !Dummy arguments
     LOGICAL :: found
     INTEGER (kind=4) :: kelem
     TYPE (ele06), POINTER :: head, anter, posic

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
   END SUBROUTINE srch_ele06e

   SUBROUTINE del_ele06e (head, tail, anter, posic)

     !This subroutine deletes element pointed with posic

     TYPE (ele06), POINTER :: head, tail, anter, posic
     TYPE (ele06), POINTER :: e

     IF (.NOT.ASSOCIATED (anter)) THEN    !
       head => posic%next
     ELSE
       anter%next => posic%next
     END IF
     e => posic%next                       !keep pointer to next element
     IF( .NOT.ASSOCIATED(e) )tail => anter !last element in list
     IF(ASSOCIATED(posic%nangl))DEALLOCATE (posic%nangl,posic%jas)
     IF(ASSOCIATED(posic%strsg))DEALLOCATE (posic%strsg)
     IF(ASSOCIATED(posic%ehist))DEALLOCATE (posic%ehist,posic%strap)
     IF(ASSOCIATED(posic%stres))DEALLOCATE (posic%stres)
     DEALLOCATE (posic)                    !deallocate fixed space
     posic => e                            !point to next element
     ! NULLIFY (posic,anter)
     RETURN
   END SUBROUTINE del_ele06e

   SUBROUTINE cut_ele06e (head, anter, posic)
     !This subroutine deletes a set pointed with posic
     ! without nullifying anter    ???? what for ????
     TYPE (ele06), POINTER :: head, anter, posic

     IF (.NOT.ASSOCIATED (anter)) THEN
       head => posic%next
     ELSE
       anter%next => posic%next
     ENDIF
     NULLIFY (posic)
   END SUBROUTINE cut_ele06e

   SUBROUTINE new_ele06e(elm)
   !Create a new element of ELE06 sets

     TYPE(ele06),POINTER:: elm

     ALLOCATE(elm)
     elm%numel = 0        !Initialize label of element
     elm%matno = 0        !     "     material number
     elm%lnods(1:4) = 0   !     "     conectivities
     elm%angle = 0d0      !     "     angle between dir 1 and orthotropic dir 1
     elm%dvolu = 0d0      !     "     initial area
     elm%cartd = 0d0      !     "     cartesyan derivatives
     elm%jacin = 0d0      !     "     inverse jacobians
     elm%stra0 = 0d0      !     "     initial strains
     elm%gamm0 = 0d0      !     "     initial distorsions
     elm%ambda = 0d0      !     "     thickness ratio
     NULLIFY(elm%ehist,elm%strap,elm%stres,elm%strsg,elm%nangl,elm%jas)
     NULLIFY(elm%next)

   RETURN
   END SUBROUTINE new_ele06e


   INCLUDE 'acvdf6.fi'
   INCLUDE 'asstr6.fi'
   INCLUDE 'bmatx6.fi'
   INCLUDE 'bshem6.fi'
   INCLUDE 'codel6.fi'
   INCLUDE 'commv6.fi'
   INCLUDE 'comp_ang06.fi'
   INCLUDE 'dump06.fi'
   INCLUDE 'elmda6.fi'
   INCLUDE 'gauss6.fi'
   INCLUDE 'intrf6.fi'
   !INCLUDE 'loadp6.fi'
   INCLUDE 'lumas6.fi'
   INCLUDE 'masel6.fi'
   INCLUDE 'outdy6.fi'
   INCLUDE 'psib06.fi'
   INCLUDE 'rest06.fi'
   INCLUDE 'resvp6.fi'
   INCLUDE 'setga6.fi'
   INCLUDE 'stran6.fi'
   INCLUDE 'surf06.fi'
   INCLUDE 'updlo6.fi'

  END MODULE ele06_db
