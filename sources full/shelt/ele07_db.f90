  MODULE ele07_db
     USE param_db,ONLY: mnam,midn,mlin
     USE mat_dba, ONLY: section,sect_search,psecs,pmats,mater,postv,snn
     USE c_input
    IMPLICIT NONE
    INTEGER (kind=4), PARAMETER :: &
      nnass = 6, &  ! number of assumed shear strain points
      nbeta = 9     !
              !                             r3   = SQRT(3)
              !    r3p1 =(r3+1d0)/2d0       r3m1 = (r3-1d0)/2d0
              !    r3p3 =(r3+3d0)/2d0       r3m3 = (r3-3d0)/2d0
              !    fxa  = 1/6+1/SQRT(12)    fxb  = 1/6-1/SQRT(12)    fxc = 2/3
              !    fxd  = 1+2/SQRT(3)       fxe  =-1+2/SQRT(3)       fxf = 4/SQRT(3)
     REAL (kind=8), PARAMETER ::            r3   = 1.73205080756890d0, &
                 r3p1 = 1.36602540378445d0, r3m1 = 0.36602540378445d0, &
                 r3p3 = 2.36602540378445d0, r3m3 =-0.63397459621555d0, &
                 fxa=  0.455341801261480d0, fxb= -0.122008467928146d0, &
                 fxc=  0.666666666666667d0, fxd=  2.154700538379250d0, &
                 fxe=  0.154700538379252d0, fxf=  2.309401076758500d0
     INTEGER( kind=4), PARAMETER :: kk(3,6) =(/ 1,4,2, 2,4,1, 2,5,3, 3,5,2, 3,6,1, 1,6,3 /)
     REAL( kind=8), PARAMETER :: nf(3) = (/ fxa, fxc, fxb /), nd(3) = (/ -fxd, fxf,-fxe /)
     !   hh= side element connectivities
     INTEGER(kind=4), PARAMETER :: hh(3,3) = RESHAPE((/ 4,3,2, 5,1,3, 6,2,1 /), (/3,3/) )

     SAVE

     REAL (kind=8) :: dn(6,2,4),        & !nodal shape functions derivatives
                      nfdas(6,2,nnass)    ! nodal shape functions derivatives

   TYPE ele07
     INTEGER (kind=4) :: numel  ! label of element
     INTEGER (kind=4) :: matno  ! Material number
     INTEGER (kind=4), POINTER  :: lnods(:)  ! Conectivities
     REAL (kind=8) :: angle                  ! angle of local system
     REAL (kind=8), POINTER :: dvolu(:),      & !(ngaus)         Gauss associated area
                               cartd(:,:,:),  & !(nnode,2,ngaus) Cartesian derivatives
                               cd(:,:,:),     & !(4,2,3)         Cartesian derivatives for QUAD approach
                               jacin(:,:,:),  & !(2,2,ngaus)     Jacobian inverse
                               stra0(:,:),    & !(6,ngaus)       initial strains
                               gamm0(:),      & !(ngamm)         initial shear strains
                               strsg(:,:),    & !(nstre,ngaus)   present forces and moments
                               qbar(:),       & !(ngamm)         present equivalent shears
                               ambda(:,:)       !(2,ngaus)       thickness ratio
     REAL (kind=8), POINTER :: beta(:)
     REAL (kind=8), POINTER :: & !plastictiy variables
                      nangl(:,:),      & !(2,nnode)   Cos - Sin of nodal angles
                      jas(:,:),        & !(2,ngamm)   local direction at assumed shear strain points
                      ehist(:,:),      & !(5,ngaus)       first set
                      strap(:,:),      & !(nstre,ngaus)   second set
                      stres(:,:,:)       !(5,nlayr,ngaus) third set
     !ehist = (/ epstr-actual,eqstr,epstr-lastconverged,consist-param(2) /)
     TYPE (ele07), POINTER :: next              !pointer to next element
   END TYPE ele07

    ! Derived type for a set of SHELT/TLLL elements
    TYPE ele07_set
      CHARACTER (len=mnam) :: sname   ! set name
      INTEGER (kind=4) :: nelem, & ! number of elements
                          nstre, & ! number of stress components (8 or 14)
                          nnode, & ! 6 9 number of nodes per element
                          stype, & ! formulation type 0:3
                          nreqs, & ! number of GP for hist. output
                          ngaus, & ! = 3  1 number of GP
                          ngamm, & ! = 6  3 number of shear strain values
                          narch    ! number of output unit
      REAL (kind=8),POINTER :: posgp(:,:),    & !(2,ngaus)        Gauss point position
                               weigp(:),      & !(ngaus)          Gauss point weigths
                               shape(:,:),    & !(nnode,ngaus)    nodal shape functions at Gauss points
                               ap1(:,:,:)       !(2,ngamm,ngaus)  shear interpolation matrix

      LOGICAL :: lside           ! .FALSE. -> topological arrays not
                                 !  defined or not updated
      LOGICAL :: gauss           ! .FALSE. -> Initial constants not
                                 !  defined or not updated
      LOGICAL :: zigzag = .FALSE.  ! .TRUE. use superimposed zig-zag function
      INTEGER :: plstr           ! compute Plastic Strain Flag
         ! -1 from Cauchy stress  0 - do not   1 from 2nd P-K
      INTEGER :: locax           ! local x definition option
      REAL (kind=8) ::  angdf    ! angle between X_1 and orthotropic dir 1
      REAL (kind=8) ::  stabq    ! stabilization factor for shear
      TYPE (ele07), POINTER    :: head, tail !pointer to first and last elm.
      INTEGER (kind=4), POINTER :: ngrqs(:)  !(nreqs) Points for history output
      TYPE (ele07_set), POINTER :: next
    END TYPE ele07_set
    TYPE (ele07_set), POINTER :: head
    TYPE (ele07_set), POINTER :: tail

  CONTAINS
    SUBROUTINE ini_ele07 (head, tail)
      !initialize a list of SHELT sets

      TYPE (ele07_set), POINTER :: head, tail

      NULLIFY (head, tail)

    END SUBROUTINE ini_ele07

    SUBROUTINE add_ele07 (new, head, tail)
      !This subroutine adds data to the end of the list
      !Dummy arguments
      TYPE (ele07_set), POINTER :: new, head, tail

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
    END SUBROUTINE add_ele07

    SUBROUTINE srch_ele07 (head, anter, posic, name, found)
      !This subroutine searches for a set named "name"
      !Dummy arguments
      LOGICAL :: found
      CHARACTER (len=*) :: name ! set name
      TYPE (ele07_set), POINTER :: head, anter, posic

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
    END SUBROUTINE srch_ele07

    SUBROUTINE del_ele07 (head, anter, posic)
      !This subroutine deletes a set pointed with posic
      TYPE (ele07_set), POINTER :: head, anter, posic

      TYPE (ele07), POINTER :: ea,ep
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
        CALL del_ele07e (posic%head,posic%tail, ea, ep )  !deletes element
      END DO

     NULLIFY (posic,anter)          !point to nothing
    END SUBROUTINE del_ele07

   SUBROUTINE new_ele07(elset)
   !Create a new element of ELE06 sets

     TYPE(ele07_set),POINTER:: elset

     ALLOCATE(elset)
     elset%sname = ''       !Initialize set name
     elset%nelem = 0        !     "     number of elements
     elset%nstre = 8        !     "     number of stresses
     elset%nreqs = 0        !     "     number of GP for hist. output
     elset%nnode = 6        !     "     number of element node
     elset%ngaus = 3        !     "     number of GP
     elset%ngamm = 6        !     "     number of shear strain values
     elset%narch = 0        !     "     number of output unit
     elset%gauss = .FALSE.  !     "     flag to compute Gauss constants
     elset%lside = .FALSE.  !     "     flag to compute side connectivities
     elset%zigzag= .FALSE.  !     "     flag to use additional in-plane displ.
     elset%plstr = 0        !     "     compute Plastic Strain Flag
     elset%locax = 3        !     "     local x definition
     elset%angdf = 0d0      !     "     angle between X_1 and orthotropic dir 1
     elset%stabq = 1d0      !     "     stabilization factor
     NULLIFY(elset%posgp,elset%weigp,elset%shape,elset%ap1)
     NULLIFY(elset%head,elset%tail,elset%ngrqs)
     NULLIFY(elset%next)

   RETURN
   END SUBROUTINE new_ele07

   SUBROUTINE ini_ele07e (head, tail)
     !initialize a list of ELE07 elements

     TYPE (ele07), POINTER :: head, tail

     NULLIFY (head, tail)       !initializes first and last pointer

   END SUBROUTINE ini_ele07e

   SUBROUTINE add_ele07e (new, head, tail)
     !This subroutine adds data to the end of the list
     !Dummy arguments
     TYPE (ele07), POINTER :: new, head, tail

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
   END SUBROUTINE add_ele07e

   SUBROUTINE srch_ele07e (head, anter, posic, kelem, found)
     !This subroutine searches for an element labeled "kelem"
     !Dummy arguments
     LOGICAL :: found
     INTEGER (kind=4) :: kelem
     TYPE (ele07), POINTER :: head, anter, posic

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
   END SUBROUTINE srch_ele07e

   SUBROUTINE del_ele07e (head, tail, anter, posic)

     !This subroutine deletes element pointed with posic

     TYPE (ele07), POINTER :: head, tail, anter, posic
     TYPE (ele07), POINTER :: e

     IF (.NOT.ASSOCIATED (anter)) THEN    !
       head => posic%next
     ELSE
       anter%next => posic%next
     END IF
     e => posic%next                       !keep pointer to next element
     IF( .NOT.ASSOCIATED(e) )tail => anter !last element in list
     DEALLOCATE (posic%dvolu,posic%cartd,posic%jacin,posic%stra0,posic%gamm0, &
                 posic%strsg,posic%qbar,posic%ambda,posic%lnods)
     IF(ASSOCIATED(posic%beta)) DEALLOCATE (posic%beta)
     IF(ASSOCIATED(posic%ehist))DEALLOCATE (posic%ehist,posic%strap)
     IF(ASSOCIATED(posic%stres))DEALLOCATE (posic%stres)
     IF(ASSOCIATED(posic%cd))DEALLOCATE (posic%cd)
     DEALLOCATE (posic)                    !deallocate fixed space
     posic => e                            !point to next element
     ! NULLIFY (posic,anter)
     RETURN
   END SUBROUTINE del_ele07e

   SUBROUTINE cut_ele07e (head, anter, posic)
     !This subroutine deletes a set pointed with posic
     ! without nullifying anter    ???? what for ????
     TYPE (ele07), POINTER :: head, anter, posic

     IF (.NOT.ASSOCIATED (anter)) THEN
       head => posic%next
     ELSE
       anter%next => posic%next
     ENDIF
     NULLIFY (posic)
   END SUBROUTINE cut_ele07e

   SUBROUTINE new_ele07e(elm)
   !Create a new element of ELE07 sets

     TYPE(ele07),POINTER:: elm

     ALLOCATE(elm)
     elm%numel = 0        !Initialize label of element
     elm%matno = 0        !     "     material number
     elm%angle = 0d0      !Initialize local angle
     NULLIFY(elm%lnods,elm%dvolu,elm%cartd,elm%jacin,elm%stra0,elm%gamm0,elm%qbar, &
             elm%ambda,elm%nangl,elm%jas,elm%ehist,elm%strap,elm%stres,elm%strsg,elm%beta)
     NULLIFY(elm%next)

   RETURN
   END SUBROUTINE new_ele07e

   INCLUDE 'acvdf7.fi'
   INCLUDE 'ap1tm7.fi'
   INCLUDE 'asstr7.fi'
   INCLUDE 'axep07.fi'
   INCLUDE 'bmat07.fi'
   INCLUDE 'bmatx7.fi'
   INCLUDE 'bmmt27.fi'
   INCLUDE 'bmmt37.fi'
   INCLUDE 'bshe07.fi'
   INCLUDE 'bshem7.fi'
   INCLUDE 'comp_ang07.fi'
   INCLUDE 'comp_angl7.fi'
   INCLUDE 'codel7.fi'
   INCLUDE 'commv7.fi'
   INCLUDE 'deriv6.fi'
   INCLUDE 'dump07.fi'
   INCLUDE 'elmda7.fi'
   INCLUDE 'expo07.fi'
   INCLUDE 'gauss7.fi'
   INCLUDE 'impo07.fi'
   INCLUDE 'intem7.fi'
   INCLUDE 'intr07.fi'
   INCLUDE 'intrf7.fi'
   !INCLUDE 'loadp7.fi'
   INCLUDE 'lumas7.fi'
   INCLUDE 'masel7.fi'
   INCLUDE 'nodxy7.fi'
   INCLUDE 'nods07.fi'
   INCLUDE 'outdy7.fi'
   INCLUDE 'rest07.fi'
   INCLUDE 'resvp7.fi'
   INCLUDE 'secd07.fi'
   INCLUDE 'setg07.fi'
   INCLUDE 'setga7.fi'
   INCLUDE 'stra07.fi'
   INCLUDE 'stran7.fi'
   INCLUDE 'surf07.fi'
   INCLUDE 'toar07.fi'
   INCLUDE 'updlo7.fi'

  END MODULE ele07_db
