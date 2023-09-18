 MODULE ele24_db
   USE param_db,ONLY: mnam,midn,mlin
   USE mat_dba, ONLY: section,sect_search,psecs,pmats,mater,curve,inte_cr,cur_point,postv,snn
   USE c_input
   USE kinc_db, ONLY : nndpd
   USE ctrl_db, ONLY : npoin, therm, ndoft, itemp, tref, tdtime
   USE npo_db, ONLY : iftmp, tempe,  dtemp
   IMPLICIT NONE

   ! Derived type for an ELE24 element
   ! Basic thin Shell Triangle (BST) for coupled thermo-mechanical analysis
   ! including Non-smooth surfaces and branching

   ! Reference: F.Flores & E.Oñate
   ! A rotation-free shell triangle for the analysis
   ! of kinked and branching shells
   ! IJNME, 69:1521-1551 (2007)

   SAVE

   TYPE sideb
     INTEGER (kind=4) :: nn
     INTEGER (kind=4), POINTER :: lnods(:)  !-1:nn = connectivities
     REAL (kind=8), POINTER  :: alph0(:)    !nn   = initial angles + side length
     REAL (kind=8), POINTER  :: gamma(:)    !nn   = distorsions
     REAL (kind=8), POINTER  :: fc(:,:)     !nn   = factors
     REAL (kind=8), POINTER  :: c(:,:)      !nn   = normal derivatives
     REAL (kind=8), POINTER  :: bb(:,:)     !ndof*(nn+2) = nodal B matrix
     TYPE (sideb), POINTER :: next
   END TYPE sideb

   TYPE pside
      TYPE (sideb), POINTER :: p
   END TYPE pside

   !   hh= side element connectivities
   INTEGER(kind=4), PARAMETER :: hh(3,3) = RESHAPE((/ 4,3,2, 5,1,3, 6,2,1 /), (/3,3/) )
!   REAL(kind=8), PARAMETER :: alp1 = 1.27d0, & !Max angle to apply membrane cuadratic approach
!                              alp2 = 1.57d0, & !angle to use standard CST for membrane
!                              alp3 = 0.3d0     !alp2-alp1 auxiliar value
   REAL(kind=8), PARAMETER :: alp1 = 0.7854d0, & !Max angle to apply membrane cuadratic approach
                              alp2 = 1.5708d0, & !angle to use standard CST for membrane
                              alp3 = 0.7854d0    !alp2-alp1 auxiliar value


   TYPE ele24
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
                      stra1(6),  & ! Actual mid-surface fundamental forms
                      pwork(3)     ! integrated plastic power
     REAL (kind=8), POINTER :: gausv(:,:)       !layer values
     LOGICAL :: mems(3)                         !use quadratic approach for that side
     TYPE (pside), POINTER :: si(:)
     TYPE (ele24), POINTER :: next              !pointer to next element
   END TYPE ele24

   ! Derived type for a set of ELE24 elements
   TYPE ele24_set
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
     LOGICAL :: origl           ! .TRUE. -> Original labels stored in connectivities
     INTEGER :: plstr           ! compute Plastic Strain Flag
         ! -1 from Cauchy stress  0 - do not   1 from 2nd P-K
     INTEGER :: locax           ! local x definition option
     REAL (kind=8) ::  angdf     ! angle between X_1 and orthotropic dir 1
     REAL (kind=8), POINTER :: stint(:,:)   !moments, forces and shears
     TYPE (ele24), POINTER    :: head, tail !pointer to first and last elm.
     TYPE (sideb), POINTER :: bhead , btail !pointers to branching data base
     INTEGER (kind=4), POINTER :: ngrqs(:)  !gauss points for output
     ! for shear evaluation
     INTEGER (kind=4) :: shear              !compute shear forces
     INTEGER (kind=4), POINTER :: ninv(:)   !inverse nodal relation (temporary, no dumping)
     REAL (kind=8), POINTER :: moments(:,:) !smoothed nodal moments (temporary, no dumping)
     REAL (kind=8), POINTER :: factors(:)   !nodal factors for smoothing (temporary, no dumping)
     TYPE (ele24_set), POINTER :: next      !pointer to next set
   END TYPE ele24_set
   TYPE (ele24_set), POINTER, SAVE :: head,tail !first and last elements sets

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


   SUBROUTINE new_ele24(elset)
   !Create a new element of ELE24 sets

     TYPE(ele24_set),POINTER:: elset

     ALLOCATE(elset)
     elset%sname = ''       !Initialize set name
     elset%nelem = 0        !     "     number of elements
     elset%nreqs = 0        !     "     number of GP for hist. output
     elset%narch = 0        !     "     number of output unit
     elset%nbs   = 0        !     "     number of branching sides
     elset%logst = .FALSE.  !     "     use logarithmic strain
     elset%lside = .FALSE.  !     "     flag to compute LSIDE
     elset%gauss = .FALSE.  !     "     flag to compute Gauss constants
     elset%shear = 0        !     "     flag to activate shear forces computation
     elset%origl = .FALSE.  !     "     flag to consider original labels
     elset%plstr = 0        !     "     compute Plastic Strain Flag
     elset%locax = 3        !     "     local x definition
     elset%angdf = 0d0      !     "     angle between X_1 and orthotropic dir 1
     NULLIFY(elset%head,elset%tail,elset%ngrqs,elset%stint,elset%bhead,elset%btail)
     NULLIFY(elset%moments,elset%factors,elset%ninv)
     NULLIFY(elset%next)

   RETURN
   END SUBROUTINE new_ele24

   SUBROUTINE ini_ele24 (head, tail)
     !initialize a list of ELE24 sets

     TYPE (ele24_set), POINTER :: head, tail

     NULLIFY (head, tail)

   END SUBROUTINE ini_ele24

   SUBROUTINE add_ele24 (new, head, tail)
     !This subroutine adds data to the end of the list
     !Dummy arguments
     TYPE (ele24_set), POINTER :: new, head, tail

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
   END SUBROUTINE add_ele24

   SUBROUTINE srch_ele24 (head, anter, posic, name, found)
     !This subroutine searches for a set named "name"
     !Dummy arguments
     LOGICAL :: found
     CHARACTER (len=*) :: name ! set name
     TYPE (ele24_set), POINTER :: head, anter, posic

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
   END SUBROUTINE srch_ele24

   SUBROUTINE del_ele24 (head, anter, posic)

     !This subroutine deletes a set pointed with posic

     TYPE (ele24_set), POINTER :: head, anter, posic

     TYPE (ele24), POINTER :: ea,ep
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
       CALL del_ele24e (posic%head,posic%tail, ea, ep )  !deletes element
     END DO

     NULLIFY (posic,anter)          !point to nothing
   END SUBROUTINE del_ele24

   ! ******* functions for a list of elements in a set ********

   SUBROUTINE ini_ele24e (head, tail)
     !initialize a list of ELE24 elements

     TYPE (ele24), POINTER :: head, tail

     NULLIFY (head, tail)       !initializes first and last pointer

   END SUBROUTINE ini_ele24e

   SUBROUTINE new_ele24e(elm)
   !Create a new element of ELE24 sets

     TYPE(ele24),POINTER:: elm

     ALLOCATE(elm)
     elm%numel = 0        !Initialize label of element
     elm%matno = 0        !     "     material number
     elm%lnods(1:6) = 0   !     "     conectivities
     elm%lside(1:3) = 0   !     "     neighbour elements
!     elm%ldv2sd(1:3) = 1  !     "     status of refined side
     elm%Area1 = 0d0      !     "     initial area
     elm%angle = 0d0      !     "     angle between dir 1 and orthotropic dir 1
     elm%lb = 0d0         !     "     initial area
     elm%a = 0d0          !     "     side projections in dir 1
     elm%b = 0d0          !     "     side projections in dir 2
     elm%c = 0d0          !     "     side projections in dir n
     elm%cd = 0d0         !     "     cartesian derivatives
     elm%a0 = 0d0         !     "     initial angles
     elm%gamma = 0d0      !     "     distorsions
     elm%stra1(1:6) = 0d0 !     "     actual mid-surface fundamental forms
     elm%mems = .FALSE.   !     "
     ALLOCATE( elm%si(3))
     NULLIFY(elm%gausv,elm%si(1)%p,elm%si(2)%p,elm%si(3)%p)
     elm%pwork(1:3) = 0d0 !     "     internal generated heat
     NULLIFY(elm%next)

   RETURN
   END SUBROUTINE new_ele24e

   SUBROUTINE add_ele24e (new, head, tail)
     !This subroutine adds data to the end of the list
     !Dummy arguments
     TYPE (ele24), POINTER :: new, head, tail

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
   END SUBROUTINE add_ele24e

   SUBROUTINE srch_ele24e (head, anter, posic, kelem, found)
     !This subroutine searches for an element labeled "kelem"
     !Dummy arguments
     LOGICAL :: found
     INTEGER (kind=4) :: kelem
     TYPE (ele24), POINTER :: head, anter, posic

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
   END SUBROUTINE srch_ele24e

   SUBROUTINE del_ele24e (head, tail, anter, posic)

     !This subroutine deletes element pointed with posic

     TYPE (ele24), POINTER :: head, tail, anter, posic

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
   END SUBROUTINE del_ele24e

   ! next two routines used by rotation-free nodal constraints
   SUBROUTINE point_ele24e (lnods, simple, kelem, sname, sfound, efound)
     ! point to a specific element in a set
     ! This subroutine searches for an element labeled "kelem" in the set SNAME

     !Dummy arguments
     INTEGER (kind=4) :: lnods(:)
     LOGICAL, OPTIONAL :: sfound,efound        ! flags
     CHARACTER (len=mnam), OPTIONAL :: sname   ! set name
     INTEGER (kind=4), OPTIONAL :: kelem
     LOGICAL, OPTIONAL :: simple
     !Local Variables
     TYPE (ele24_set), POINTER :: el_set ! elements set
     TYPE (ele24), POINTER :: posic
     INTEGER (kind=4) :: i,nset,iel

     el_set => head              !point to first set
     IF( PRESENT(kelem) ) THEN
       sfound = .FALSE.            !initializes set flag
       efound = .FALSE.            !initializes element flag
       nset = 0
       DO
         IF (.NOT.ASSOCIATED (el_set)) EXIT       !Check if a set exists
         nset = nset + 1
         IF( el_set%sname == sname )THEN
           sfound = .TRUE.
           EXIT
         END IF
         el_set => el_set%next
       END DO
       posic => el_set%head
       iel = 0
       DO
         IF (.NOT.ASSOCIATED (posic)) EXIT        !Check if a element exists
         iel = iel + 1
         IF( posic%numel == kelem )THEN
           efound = .TRUE.
           EXIT
         END IF
         posic => posic%next
       END DO
       IF( efound ) THEN
         IF( simple )THEN
           lnods(1:3) = posic%lnods(1:3)
         ELSE
           lnods(1:6) = (/ posic%lnods(1:3), -1, nset, iel /)
         END IF
       END IF
     ELSE
       DO i=1,lnods(5)-1
         el_set => el_set%next
       END DO
       posic => el_set%head
       DO i=1,lnods(6)-1
         posic => posic%next
       END DO
       lnods(4:6) =  posic%lnods(4:6)

     END IF
     RETURN
   END SUBROUTINE point_ele24e

   SUBROUTINE search_ele24 (lnods, sname, kelem, sfound, efound, sl, x, simple)
     ! point to a specific element in a set
     ! This subroutine searches for the nearest element to node SL in the set SNAME

     !Dummy arguments
     LOGICAL, INTENT(OUT) :: sfound,efound        ! flags
     CHARACTER (len=mnam), INTENT(IN) :: sname   ! set name
     INTEGER (kind=4), INTENT(IN) :: sl
     INTEGER (kind=4), INTENT(OUT) :: kelem,lnods(:)
     REAL(kind=8), INTENT(IN) :: x(:,:)
     LOGICAL, INTENT(IN) :: simple

     !Local Variables
     TYPE (ele24_set), POINTER :: el_set ! element set
     TYPE (ele24), POINTER :: new        ! element
     TYPE (ele24), POINTER :: posic
     REAL(kind=8) :: y(3,3),xc(3),xs(3),t1(3),t2(3),t3(3),dmin,d(3),l,s,t,ls,pmin
     INTEGER (kind=4) :: i,nset,iel,mnode
     REAL (kind=8), PARAMETER :: tol1 = 1.2d0, tol2 = -0.2 !acceptable proyection

     sfound = .FALSE.            !initializes set flag
     efound = .FALSE.            !initializes element flag

     mnode = kelem               !master node
     el_set => head              !point to first set
     nset = 0
     sets : DO
       IF (.NOT.ASSOCIATED (el_set)) EXIT sets     !Check if a set exists
       nset = nset + 1
       IF( el_set%sname == sname )THEN             !check set name
         sfound = .TRUE.                             !element set found
         new => el_set%head                          !point auxiliar to first element
         dmin = HUGE(dmin)                           !initializes minimal distance
         pmin = 1d0                                  !initializes minimal proyection
         xs = x(:,sl)                                !slave node coordinates
         DO i=1,el_set%nelem                           !loop over elements in the set
           IF( mnode > 0)THEN
             IF( ALL(new%lnods(1:3) /= mnode) )THEN
               new => new%next                                !point to next element
               CYCLE
             END IF
           END IF
           y = x(:,new%lnods(1:3))                     !element nodal coordinates
           !first check is with central point
           xc = (y(:,1)+y(:,2)+y(:,3))/3d0-xs
           ls = DOT_PRODUCT(xc,xc)
           IF( ls <= dmin )THEN
             dmin = 2d0*ls
             !second check includes area coordinates
             d  = xs - y(:,1)                         !distance vector to first node
             t1 = y(:,2) - y(:,1)                     !side vector 3
             t2 = y(:,3) - y(:,1)                     !-side vector 2
             CALL vecpro(t1(1),t2(1),t3(1))           !normal vector
             CALL vecuni(3,t3,ls)                     !unit normal vector and 2A
             l = DOT_PRODUCT(d,t3)                    !distance along normal
             d = d - l*t3                             !in plane proyection
             CALL vecpro(t1(1),d(1),xc(1))            !normal vector
             s = DOT_PRODUCT(xc,t3)/ls
             IF( s < tol1 .AND. s > tol2)THEN           !is proyection acceptable?
               CALL vecpro(d(1),t2(1),xc(1))            !normal vector
               t = DOT_PRODUCT(xc,t3)/ls
               IF( t < tol1 .AND. t > tol2)THEN             !is proyection acceptable?
                 l = 1d0-s-t
                 IF( l < tol1 .AND. l > tol2)THEN           !is proyection acceptable?
                   ls = 0d0                                    !initializes proyection
                   IF( s < 0d0 )ls = ls-s
                   IF( s > 1d0 )ls = ls+s-1d0
                   IF( t < 0d0 )ls = ls-t
                   IF( t > 1d0 )ls = ls+s-1d0
                   IF( l < 0d0 )ls = ls-l
                   IF( l > 1d0 )ls = ls+l-1d0
                   IF( l < pmin )THEN
                     efound = .TRUE.                            !element found
                     pmin = l                                   !update minimal proyection error
                     posic => new                               !keep element pointer
                     kelem = posic%numel                        !keep element label
                     iel = i
                     IF( l == 0d0 )   EXIT
                   END IF
                 END IF
               END IF
             END IF
           END IF
           new => new%next                                !point to next element
         END DO
         EXIT sets
       END IF
       el_set => el_set%next                         !point to next set
     END DO sets
     IF( efound ) THEN
       IF( simple )THEN
         lnods(1:3) = posic%lnods(1:3)
       ELSE
         lnods(1:6) = (/ posic%lnods(1:3), -1, nset, iel /)
       END IF
     END IF
     RETURN
   END SUBROUTINE search_ele24

    INCLUDE 'acvd24.fi'
    INCLUDE 'axep24.fi'
    INCLUDE 'bbra24.fi'
    INCLUDE 'bfle24.fi'
    INCLUDE 'bmem24.fi'
    INCLUDE 'boun24.fi'
    INCLUDE 'bran24.fi'
    INCLUDE 'code24.fi'
    INCLUDE 'comm24.fi'
    INCLUDE 'delt24.fi'
    INCLUDE 'dump24.fi'
    INCLUDE 'elmd24.fi'
    INCLUDE 'expo24.fi'
    INCLUDE 'gaus24.fi'
    INCLUDE 'impo24.fi'
    INCLUDE 'inig24.fi'
    INCLUDE 'luma24.fi'
    INCLUDE 'mase24.fi'
    INCLUDE 'nods24.fi'
    INCLUDE 'outd24.fi'
    INCLUDE 'rest24.fi'
    INCLUDE 'resv24.fi'
    INCLUDE 'secd24.fi'
    INCLUDE 'slum24.fi'
    INCLUDE 'stra24.fi'
    INCLUDE 'surf24.fi'
    INCLUDE 'tlma24.fi'
    INCLUDE 'toar24.fi'
    INCLUDE 'tres24.fi'
    INCLUDE 'updl24.fi'

 END MODULE ele24_db
