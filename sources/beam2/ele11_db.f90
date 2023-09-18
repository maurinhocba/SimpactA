  MODULE ele11_db
    USE param_db,ONLY: mnam,midn
    USE mat_dba, ONLY : section,mater,sect_search,psecs,inte_cr
    USE c_input
    USE kinc_db, ONLY : nndpd
    IMPLICIT NONE
    INTEGER (kind=4), PARAMETER :: &
      ndofe = 2 , &   ! number of DOFs per node
      nnode = 4       ! number of nodes per element
    INTEGER (kind=4), PARAMETER :: nxn(2:3) = (/1,4/),nyn(2:3) = (/3,2/) !cycling list

    SAVE

   TYPE nodei
     INTEGER (kind=4) :: nn
     INTEGER (kind=4), POINTER :: lnods(:)  !0:nn = connectivities
     REAL (kind=8) :: ba                    !average rotation
     REAL (kind=8), POINTER  :: alph0(:)    !nn   = initial angles
     REAL (kind=8), POINTER  :: fc(:)       !nn   = factors
     REAL (kind=8), POINTER  :: bb(:,:)     !ndof*(nn+1) = nodal B matrix
     TYPE (nodei), POINTER :: next
   END TYPE nodei

   TYPE ele11
     INTEGER (kind=4) :: numel  ! label of element
     INTEGER (kind=4) :: matno  ! section number
     INTEGER (kind=4) :: lnods(nnode+2)  ! Conectivities  nodes + neighbor elements
     REAL (kind=8)  :: l0(3),     & ! length
                       cab(4,2),  & ! a0, li, gamma, c
                       lambd(3)     ! stretch ratios
     REAL (kind=8), POINTER ::     & !
                       ro(:,:),    & ! initial radius and curvatures at Gauss points
                       cur(:),     & ! present meridian curvatures (both nodes)
                       gausv(:,:,:)  !(3,nlayr,ngaus)  internal variables (plasticity variables)
     TYPE (ele11), POINTER :: next              !pointer to next element
   END TYPE ele11

    ! Derived type for a set of BEAM2 elements
    TYPE ele11_set
      CHARACTER (len=mnam) :: sname   ! set name
      INTEGER (kind=4) ::   &
                          nelem, & ! number of elements
                          ngaus, & ! number of stresses/GP
                          nstre, & ! number of stresses/GP
                          nreqs, & ! number of GP for hist. output
                          nbn,   & ! number of branching nodes
                          narch    ! number of output unit

      LOGICAL :: lside, &        !.TRUE. if side nodes already defined
                 strai, &        !.TRUE. if initial curvatures are null
                 gauss           ! .FALSE. -> Initial constants not defined or not updated
      INTEGER :: plstr           ! compute Plastic Strain Flag
         ! -1 from Cauchy stress  0 - do not   1 from 2nd P-K
      REAL (kind=8) :: shap(2,2)   ! shape functions
      REAL (kind=8) :: stabs       ! stabilization factor


      TYPE (ele11), POINTER :: head, tail    !pointer to first and last elm.
      INTEGER (kind=4), POINTER :: ngrqs(:)  !(nreqs) Points for history output
      REAL (kind=8), POINTER :: stint(:,:)   !(nstre*ngaus,nelem) present forces and moments
      TYPE (nodei), POINTER :: nhead , ntail !pointers to branching data base
      TYPE (ele11_set), POINTER :: next
    END TYPE ele11_set

    TYPE (ele11_set), POINTER :: head => null(), tail => null() !first and last elements sets

  CONTAINS

    FUNCTION atan3(a,b,p)
    IMPLICIT NONE
    REAL(kind=8) :: atan3
    REAL(kind=8), INTENT(IN) :: a,b,p
    REAL (kind=8),PARAMETER :: twopi=6.283185307179586d0
      atan3 = ATAN2(a,b) - p  !angle change
      !  limit angle change to Pi (180 degrees)
      IF( atan3 > 3.14d0 )atan3 = atan3 - twopi
      IF( atan3 <-3.14d0 )atan3 = atan3 + twopi
    END FUNCTION

   SUBROUTINE srch_ele11 (head, posic, name, found)
     !This subroutine searches for a set named "name"
     !Dummy arguments
     LOGICAL, INTENT(OUT) :: found
     CHARACTER (len=*), INTENT(IN) :: name ! set name
     TYPE (ele11_set), POINTER, INTENT(IN) :: head
     TYPE (ele11_set), POINTER, INTENT(OUT) :: posic

     found = .FALSE.
     NULLIFY (posic)
     !Check if a list is empty
     IF (ASSOCIATED (head)) THEN
       posic => head
       DO
         IF(TRIM(posic%sname) == TRIM(name)) THEN
           found = .TRUE.
           EXIT
         ELSE IF (ASSOCIATED(posic%next) ) THEN
           posic => posic%next
         ELSE
           NULLIFY (posic)
           EXIT
         END IF
       END DO
     END IF

   END SUBROUTINE srch_ele11

   SUBROUTINE new_ele11(elset)

   !Create a new element of ELE11 sets

 !$ USE omp_lib
     TYPE(ele11_set),POINTER :: elset
 !$ INTEGER (kind=4) :: i
 !$ INTEGER(kind=4),PARAMETER:: MAXLOCK=2048              !number of keys
 !$ INTEGER(kind=4),SAVE     :: lock1(0:MAXLOCK-1)        !keys
 !$ LOGICAL,SAVE             :: inilock=.FALSE.           !initialization flag

 !Initialize exclusion variables for parallel code (only first time)
 !$ IF (.NOT.inilock) THEN                   !first time only
 !$   print *,'first time'
 !$   DO i=0,MAXLOCK-1                     !for each key
 !$     PRINT *,i,lock1(i)
 !$     CALL omp_init_lock(lock1(i))          !unlock
 !$     PRINT *,i+1,lock1(i)
 !$   END DO
 !$   inilock = .TRUE.                         !change flag
 !$   print *,'first time ok'
 !$ END IF


     ALLOCATE(elset)
     elset%sname = ''       !Initialize set name
     elset%nelem = 0        !     "     number of elements
     elset%ngaus = 2        !     "     number of GaussPoints
     elset%nstre = 0        !     "     number of generalized stresses
     elset%nreqs = 0        !     "     number of GP for hist. output
     elset%nbn   = 0        !     "     number of branching nodes
     elset%narch = 0        !     "     number of output unit
     elset%strai = .FALSE.  !     "     initial curvatures flag
     elset%lside = .FALSE.  !     "     flag to compute LSIDE
     elset%gauss = .FALSE.  !     "     flag to compute Gauss constants
     elset%plstr = 0        !     "     compute Plastic Strain Flag
     elset%shap  = 0d0
     elset%stabs = 0.1d0
     NULLIFY(elset%head,elset%tail,elset%ngrqs,elset%stint,elset%nhead,elset%ntail)
     NULLIFY(elset%next)

   RETURN
   END SUBROUTINE new_ele11

    SUBROUTINE add_ele11 (new, head, tail)
      !This subroutine adds data to the end of the list
      !Dummy arguments
      TYPE (ele11_set), POINTER :: new, head, tail

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
    END SUBROUTINE add_ele11

    SUBROUTINE del_ele11 (head, anter, posic)
      !This subroutine deletes a set pointed with posic
      TYPE (ele11_set), POINTER :: head, anter, posic

      TYPE (ele11), POINTER :: ea,ep
      TYPE (nodei), POINTER :: nb,nba     !pointers to branching node
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
        CALL del_ele11e (posic%head,posic%tail, ea, ep )  !deletes element
      END DO
      ! deallocate branching nodes database
      IF( posic%nbn > 0 )THEN
        nba => posic%nhead
        DO iel=1,posic%nbn
          nb => nba
          DEALLOCATE(nb%lnods,nb%alph0,nb%fc,nb%bb,nb)
          nba => nb%next
          DEALLOCATE(nb)
        END DO
      END IF
     NULLIFY (posic,anter)          !point to nothing
    END SUBROUTINE del_ele11

   ! ******* functions for a list of elements in a set ********

   SUBROUTINE new_ele11e(elm)
   !Create a new element of ELE11 type

     TYPE(ele11),POINTER:: elm
     ALLOCATE(elm)
     elm%numel = 0        !Initialize label of element
     elm%matno = 0        !     "     material number
     elm%lnods = 0        !     "     conectivities
     elm%l0    = 0d0      !     "     initial lengths
     elm%cab   = 0d0      !     "     auxiliar values
     elm%lambd = 1d0      !     "     stretch ratios
     NULLIFY(elm%ro)
     NULLIFY(elm%cur)
     NULLIFY(elm%gausv)
     NULLIFY(elm%next)

   RETURN
   END SUBROUTINE new_ele11e

   SUBROUTINE add_ele11e (new, head, tail)
     !This subroutine adds data to the end of the list
     !Dummy arguments
     TYPE (ele11), POINTER :: new, head, tail

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
   END SUBROUTINE add_ele11e

   SUBROUTINE srch_ele11e (head, anter, posic, kelem, found)
     !This subroutine searches for an element labeled "kelem"
     !Dummy arguments
     LOGICAL :: found
     INTEGER (kind=4) :: kelem
     TYPE (ele11), POINTER :: head, anter, posic

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
   END SUBROUTINE srch_ele11e

   SUBROUTINE del_ele11e (head, tail, anter, posic)

     !This subroutine deletes element pointed with posic

     TYPE (ele11), POINTER :: head, tail, anter, posic

     IF (.NOT.ASSOCIATED (anter)) THEN    !
       head => posic%next
     ELSE
       anter%next => posic%next
     END IF
     IF( .NOT.ASSOCIATED(posic%next) )tail => anter !last element in list
     DEALLOCATE(posic%gausv,posic%ro,posic%cur)
     DEALLOCATE (posic)                    !deallocate fixed space
     IF( ASSOCIATED( anter) )posic => anter%next                   !point to next element
     RETURN
   END SUBROUTINE del_ele11e

   ! next two routines used by rotation-free nodal constraints
   SUBROUTINE point_ele11e (lnods, simple, kelem, sname, sfound, efound)
     ! point to a specific element in a set
     ! This subroutine searches for an element labeled "kelem" in the set SNAME

     !Dummy arguments
     INTEGER (kind=4), INTENT(IN OUT) :: lnods(:)

     LOGICAL, INTENT(IN), OPTIONAL :: simple
     LOGICAL, OPTIONAL :: sfound,efound        ! flags
     CHARACTER (len=mnam), OPTIONAL :: sname   ! set name
     INTEGER (kind=4), OPTIONAL :: kelem

     !Local Variables
     TYPE (ele11_set), POINTER :: el_set ! elements set
     TYPE (ele11), POINTER :: posic
     INTEGER (kind=4) :: nset,iel

     el_set => head              !point to first set

     IF( PRESENT(kelem) )THEN
       sfound = .FALSE.            !initializes set flag
       efound = .FALSE.            !initializes element flag
       nset = 0
       DO
         IF (.NOT.ASSOCIATED (el_set)) EXIT       !Check if a set exists
         nset = nset+1
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
           IF( simple )THEN
             lnods(1:2) = posic%lnods(2:3)
           ELSE
             lnods(1:3) = (/-1, nset, iel /)
           END IF
           EXIT
         END IF
         posic => posic%next
       END DO
     ELSE
       DO nset=1,lnods(2)-1
         el_set => el_set%next
       END DO
       posic => el_set%head
       DO iel=1,lnods(3)-1
         posic => posic%next
       END DO
       lnods(1:4) = posic%lnods(1:4)
     END IF
     RETURN
   END SUBROUTINE point_ele11e

   SUBROUTINE search_ele11 (lnods, sname, kelem, sfound, efound, sl, x, simple)
     !
     ! This subroutine searches for the nearest element to node SL in the set SNAME

     !Dummy arguments
     LOGICAL, INTENT(OUT) :: sfound,efound       ! flags
     CHARACTER (len=mnam), INTENT(IN) :: sname   ! set name
     INTEGER (kind=4), INTENT(IN) :: sl
     INTEGER (kind=4), INTENT(OUT) :: kelem,lnods(:)
     REAL(kind=8), POINTER :: x(:,:)
     LOGICAL, INTENT(IN) :: simple

     !Local Variables
     TYPE (ele11_set), POINTER :: el_set ! element set
     TYPE (ele11), POINTER :: new        ! element
     TYPE (ele11), POINTER :: posic
     REAL(kind=8) :: y(2,2),xc(2),xs(2),t(2),dmin,d(2),l,s,ls
     INTEGER (kind=4) :: i,nset,iel,mnode
     REAL (kind=8), PARAMETER :: toler = 1.05d0 !acceptable proyection

     sfound = .FALSE.            !initializes set flag
     efound = .FALSE.            !initializes element flag
     mnode = kelem

     el_set => head              !point to first set
     nset = 0
     sets : DO
       IF (.NOT.ASSOCIATED (el_set)) EXIT sets      !Check if a set exists
       nset = nset + 1
       IF( el_set%sname == sname )THEN              !check set name
         sfound = .TRUE.                               !element set found
         new => el_set%head                            !point auxiliar to first element
         dmin = HUGE(dmin)                             !initializes minimal distance
         xs = x(:,sl)                              !slave node coordinates

         DO i=1,el_set%nelem                              !loop over elements in the set
           IF( mnode > 0 )THEN                            !if a master node is given
             IF( ALL(new%lnods(2:3) /= mnode) )THEN         !if master node is not in element
               new => new%next                                !point to next element
               CYCLE                                          !and go to next element
             END IF
           END IF
           y = x(:,new%lnods(2:3))                    !element nodal coordinates
           xc = (y(:,1)+y(:,2))/2d0                       !element center coordinates
           d  = xs - xc                                   !distance vector to element center
           l  = SQRT(DOT_PRODUCT(d,d))                    !distance
           IF( l < dmin )THEN                             !compare with minimal distance
             t = y(:,2) - y(:,1)                            !tangent vector
             CALL vecuni(2,t,ls)                            !unit tangent vector and length
             s = ABS(2d0*DOT_PRODUCT(d,t)/ls)               !normalized coordinate
             IF( s < toler )THEN                            !is proyection acceptable?
               efound = .TRUE.                                !element found
               dmin = ls                                      !update minimal distance
               posic => new                                   !keep element pointer
               kelem = posic%numel                            !keep element label
               iel = i
             END IF
           END IF
           new => new%next                                !point to next element
         END DO
         EXIT sets
       END IF
       el_set => el_set%next                         !point to next set
     END DO sets
     IF( efound )THEN
       IF( simple )THEN
         lnods(1:2) = posic%lnods(2:3)
       ELSE
         lnods(1:3) = (/ -1, nset, iel /)
       END IF
     END IF
     RETURN
   END SUBROUTINE search_ele11

   INCLUDE 'acvd11.fi'
   INCLUDE 'bfle11.fi'
   INCLUDE 'bmem11.fi'
   INCLUDE 'comm11.fi'
   INCLUDE 'delt11.fi'
   INCLUDE 'dump11.fi'
   INCLUDE 'gaus11.fi'
   INCLUDE 'elmd11.fi'
   INCLUDE 'luma11.fi'
   INCLUDE 'mase11.fi'
   INCLUDE 'outd11.fi'
   INCLUDE 'rest11.fi'
   INCLUDE 'resv11.fi'
   INCLUDE 'stra11.fi'
   INCLUDE 'streb2.fi'
   INCLUDE 'surf11.fi'
   INCLUDE 'toar11.fi'
   INCLUDE 'updl11.fi'

  END MODULE ele11_db
