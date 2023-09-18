 MODULE ele10_db
   USE param_db,ONLY : mnam, midn
   USE mat_dba, ONLY: section,sect_search,psecs,pmats,mater,curve,inte_cr,cur_point,postv
   USE c_input
   IMPLICIT NONE
   SAVE

   ! Derived type for a set of RIGID elements
   TYPE ele10_set
     CHARACTER (len=mnam) :: sname ! set name
     INTEGER (kind=4) :: nelem   ! number of elements
     INTEGER (kind=4) :: nnode   ! maximum number of nodes per element 2D=4, 3D=8
     INTEGER (kind=4) :: rbnod   ! label of master node
     INTEGER (kind=4) :: nmast   ! master node (internal number)
     REAL (kind=8) :: tmass                  !Total translational mass
     INTEGER (kind=4) :: ntype   ! problem type
             ! 2-D 0 rigid body defined by particles
             !     1 plane solid with thickness
             !     2 plane solid with thickness = 1
             !     3 axilsymmetric solid
             !     4 surface with transversal thickness
             !     2 surface with transversal thickness = 1
             !     3 axilsymmetric surface
             ! 3-D 0 rigid body defined by particles
             !     1-3 solid
             !     4-6 surface
     INTEGER (kind=4), POINTER :: lnods(:,:) !connectivities
     INTEGER (kind=4), POINTER :: matno(:)   !associated material
     ! variable exclusively for thermal problems
     LOGICAL :: heat = .FALSE.               !.TRUE. if conduction is considered
     INTEGER (kind=4) :: ngaus   ! number of integration points for thermal problems
     REAL (kind=8), POINTER :: cartd(:,:,:,:)   !cartesian derivatives of shape functions
     REAL (kind=8), POINTER :: dvolu(:,:)       !gauss point associated volume
     REAL (kind=8), POINTER :: shape(:,:)       !gauss point shape functions
     !
     TYPE (ele10_set), POINTER :: next       !pointer to next set

   END TYPE ele10_set
   TYPE (ele10_set), POINTER :: head
   TYPE (ele10_set), POINTER :: tail

 CONTAINS
   SUBROUTINE ini_ele10 (head, tail)
     !initialize a list of RIGID sets

     TYPE (ele10_set), POINTER :: head, tail

     NULLIFY (head, tail)

   END SUBROUTINE ini_ele10

   SUBROUTINE add_ele10 (new, head, tail)
     !This subroutine adds data to the end of the list
     !Dummy arguments
     TYPE (ele10_set), POINTER :: new, head, tail

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
   END SUBROUTINE add_ele10

   SUBROUTINE srch_ele10 (head, anter, posic, name, found)
     !This subroutine searches for a set named "name"
     !Dummy arguments
     LOGICAL :: found
     CHARACTER (len=*) :: name ! set name
     TYPE (ele10_set), POINTER :: head, anter, posic

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
   END SUBROUTINE srch_ele10

   SUBROUTINE del_ele10 (head, tail, anter, posic)
     !This subroutine deletes a set pointed with posic
     TYPE (ele10_set), POINTER :: head, tail, anter, posic

     IF (.NOT.ASSOCIATED (anter)) THEN
         IF( ASSOCIATED(posic%next) )THEN
           head => posic%next
         ELSE
           NULLIFY( head )
         END IF
     ELSE
       anter%next => posic%next
     ENDIF
     ! if posic == tail
     IF (.NOT.ASSOCIATED (posic%next) ) tail => anter
     CALL dalloc10 (posic)
     !NULLIFY (anter)
   END SUBROUTINE del_ele10

   SUBROUTINE dalloc10 (elset)
     ! deallocates a set
     TYPE (ele10_set), POINTER :: elset
     DEALLOCATE ( elset%lnods )
     IF(ASSOCIATED (elset%matno) ) DEALLOCATE ( elset%matno)
     IF(elset%heat) DEALLOCATE ( elset%shape, elset%cartd, elset%dvolu )
     DEALLOCATE (elset)
   END SUBROUTINE dalloc10

   INCLUDE 'acvd10.fi'
   INCLUDE 'code10.fi'
   INCLUDE 'comm10.fi'
   INCLUDE 'dump10.fi'
   INCLUDE 'elmd10.fi'
   INCLUDE 'gaus10.fi'
   INCLUDE 'luma10.fi'
   INCLUDE 'mase10.fi'
   INCLUDE 'poin10.fi'
   INCLUDE 'rearb0.fi'
   INCLUDE 'rest10.fi'
   INCLUDE 'rigb10.fi'
   INCLUDE 'rigbdy.fi'
   INCLUDE 'surf10.fi'
   INCLUDE 'tlma10.fi'
   INCLUDE 'tres10.fi'
   INCLUDE 'updl10.fi'

 END MODULE ele10_db
