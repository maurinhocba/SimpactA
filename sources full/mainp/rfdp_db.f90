 MODULE rfdp_db
   ! module to manage rotation-free slave-master constraints

   ! *** to obtain connectivities from element sets
   USE ele11_db, ONLY : ele11, point_ele11e, search_ele11
   USE ele13_db, ONLY : ele13, point_ele13e, search_ele13
   ! ***
   IMPLICIT NONE
   SAVE

   ! Derived type for a list containing dependent nodes

   ! list
   TYPE rfdp_nod
     INTEGER (kind=4) :: slave               !slave node (internal)
     REAL (kind=8), POINTER :: rfpar(:)      !local fixed coordinates
     INTEGER (kind=4), ALLOCATABLE :: lnods(:)   !master connectivities
     TYPE (rfdp_nod), POINTER :: next        !pointer to next element in the list
   END TYPE rfdp_nod

   ! Variables
   INTEGER (kind=4) :: nrfdp=0  !number of master-slave pairs
   LOGICAL :: planar            !if planar approach used
   TYPE (rfdp_nod), POINTER  :: rfd_head, rfd_tail

 CONTAINS

   SUBROUTINE ini_rfdp (head, tail)
     !initialize a dependent nodes list

     !Dummy arguments
     TYPE (rfdp_nod), POINTER :: head, tail

     NULLIFY (head, tail)
     nrfdp = 0
     RETURN

   END SUBROUTINE ini_rfdp

   SUBROUTINE add_rfdp (new, head, tail)
     !This subroutine adds data to the end of the list
     !Dummy arguments
     TYPE (rfdp_nod), POINTER :: new, head, tail

     !Check if a list is empty
     IF (.NOT. ASSOCIATED (head)) THEN
       !list is empty, start it
       head => new
       tail => new
       NULLIFY (tail%next)
       nrfdp = 1
     ELSE
       !add data to the list
       tail%next => new
       NULLIFY (new%next)
       tail => new
       nrfdp = nrfdp + 1

     END IF

   END SUBROUTINE add_rfdp

 END MODULE rfdp_db
