 SUBROUTINE locla8(nnode,axesc,heade,iwrit,coord)
 !**********************************************************************
 !
 !      READ local systems for beam elements
 !
 !**********************************************************************
 IMPLICIT NONE
 INTEGER (kind=4) :: nnode, & !number of nodes per element
                     axesc, & !local axes code
                     iwrit    !flag to write to disk
 REAL (kind=8) :: coord(3,*)  !global coordinates
 TYPE (ele08), POINTER :: heade  !pointer to first element

 INTEGER (kind=4) iel,i,j,k,rn(3)
 REAL (kind=8) a(3,3),aux,euler(9)
 REAL (kind=8),PARAMETER :: facto = 1.74532925199433D-02
 LOGICAL :: flag = .FALSE.
 TYPE(ele08), POINTER :: elem    !pointer to an element

 INTERFACE
   INCLUDE 'inrotm.h'
 END INTERFACE

 CALL listen('LOCLA8')           !read a line
 IF (exists('LOCALS')) THEN      !If local systems are to be read

   IF(iwrit == 1)WRITE(lures,"(/,'   Local cartesian systems for ',     &
                      &    'beam elements ')",ERR=9999)
   DO
     CALL listen('LOCLA8')       !read a line
     IF (exists('ENDLOC'))EXIT   !end of data
     iel = getint('IEL   ',0,'!Element number ...................')
     ! search element in list
     i = 1                !initializes
     elem => heade        !point to first
     DO
       IF(.NOT.ASSOCIATED(elem))THEN     !check
         WRITE(lures,"(' Error element No:',i5, ' not in list')",ERR=9999) iel
         STOP
       END IF
       IF( elem%numel == iel )EXIT       !element found
       elem => elem%next                 !point to next
     END DO
     IF( axesc <= 0)THEN                 !error
       WRITE(lures,"(' ERROR, Element',i7,' Does NOT expect Euler',     &
                     ' Angles for local system')",ERR=9999) iel
       CALL runend('LOCLA8: ELEMENT SYSTEM not EXPECTED')
     END IF

     DO j=1,axesc
       ! if euler angles are read, form local transformation matrix
       CALL rdfrre('LOCLA8',a(1:3,2),k,3,flag)
       IF(iwrit == 1)WRITE(lures,"(10x,'Angles',3f10.5)",ERR=9999) a(1:3,2)
       a(1:3,2) = a(1:3,2)*facto
       CALL inrotm (a(:,2),euler)
       elem%llbd(:,:,j) = RESHAPE(euler(1:9),(/3,3/))
     END DO

   END DO
 ELSE
   backs = .TRUE.
 END IF
 elem => heade
 DO
   IF( .NOT.ASSOCIATED(elem)) EXIT
   SELECT CASE (axesc)
   CASE (:-1)
     rn(1:-axesc) = INT(elem%llbd(1:-axesc,1,1))  !auxiliar nodes
     DO j=1,-axesc                             !loop over each node
       ! if auxiliar node used, form global transformation matrix
       a(1:3,2) = coord(1:3,rn(j))     !auxiliar node coordinates
       a(1:3,1) = coord(1:3,elem%lnods(nnode))-coord(1:3,elem%lnods(1))  !element axis
       CALL vecuni(3,a(1,1),aux)                      !Element direction
       a(1:3,3) = a(1:3,2) - coord(1:3,elem%lnods(1)) !direction in plane x1-x3
       aux = DOT_PRODUCT(a(1:3,1),a(1:3,3))           !proyection over element axix
       a(1:3,3) = a(1:3,3) - aux*a(1:3,1)             !x3 direction
       CALL vecuni(3,a(1,3),aux)                      !unit vector
       CALL vecpro(a(1,3),a(1,1),a(1,2))              !x2 direction
       elem%llbd(:,:,j) = a                           !assign local system
     END DO
   CASE (0)
     !nothing
   CASE (1:)
     !check all local system have been read
     IF( elem%llbd(1,1,1) < -1.5)THEN
       WRITE(lures,"(' ERROR, Element',i7,' expects Euler Angles')",ERR=9999) &
                      elem%numel
       CALL runend('LOCLA8: ELEMENT SYSTEM MISSING ....')
     END IF
   END SELECT
   elem => elem%next
 END DO

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE locla8
