 MODULE nes_db

   IMPLICIT NONE
   SAVE

   ! Derived type for a list containing Multi-Point-Constraint definition

   ! list
   TYPE nes_nod
     INTEGER (kind=4) :: node,ndof       !master (or SLAVE) node and DOF
     REAL (kind=8) :: factor             !factor (if 0d0 is a slave)
     TYPE (nes_nod), POINTER :: next     !pointer to next
   END TYPE nes_nod
   ! array
   TYPE nes_data
     INTEGER (kind=4) :: node,ndof       !master (or SLAVE) node and DOF
     REAL (kind=8) :: factor             !factor (if 0d0 is a slave)
   END TYPE nes_data
   TYPE (nes_data), POINTER :: nesdat(:)
   ! factor = 0d0 => indicates the SLAVE node
   INTEGER (kind=4) :: nnes=0            !number of data in the array

 CONTAINS

   SUBROUTINE ini_nes (head, tail)
     !initialize a MultiPointConstraint list
     !a unique list for all the data

     !Dummy arguments
     TYPE (nes_nod), POINTER :: head, tail

     NULLIFY (head, tail)

   END SUBROUTINE ini_nes

   SUBROUTINE add_nes (new, head, tail)
     !This subroutine adds a MPC data to the end of the list
     !Dummy arguments
     TYPE (nes_nod), POINTER :: new, head, tail

     !Check if a list is empty
     IF (.NOT. ASSOCIATED (head)) THEN
       !list is empty, start it
       head => new
       tail => new
       NULLIFY (tail%next)

     ELSE
       !add data to the list
       tail%next => new
       NULLIFY (new%next)
       tail => new

     ENDIF
   END SUBROUTINE add_nes

   SUBROUTINE store_nes (head)
     !This subroutine writes the MPC list into the array nesdat
     ! and release memory associated to list

     !Dummy argument
     TYPE (nes_nod), POINTER :: head  !begining of the list

     !Local variables
     TYPE (nes_nod), POINTER :: ptr,prev
     INTEGER :: n

     ptr => head                !point to first
     n = 0                      !initializes counter
     DO
       IF (.NOT.ASSOCIATED(ptr)) EXIT  !exit if at end of the list found
       n = n + 1                       !increment counter
       nesdat(n)%node   = ptr%node     !copy data in list into the array
       nesdat(n)%ndof   = ptr%ndof
       nesdat(n)%factor = ptr%factor
       prev => ptr
       ptr => ptr%next                 !point to next
       DEALLOCATE (prev)               !release memory
     END DO
   END SUBROUTINE store_nes

   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   SUBROUTINE reasd0 (iwrit,actio)

   !This routine reads Multi Point Constraint definitions

   USE param_db,ONLY: mnam
   USE c_input
   USE nsets_db
   IMPLICIT NONE

   CHARACTER(len=*),INTENT(INOUT):: actio
   INTEGER(kind=4):: iwrit

   ! Local
   TYPE (nes_nod), POINTER :: mpcs,head,tail,first
   TYPE (nes_data), POINTER :: ms_list(:)
   TYPE (nset), POINTER :: ns
   INTEGER :: i,j,n,idof,nms
   LOGICAL ::  aset, slave, firsn, found
   CHARACTER(len=mnam) :: sname

   IF (iwrit == 1) WRITE(lures,"(//,5x,'***** SLAVE D.O.F. DATA *****',/)",ERR=9)

   !Initialize empty list
   CALL ini_nes(head,tail)
   IF (TRIM(actio) == 'NEW' .OR. .NOT.exists('ADD   ')) THEN
     nnes=0                          !initializes number of slave DOFs
   ELSE IF (exists('ADD   ')) THEN
     IF (iwrit == 1) WRITE(lures,"(//' Previous NESCV definitions are kept')",ERR=9)
     DO n=1,nnes                         !for existing slave DOFs
       ALLOCATE (mpcs)                   !get memory
       mpcs%node = nesdat(n)%node        !master (slave) node
       mpcs%ndof = nesdat(n)%ndof        !master (slave) DOF
       mpcs%factor = nesdat(n)%factor    !factor (0 indicates a slave)
       CALL add_nes (mpcs, head, tail)   !add to list
     END DO
   END IF

   IF (TRIM(actio) == 'NSTRA0') actio='NSTRA1'
   IF (ASSOCIATED (nesdat) ) DEALLOCATE (nesdat) !release memory

   IF (iwrit == 1) WRITE(lures,"(//'   Slave  DOF  Master  DOF  Factor')",ERR=9)

   !Loop to read constraint data and add them to the list
   CALL listen('REASD0')      !read a line
   DO
     IF (exists('ENDNES')) EXIT     !end of MPC data => EXIT
     !    Read a new MPC
     IF(exists('CONNOD',i))THEN     !IF a new slave node read
       ALLOCATE (mpcs)              !get memory
       mpcs%node = INT(param(i))    !associated node
       IF(exists('CONDOF',j))THEN   !associated DOF
         mpcs%ndof = INT(param(j))
         idof = mpcs%ndof           !keep slave DOF
       ELSE
         CALL runend('REASD0: CONDOF is compulsory  .... ')
       END IF
       mpcs%factor = 0d0            ! indicator of a constrained node
       nnes=nnes+1                  ! increase counter
       CALL add_nes( mpcs, head, tail )    ! add to list
       IF( mpcs%node /= 0 )THEN
         IF (iwrit == 1) WRITE(lures, '(i8,i5)',ERR=9) mpcs%node,mpcs%ndof
         aset = .FALSE.
       ELSE
         IF (TRIM(words(i+1)(1:midn)) == 'SET') THEN
           aset = .TRUE.
           sname = get_name(posin=i+1,stype='NSET')      !set name
           CALL nsdb_search (sname, found, ns)   !search if node set exists
           IF (found) THEN              !set found, go ahead
             IF (iwrit == 1) WRITE(lures,"('SET ',A,i5)",ERR=9)  &
               TRIM(sname),mpcs%ndof
             CALL ns_head (ns)          !point to first value in the set
           ELSE          !set not found => ERROR
             WRITE (lures, "(' Set ',A,' not found')",ERR=9) TRIM(sname)
             CALL runend('Intime: Set not found ')
           END IF
           first => tail
         END IF
       END IF
     ELSE
       CALL runend('REASD0: CONNOD expected .......... ')
     END IF
     !     Read the dependencies
     IF(.NOT.exists('INDNOD',i)) THEN !check if first line include information
       CALL listen('REASD0')     !read a line
       IF(.NOT.exists('INDNOD',i))CALL runend('REASD0: INDNOD expected .......... ')
     END IF
     nms = 0
     DO
       ALLOCATE (mpcs)                           !get memory
       mpcs%node = INT(param(i))                 !Master node
       IF(exists('INDDOF',i))THEN                !Master DOF
         mpcs%ndof = INT(param(i))
       ELSE
         mpcs%ndof = idof                        !default dof
       END IF
       IF(exists('FACTOR',i))THEN                !factor
         IF(param(i) == 0d0) CALL runend('REASD0: FACTOR must not be 0 ..... ')
         mpcs%factor =  param(i)
       ELSE
         mpcs%factor = 1d0                       !default value
       END IF
       nnes=nnes+1                               !increase counter
       nms = nms + 1                             !number of master dofs

       IF (iwrit == 1) WRITE(lures, '(13x,i8,i5,e12.4)',ERR=9) &
                       &       mpcs%node,mpcs%ndof,mpcs%factor

       CALL add_nes( mpcs, head, tail )          !add to list
       CALL listen('REASD0')                     !read new line
       IF(.NOT.exists('INDNOD',i).OR.exists('CONNOD'))EXIT           !exit loop
     END DO  !end of dependencies for a DOF

     !*****  Add data for a set of nodes
     IF( aset )THEN             !if slave data is associated to a set
       ALLOCATE( ms_list(nms) ) !get memory for the list of master nodes
       mpcs => first%next       !point to the first master value
       DO i=1,nms               !for each master value
         ms_list(i)%node   = mpcs%node  !copy data in list into the array
         ms_list(i)%ndof   = mpcs%ndof
         ms_list(i)%factor = mpcs%factor
         mpcs => mpcs%next
       END DO
       firsn = .TRUE.             !initializes to first node in the set
       DO                         !loop over the nodes in the set
         IF ( end_ns(ns) ) EXIT   !last node processed, EXIT
         n = get_label(ns)        !get nodal label
         slave = .TRUE.           !initializes
         DO i=1,nms               !check if the node is defined master
           IF(    n == ms_list(i)%node .AND.   &
               idof == ms_list(i)%ndof )THEN !if node is the master node list
             slave = .FALSE.      !this node must be skipped
             EXIT
           END IF
         END DO
         IF( slave )THEN          !if node is slave
           IF( firsn )THEN        !if it is the first node
             first%node   = n     !copy node data in first position
             firsn = .FALSE.      !not first node any longer
           ELSE
             ALLOCATE (mpcs)      !get memory
             mpcs%node   = n      !slave node
             mpcs%ndof   = idof   !slave dof
             mpcs%factor = 0d0    ! indicator of a constrained node
             nnes=nnes+1          ! increase counter
             CALL add_nes( mpcs, head, tail )    ! add to list
             DO i=1,nms           !for each master node
               ALLOCATE (mpcs)    !get memory
               mpcs%node   = ms_list(i)%node   !copy data in list
               mpcs%ndof   = ms_list(i)%ndof
               mpcs%factor = ms_list(i)%factor
               nnes=nnes+1                  ! increase counter
               CALL add_nes( mpcs, head, tail )    ! add to list
             END DO               !list of master nodes
           END IF
         END IF
         CALL ns_next (ns)        !point to next node in the list
       END DO
       DEALLOCATE(ms_list)
     END IF
     !*****
   END DO    !end of constraints

   !Store MPC data in an array for easy reference
   ALLOCATE (nesdat(nnes))    !get memory for the array
   CALL store_nes(head)       !transfer data from list to array
                              ! and release memory
   CALL ini_nes(head,tail)    !nullify pointers

   RETURN
   9 CALL runen2(' error while writing data to the disk')
   END SUBROUTINE reasd0

 END MODULE nes_db
