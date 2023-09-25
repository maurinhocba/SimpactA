 SUBROUTINE readp0 (iwrit,actio)

 !This routine reads dependent nodes (nndp)

 USE param_db,ONLY: mnam
 USE c_input
 USE ndp_db
 USE nsets_db
 IMPLICIT NONE
 CHARACTER(len=*),INTENT(INOUT):: actio
 INTEGER :: iwrit

 ! Local
 TYPE (ndp_nod), POINTER :: ndp,head,tail
 TYPE (nset), POINTER :: ns
 INTEGER :: i,j,n,slave,master
 LOGICAL :: dummy,found
 CHARACTER(len=mnam) :: sname

 IF (iwrit == 1) WRITE(lures,"(//,3x,'Slave/Master nodes Data',/)",ERR=9)

 !Initialize empty list
 CALL ini_ndp(head,tail)
 IF (TRIM(actio) == 'NEW' .OR. .NOT.exists('ADD   ')) THEN
   nndp = 0                       !initializes list
 ELSE IF (exists('ADD   ')) THEN
   IF (iwrit == 1) WRITE(lures,                                 &
      "(//,3x,'Kept Slave/Master nodes Data from the previous strategy',/)",ERR=9)
   DO n=1,nndp                       !for previous data
     ALLOCATE (ndp)                  !get memory
     ndp%slave = ndpdat(1,n)         !slave node
     ndp%master = ndpdat(2,n)        !master node
     CALL add_ndp (ndp, head, tail)  !add to list
   END DO
 END IF

 IF (TRIM(actio) == 'NSTRA0') actio='NSTRA1'
 IF (ASSOCIATED (ndpdat) ) DEALLOCATE (ndpdat) !release memory

 IF (iwrit == 1) WRITE(lures,"(/'   Slave Node     Master Node')",ERR=9)

 !Loop to read dat and add them to the list
 DO
   CALL listen('READP0')             !read a line
   IF (exists('ENDNDE')) EXIT        !end of data found => exit
   dummy = exists('SLAVE ',i) .AND. exists('MASTER',j)
   IF( dummy ) THEN
     slave = INT(param(i))
     master= INT(param(j))
   ELSE
     slave = INT(param(1))
     master= INT(param(2))
     !master = getint('MASTER',0,'!-Number of a Master node .........') !get master node
     !slave  = getint('SLAVE ',0,'!Number of a Slave node ...........')  !get slave node
   END IF
   IF( slave /= 0 )THEN
     ALLOCATE (ndp)                    !get memory for a new constraint
     ndp%slave  = slave                !assign
     ndp%master = master               !assign
     IF (iwrit == 1) WRITE(lures, '(2i10)',ERR=9) ndp%slave,ndp%master
     nndp=nndp+1                       !increase number of dependencies
     CALL add_ndp( ndp, head, tail )   !add to list
     IF( iwrit == 1 )WRITE(lures,"(5x,' master node:',i8,'  slave node:',i8)")master,slave
   ELSE IF( TRIM(words(i+1)(1:midn)) == 'SET')THEN
     IF( iwrit == 1 )WRITE(lures,"(5x,' master node:',i8)")master
     sname = get_name(posin=i+1,stype='NSET')      !set name
     CALL nsdb_search (sname, found, ns)   !search if node set exists
     IF (found) THEN              !set found, go ahead
       IF (iwrit == 1) WRITE(lures,"('SLAVE SET ',A)",ERR=9)  TRIM(sname)
       CALL ns_head (ns)          !point to first value in the set
     ELSE          !set not found => ERROR
       WRITE (lures, "(' Set ',A,' not found')",ERR=9) TRIM(sname)
       CALL runend('Kinematic constraints: Set not found ')
     END IF

     DO                         !loop over the nodes in the set
       IF ( end_ns(ns) ) EXIT   !last node processed, EXIT
       slave = get_label(ns)    !get nodal label
       dummy =  slave /= master !check it is not the maste
       IF( dummy )THEN          !if node is slave
         ALLOCATE (ndp)                    !get memory for a new constraint
         ndp%slave  = slave                !assign
         ndp%master = master               !assign
         nndp=nndp+1                       !increase number of dependencies
         CALL add_ndp( ndp, head, tail )   !add to list
       END IF
       CALL ns_next (ns)        !point to next node in the list
     END DO
   ELSE          !set not found => ERROR
     WRITE (lures, "(' Invalid data for SLAVE node/s')",ERR=9)
     CALL runend('Kinematic constraints: invalid data ')
   END IF

 END DO

 !Store data in an array
 ALLOCATE (ndpdat(2,nndp))           !get memory for array
 CALL store_ndp(head)                !transfer from list to array
                                     !and release memory
 CALL ini_ndp(head,tail)             !nullifies pointer

 RETURN
 9 CALL runen2(' error while writing data to the disk')
 END SUBROUTINE readp0
