 SUBROUTINE rdrfd0 (iwrit)

 !This routine reads rotation-free dependent nodes (nrfdp)

 USE param_db,ONLY: mnam
 USE ctrl_db,ONLY: ndime
 USE c_input
 USE rfdp_db
 USE npo_db, ONLY : coord,label
 IMPLICIT NONE
 INTEGER :: iwrit

 ! Local
 TYPE (rfdp_nod), POINTER :: rfd_new
 INTEGER :: slave,i,j,k,nnode,master,lnods(6)
 LOGICAL :: sfound,efound,found
 CHARACTER(len=mnam) :: set_name
 INTEGER(kind=4) chnode


 IF (iwrit == 1) WRITE(lures,"(//,3x,'Rotation free Slave/Master Data',/)",ERR=9)

 !Initialize empty list
 IF (nrfdp == 0 .OR. .NOT.exists('ADD   ')) THEN
   CALL ini_rfdp(rfd_head,rfd_tail)                     !initializes list
   planar = exists('PLANAR')
 ELSE IF (exists('ADD   ')) THEN
   IF (iwrit == 1) WRITE(lures,                                 &
      "(//,3x,'Kept Slave/Master nodes Data from the previous strategy',/)",ERR=9)
 END IF

 IF (iwrit == 1) WRITE(lures,"(/'   Slave Node  Master Set  Element   Nodes')",ERR=9)

 nnode = ndime        !number of master nodes
 IF(.NOT.planar) nnode = 3*(ndime-1)

 !Loop to read dat and add them to the list
 DO
   CALL listen('RDRFD0')             !read a line
   IF (exists('ENDRFN')) EXIT        !end of data found => exit
   ALLOCATE( rfd_new)                !get memory for new par
   ALLOCATE( rfd_new%lnods(nnode),rfd_new%rfpar(ndime))  !connectivities and local coordinates
   IF( .NOT.exists('SLAVE',i))i=1    !search for key-word
   rfd_new%slave = INT(param(i))     !get slave node
   ! stop with error if not a valid integer
   IF(rfd_new%slave <= 0) rfd_new%slave= getint('SLAVE ',0,'!Number of a Slave node')
   slave = chnode(rfd_new%slave)                                      !slave internal node
   IF( exists('ELMSET',j))THEN         !MASTER segment is an element in a SET
     set_name = get_name('ELMSET',found,'!-Master Element Set:',stype='ESET',posin=j) !Element set name
     IF( exists('ELEMNT',j))THEN       !master element is given
       master = INT(param(j))          !get master element
       IF( ndime == 2)THEN
         CALL point_ele11e (rfd_new%lnods, planar, master, set_name, sfound, efound)
       ELSE !IF(ndime == 3)THEN
         CALL point_ele13e (rfd_new%lnods, planar, master, set_name, sfound, efound)
       END IF
     ELSE !Search for nearest MASTER element on element set
       master = getint('MNODE',0,' -Number of the Master node')           !get master node
       IF( master /= 0 )master = chnode(master)
       IF( ndime == 2)THEN
         CALL search_ele11 (rfd_new%lnods, set_name, master, sfound, efound, slave, coord, planar)
       ELSE
         CALL search_ele13 (rfd_new%lnods, set_name, master, sfound, efound, slave, coord, planar)
       END IF
     END IF
   ELSE ! MASTER segment connectivities are given
     efound = exists('MNODES',j)     !see if Key-word is given
     IF( .NOT.efound )j = i+1        !point to next paramenter in the list
     DO i=1,nnode
       k = INT(param(j))              !read each node
       rfd_new%lnods(i) = chnode(k)   !convert to internal position
       j = j+1
     END DO
     efound = ANY(rfd_new%lnods > 0 )!minimum verification
     set_name = ''                   !set name is null
     master = 0                      !no master segment
     sfound = .TRUE.                 !set is NOT given because unnecessary
   END IF

   IF( .NOT.sfound )THEN             !set given but not found
     WRITE (lures, "(' Element set not found ')",ERR=9)
     CALL runend('Kinematic constraints: invalid data ')
   ELSE IF( .NOT.efound )THEN        !element not found or nodes not given
     WRITE (lures, "(' Element not found in the set ')",ERR=9)
     CALL runend('Kinematic constraints: invalid data ')
   ELSE                              !print data read
     IF (iwrit == 1) THEN
       IF( LEN_TRIM(set_name) == 0 )THEN   !set not given, print nodes
         lnods = 0
         DO i=1,nnode
           IF( rfd_new%lnods(i) > 0) lnods(i) = label(rfd_new%lnods(i))
         END DO
         WRITE(lures,"(i10,20x,6i6)",ERR=9) lnods(1:nnode)
       ELSE                                !set given, print set and element
         WRITE(lures,"(i10,3x,a,i7)",ERR=9) rfd_new%slave,TRIM(set_name),master
       END IF
     END IF
     rfd_new%slave = slave                          !keep slave internal node
     CALL add_rfdp( rfd_new, rfd_head, rfd_tail )   !add to list
   END IF

 END DO

 RETURN
 9 CALL runen2(' error while writing data to the disk')
 END SUBROUTINE rdrfd0
