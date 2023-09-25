 SUBROUTINE rearb0 (nnode,lnods,iwrit,rbnod,nmast)

 !This routine reads rigid bodies defined from particles
 USE param_db, ONLY : mnam
 USE c_input, ONLY : lures,listen,exists,get_name,nnpar,param
 USE nsets_db
 IMPLICIT NONE
 INTEGER (kind=4), INTENT(IN) :: iwrit, &   !flag to print
                                 rbnod      !master node label
 INTEGER(kind=4), INTENT(OUT) :: nnode, &   !number of slave nodes
                                 nmast      !master node
 INTEGER(kind=4), POINTER :: lnods(:,:)
 ! Local
 TYPE (node), POINTER :: nhead,ntail,nn
 TYPE (nset), POINTER :: ns
 INTEGER :: n,nns,kf,nl
 LOGICAL :: set, found
 CHARACTER(len=mnam) :: sname
 INTEGER (kind=4) chnode

 IF (iwrit == 1) WRITE(lures,"(//,3x,'Rigid bodies data',/)",ERR=9999)

   nns = 0        !initializes number of slaves nodes
   DO             !loop to read slave nodes
     CALL listen('REARB0')             !read a line
     IF (exists('SET   ',n))THEN       !a set of nodes
       sname = get_name(posin=n,stype='NSET')       !set name
       CALL nsdb_search (sname, found, ns)
       IF (found) THEN              !if set found get node labels
         kf = get_length (ns)       !number of nodes in the set
         set = .TRUE.
         CALL ns_head (ns)          !go to top of the listt
       ELSE                         !error in input data
         WRITE (lures,"(' Set ',A,' not found')",ERR=9999) TRIM(sname)
         CALL runend('Intime: Set not found ')
       END IF
     ELSE IF(exists('ENDELE'))THEN     !last line read
       ALLOCATE( lnods(nns,1) )        !get memory for slave nodes
       nnode = nns                     !number of slave nodes
       nn => nhead                     !point to first node of the list
       DO n=1,nns                      !loop over all the nodes
         lnods(n,1) = chnode(nn%label) !assign number to slave node
         nhead => nn                   !point head to the node
         nn => nn%next                 !point to next node in the list
         DEALLOCATE(nhead)             !release transfered node
       END DO
       nmast = chnode(rbnod)           !assign number to master node
       EXIT                            !exit reading loop

     ELSE                              !a list of nodes
       kf = nnpar                      !number of nodes read
       set = .FALSE.                   !
     END IF

     DO n = 1,kf                      !for each node
       IF (set) THEN
         nl = get_label(ns)              !get node label
         CALL ns_next (ns)               !go to next node
       ELSE
         nl = INT(param(n))              !read label from list
       END IF
       IF( nl /= rbnod )THEN             !not the master node
         nns = nns+1                     !increment counter
         ALLOCATE( nn )                  !get memory for node
         nn%label = nl                   !keep node label
         IF (iwrit == 1) WRITE(lures, '(8i9)',ERR=9999) nl
         IF( nns == 1 )THEN              !first node of the list
           nhead => nn                   !initializes list
           ntail => nn
         ELSE                            !any other node
           ntail%next => nn              !add to list
           ntail => ntail%next
         END IF
         NULLIFY (nn)
       END IF
     END DO

   END DO

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE rearb0
