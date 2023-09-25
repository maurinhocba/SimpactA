 SUBROUTINE elmd10(nnode,ndime,nelem,lnods,ntype,&
                   matno,iwrit,rbnod,nmast,flag,heat)
 !***************************************************************** *
 !
 !***  READ element DATA for Rigid bodies
 !
 !***************************************************************** *
 IMPLICIT NONE
 INTEGER (kind=4), INTENT(IN) :: nnode, & ! number of nodes per element
                                 ndime, & ! problem dimension
                                 nelem, & ! number of elements
                                 iwrit, & ! write flag
                                 ntype, & ! problem type for 2-D case
                                 rbnod    ! label of master node
 INTEGER (kind=4), INTENT(OUT) :: lnods(:,:), & !element connectivities
                                  nmast         !Associated Node to master
 INTEGER (kind=4), POINTER  ::  matno(:)      !associated materials
 LOGICAL (kind=4), INTENT(IN) ::  flag,    &  !associated material
                                  heat        !coupled problem

 LOGICAL found,diff
 INTEGER (kind=4) i,numel,j,isec,stype,osec,ip
 INTEGER (kind=4) chnode
 TYPE(section), POINTER :: sec

 !     *** READ the element nodal connections, and the property numbers.

 IF(iwrit == 1) WRITE (lures,"(//5x,' Elmnt Mat        Node Numbers')",ERR=9999)

 osec = -1
 diff = .FALSE.
 DO i=1,nelem
   CALL rdfrin('ELMD10',intrd,nintg,maxni)
   numel = intrd(1)                         !element label (not stored)
   lnods(1:nnode,i) = intrd(3:nnode+2)  !element nodes

   IF(iwrit == 1)WRITE(lures,"(4x,i7,1i5,10i8)",ERR=9999) numel,intrd(2),    &
                               lnods(1:nnode,i)
   isec = intrd(2)                  !associated material
   diff = diff .OR. (numel /= i)
   IF( flag )THEN
     IF( isec /= osec )THEN
       CALL sect_search (isec,found,sec,ip)
       IF(.NOT.found) CALL runend('RIGID : Invalid Section Number used')
       stype = sec%secty        !check section type for element
       IF( ntype <=3 )THEN !solid
         IF(stype /= 0) CALL runend('RIGID: Invalid sec. type for RIGID')
         IF( heat )THEN
           IF( sec%mtbas%mtype /= 4 ) &
             CALL runend('RIGID: Invalid mat. type for RIGID-HEAT')
         END IF
       ELSE                                 !shells
         IF( ndime == 2 )THEN      !2-D shells
           IF( stype /= 43 ) CALL runend('RIGID: Invalid sec. type for RIGID')
         ELSE                      !3-D shelss
           IF( stype /= 11 .AND. stype /= 12 .AND. stype /= 14 ) &
           CALL runend('RIGID: Invalid sec. type for RIGID')
         END IF
       END IF
       osec = isec
     END IF
     matno(i) = ip     !associate Rigid Bodies with a material
   END IF
 END DO
 IF( diff )WRITE(lures,"(' WARNING: elements not read in order,', &
           &              ' original element labels modified')")
 !     original nodes are changed
 IF( rbnod > 0 )THEN
   nmast = chnode(rbnod)      !stores Master node
 ELSE
   nmast = 0
 END IF
 DO j=1,nelem
   DO i=1,nnode
     IF( lnods(i,j) == 0 )THEN    !some other options may exist
       lnods(i,j) = lnods(i-1,j)
     ELSE
       lnods(i,j) = chnode(lnods(i,j))
     END IF
   END DO
 END DO

 CALL listen('ELMD10')
 IF(.NOT.exists('ENDELE'))CALL runend('ELMD10: END_ELEMENT expected')

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE elmd10
