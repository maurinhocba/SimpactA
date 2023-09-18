 SUBROUTINE elmda1(nnode,ndime,neulr,nelem,heade,taile,iwrit)
 !******************************************************************
 !
 !*** READ element DATA for truss elements
 !
 !******************************************************************
 IMPLICIT NONE
 INTEGER (kind=4), INTENT(IN) :: nnode,ndime,neulr,iwrit
 INTEGER (kind=4), INTENT(IN OUT) :: nelem
 TYPE(ele01), POINTER :: heade,taile

 TYPE(ele01), POINTER :: elem
 TYPE(section), POINTER :: sec
 INTEGER (kind=4) stype,i,isec,osec,nvare,ip
 LOGICAL :: found,rot,same
 INTEGER (kind=4) chnode

 ! *** READ the element nodal connections, and the property numbers.

 IF(iwrit == 1)WRITE(lures,"(/4x,'element  mat        node numbers'/)",ERR=9999)

 osec = -1
 DO
   CALL listen('ELMD01')          !read a card
   IF (exists('ENDELE')) EXIT     !if key word END_ELEMENT read, exit loop

   nelem = nelem + 1         !increase element counter
   ALLOCATE (elem)           !reserve memory for data
   ALLOCATE (elem%lnods(nnode))  !reserve memory for nodes

   elem%numel = INT (param(1))    !element label
   isec = INT (param(2))         !associated material
   !     check if associated material exists
   IF( isec /= osec )THEN
     CALL sect_search (isec,found,sec,ip)
     IF(.NOT.found)CALL runend('ELMDA1: Invalid Section Number used')
     stype = sec%secty       !check section type for element
     IF( stype /= 2 ) CALL runend('ELMDA1: Invalid sect. type for SPOT ')
     osec = isec
     IF( nnode == 2 )THEN
       rot = ( sec%rprop(2) > 0d0 .OR. sec%rprop(4) > 0d0 .OR.   &
               sec%iprop(2) > 0   .OR. sec%iprop(2) > 0 ) .AND. neulr > 0
     ELSE
       rot = .FALSE.
     END IF
   END IF
   elem%matno = ip        !store section Label

   IF( nnode == 2) THEN
     elem%lnods(1:2) = INT (param(3:4))  !connectivities
     same =  (elem%lnods(1) == elem%lnods(2)) .OR. (elem%lnods(2) == 0)
     IF( same ) THEN
       IF(iwrit == 1) WRITE(lures,"(3x,i8,i5,6i7)",ERR=9999) elem%numel,     &
                      isec,elem%lnods(1)
       elem%lnods(2) = 0

       IF( ndime == 2 ) THEN
         IF( rot )THEN
           nvare = 5
         ELSE      ! no Euler angles
           nvare = 3
         END IF
       ELSE !ndime ==3
         IF( rot )THEN
           nvare = 14
         ELSE      ! no Euler angles
           nvare = 4
         END IF
       END IF
     ELSE    !different nodes
       IF(iwrit == 1) WRITE(lures,"(3x,i8,i5,6i7)",ERR=9999) elem%numel,     &
                      isec,elem%lnods(1:2)
       IF( rot )THEN
         nvare = 2
       ELSE      ! no Euler angles
         nvare = 1
       END IF
     END IF
   ELSE
     nvare = 22
     SELECT CASE (nnpar)
     CASE(3)
       elem%lnods = (/ INT (param(3)), 0,0,0,0,0 /)  !auxiliar node
       IF(iwrit == 1) WRITE(lures,"(3x,i8,i5,i10)",ERR=9999) elem%numel,isec,elem%lnods(1)
     CASE(8)
       elem%lnods(1:nnode) = INT (param(3:8))  !connectivities
       IF(iwrit == 1) WRITE(lures,"(3x,i8,i5,6i7)",ERR=9999) elem%numel,isec,elem%lnods
     CASE(9)
       elem%lnods(1:nnode) = INT (param(3:8))  !connectivities
       i = INT(param(9))
       IF(iwrit == 1) WRITE(lures,"(3x,i8,i5,6i7,i10)",ERR=9999) elem%numel,isec,elem%lnods,i
       i = chnode(i)
       elem%numel = - i
     END SELECT
   END IF
   ALLOCATE (elem%gausv(nvare))

   elem%gausv = 0d0                     !initializes
   elem%rot   = rot

   !     original nodes are changed
   DO i=1,nnode
     IF( elem%lnods(i) > 0 ) elem%lnods(i) = chnode(elem%lnods(i))
   END DO

   CALL add_ele01e (elem, heade, taile) !add element to list
 END DO

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE elmda1
