 SUBROUTINE elmda2(nelem,heade,taile,iwrit)
 !******************************************************************
 !
 !*** READ element DATA for truss elements
 !
 !******************************************************************
 IMPLICIT NONE
 INTEGER (kind=4), INTENT(IN) :: iwrit
 INTEGER (kind=4), INTENT(IN OUT) :: nelem
 TYPE(ele02), POINTER :: heade,taile

 TYPE(ele02), POINTER :: elem
 TYPE(section), POINTER :: sec
 INTEGER (kind=4) mtype,i,isec,osec,ip
 LOGICAL :: found,plast
 INTEGER (kind=4) chnode
 REAL (kind=8) :: pinit,area

 ! *** READ the element nodal connections, and the property numbers.

 IF(iwrit == 1)WRITE(lures,"(/4x,'element  mat        node numbers' &
      &                          ,'  Prestress' ,/)",ERR=9999)

 osec = -1
 DO
   CALL listen('ELMD02')          !read a card
   IF (exists('ENDELE')) EXIT     !if key word END_ELEMENT read, exit loop

   nelem = nelem + 1         !increase element counter
   ALLOCATE (elem)           !reserve memory for data

   elem%numel = INT (param(1))    !element label
   isec = INT (param(2))         !associated material
   !     check if associated material exists
   IF( isec /= osec )THEN
     CALL sect_search (isec,found,sec,ip)
     IF(.NOT.found)CALL runend('ELMDA2: Invalid Section Number used')
     mtype = sec%secty       !check section type for element
     IF( mtype /= 1 ) CALL runend('ELMDA2: Invalid sect. type for TRUSS')
     plast =  sec%mtbas%matdef(3) > 1  !plasticiy
     area =  sec%rprop(1)              !transversal area
     osec = isec
   END IF
   elem%matno = ip        !store section position in array
   elem%lnods(1:nnode) = INT (param(3:4))  !connectivities
   pinit = param(5)

   IF(iwrit == 1) WRITE(lures,"(3x,i8,i5,2i10,e12.4)",ERR=9999) elem%numel,     &
                  isec,elem%lnods,pinit

   IF( plast )THEN
     ALLOCATE (elem%gausv(8))
   ELSE
     ALLOCATE (elem%gausv(4))
   END IF

   elem%gausv = 0d0                     !initializes
   elem%gausv(2) = pinit/area           !initializes

   !     original nodes are changed
   DO i=1,nnode
     elem%lnods(i) = chnode(elem%lnods(i))
   END DO

   CALL add_ele02e (elem, heade, taile) !add element to list
 END DO

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE elmda2
