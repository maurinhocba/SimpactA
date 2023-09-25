 SUBROUTINE elmd11(nelem,nstre,ntype,ngaus,heade,taile,iwrit)
 !******************************************************************
 !
 !*** READ element DATA for 2-D beam/shell elements
 !
 !******************************************************************
 IMPLICIT NONE

 ! dummy arguments
 INTEGER (kind=4), INTENT(IN) :: nstre,ntype,ngaus,iwrit
 INTEGER (kind=4), INTENT(IN OUT) :: nelem
 TYPE (ele11), POINTER :: heade,taile !pointers to first and last elements

 ! local variables
 INTEGER (kind=4) i,isec,mtype,osec,nvare,nlayr,ip
 TYPE(ele11), POINTER :: elem
 TYPE (section), POINTER :: sec      !auxiliar pointer
 LOGICAL found,plast
 INTEGER(kind=4) chnode
 INTEGER (kind=4) :: secs = 43

! *** READ the element nodal connections, and the property numbers.
 osec = -1
 IF(iwrit == 1) WRITE (lures,"(/' ELEMENT  SECC      NODE NUMBERS ')",ERR=9999)
 DO
   ! loop over elements

   CALL listen('ELMD11')          !read a card
   IF (exists('ENDELE')) EXIT     !if key word END_ELEMENT read, exit loop

   nelem = nelem + 1         !increase element counter
   CALL new_ele11e(elem)     !reserve memory for data
   IF( ntype == 3 ) ALLOCATE ( elem%ro(2,ngaus) )
   elem%lambd = 1d0
   ! initializes

   elem%numel = INT (param(1))    !element label
   isec = INT (param(2))         !associated material
   !     check if associated material exists
   IF( isec /= osec )THEN
     CALL sect_search (isec,found,sec,ip)
     IF(.NOT.found)CALL runend('SHREV:  Invalid Section Number used')
     mtype = sec%secty       !check section type for element
     IF( mtype /= secs ) CALL runend('SHREV : Invalid sect. type for SHREV')

     nlayr = sec%iprop(1)     !number of layers
     plast = sec%mtbas%matdef(3) > 1   !plasticiy
     nvare = nstre/2 + 1       !number of internal variables
     IF( ntype == 1 .AND. sec%mtbas%matdef(5) > 1 )nvare = nvare+1
     osec = isec
   END IF
   elem%matno = ip        !store section position in array
   elem%lnods(2:3) = INT (param(3:4))  !connectivities
   elem%lnods(1) = 0      !initializes
   elem%lnods(4) = 0
   !    echo according to flag
   IF(iwrit == 1) WRITE(lures,"(i8,i5,3i7)",ERR=9999) elem%numel,isec,elem%lnods(2:3)
   !     original nodes are changed
   DO i=2,3
     elem%lnods(i) = chnode(elem%lnods(i))
   END DO

   IF( plast )THEN   !reserve space for Gauss variables
     ALLOCATE( elem%gausv(nvare,nlayr,ngaus) )
     elem%gausv = 0d0
   ELSE
     NULLIFY( elem%gausv )
   END IF
   CALL add_ele11e (elem, heade, taile) !add element to list
 END DO

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE elmd11
