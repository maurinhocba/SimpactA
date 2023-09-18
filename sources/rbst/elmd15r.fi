 SUBROUTINE elmd15r ( nelem, heade, taile, iwrit )
 !******************************************************************
 !
 !*** READ element DATA for shell element RBST
 !
 !******************************************************************
 IMPLICIT NONE

 INTEGER (kind=4), PARAMETER ::  nnode=2
 INTEGER (kind=4) :: nelem, iwrit
 TYPE (reinf), POINTER :: heade,taile !pointers to first and last elements

 INTEGER (kind=4) :: i,chnode,isec,secty,osec,ip
 INTEGER (kind=4) :: secs = 42      !only acceptable section type
 TYPE (reinf), POINTER :: elem      !auxiliar pointer
 TYPE(section), POINTER :: sec
 LOGICAL :: found


 CALL listen('ELMD15')          !read a card
 IF( .NOT.exists('REINFO') )THEN
   WRITE(lures,"('  No reinforcement for the set of elements ')")
   backs = .TRUE.
   RETURN
 END IF

 ! *** READ the element nodal connections, and the property numbers.

 IF(iwrit == 1) WRITE (lures,"(//6X,'ELEMENT   MAT    NODE NUMBERS ' )",ERR=9999)

 osec = -1
 DO
   ! loop over elements

   CALL listen('ELMD15')          !read a card
   IF (exists('ENDREI')) EXIT     !if key word END_ELEMENT read, exit loop

   nelem = nelem + 1         !increase element counter
   CALL new_reinf(elem)      !reserve memory and initializes data

   elem%numel = INT (param(1))    !element label
   isec = INT (param(2))          !associated material
   !     check if associated material exists
   IF( isec /= osec )THEN
     CALL sect_search (isec,found,sec,ip)
     IF(.NOT.found)CALL runend('ELRBST: Invalid Section Number used for reinforcement')
     secty = sec%secty       !check section type for element
     IF( secty /= secs ) CALL runend('ELRBST: Invalid sect. type for reinforcement')
     IF( sec%mtbas%mtype == 8 ) CALL runend('ELRBST: Invalid mat. type for reinforcement')
     IF( sec%mtbas%matdef(3) > 1 ) CALL runend('ELRBST: Invalid mat. type for reinforcement')
     IF( sec%iprop(1) == 1 ) CALL runend('ELRBST: Invalid mat. type for reinforcement')
     osec = isec
   END IF
   elem%secno = ip        !store section position

   elem%lnods(1:nnode) = INT (param(3:4))  !connectivities
   elem%lnods(3:4) = 0                     !initializes elemetns
   !    echo according to flag
   IF(iwrit == 1) WRITE(lures,"(i15,i5,2i7)",ERR=9999) elem%numel,isec,elem%lnods(1:2)

   !     original nodes are changed
   DO i=1,nnode
     elem%lnods(i) = chnode(elem%lnods(i))
   END DO
   elem%dc = 0d0
   IF( nelem == 1 )THEN
     heade => elem
   ELSE
     taile%next => elem
   END IF
   taile => elem

 END DO

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE elmd15r
